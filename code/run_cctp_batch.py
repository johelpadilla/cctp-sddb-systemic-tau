#!/usr/bin/env python3
"""
run_cctp_batch.py

Surgical batch orchestrator for CCTP v0.1 full pipeline on multiple SDDB records.

Reuses without modification:
  - analyze_cctp_pilot.py
  - run_cctp_surrogates.py
  - run_recd_on_rr.py
  - run_recd_weighted_on_rr.py
  - extract_rr.py (if clean npz missing)

For each record:
  1. Ensure rr_{rec}_clean.npz (calls extract if needed)
  2. Run full analysis stack (idempotent-ish: skip if outputs exist unless --force)
  3. Collect metrics

Outputs:
  - results/cctp_batch_summary.csv   (main consolidated table)
  - figures/batch/                   (comparative plots)
  - Per-record results/ + figures/   (already produced by the 4 scripts)

Usage (surgical):
  python3 code/run_cctp_batch.py --records 30,31,35
  python3 code/run_cctp_batch.py --list selected_records.txt --force
  python3 code/run_cctp_batch.py   # defaults to selected_records.txt
"""

import argparse
import os
import sys
import subprocess
import json
import csv
from pathlib import Path
from datetime import datetime
import numpy as np

try:
    import pandas as pd
    import matplotlib.pyplot as plt
    HAS_PANDAS_MPL = True
except Exception:
    HAS_PANDAS_MPL = False
    pd = None
    plt = None

BASE = Path(__file__).resolve().parent.parent
CODE_DIR = BASE / "code"
DATA_DIR = BASE / "data"
RES_DIR = BASE / "results"
FIG_BATCH = BASE / "figures" / "batch"
RES_DIR.mkdir(parents=True, exist_ok=True)
FIG_BATCH.mkdir(parents=True, exist_ok=True)

SCRIPTS = [
    ("analyze", "analyze_cctp_pilot.py"),
    ("surrogates", "run_cctp_surrogates.py"),
    ("recd", "run_recd_on_rr.py"),
    ("weighted", "run_recd_weighted_on_rr.py"),
]

def load_records(arg_records: str, arg_list: str) -> list:
    if arg_records:
        return [r.strip() for r in arg_records.split(",") if r.strip()]
    if arg_list:
        p = Path(arg_list)
    else:
        p = BASE / "selected_records.txt"
    if p.exists():
        recs = []
        for line in p.read_text().splitlines():
            line = line.strip()
            if line and not line.startswith("#"):
                recs.append(line)
        return recs
    return ["30", "35"]

def clean_npz_exists(rec: str) -> bool:
    return (DATA_DIR / f"rr_{rec}_clean.npz").exists()

def ensure_clean_npz(rec: str):
    if clean_npz_exists(rec):
        return
    print(f"  [{rec}] clean npz missing → running extraction...")
    cmd = [sys.executable, str(CODE_DIR / "extract_rr.py"), "--record", rec]
    subprocess.run(cmd, check=True, cwd=BASE)

def outputs_exist(rec: str) -> dict:
    """Check which pipeline stages have outputs."""
    status = {}
    # analyze summary
    status["analyze"] = (RES_DIR / f"cctp_pilot_summary_{rec}.json").exists()
    # surrogates
    status["surrogates"] = (RES_DIR / f"surrogate_cctp_{rec}.json").exists()
    # recd levels
    status["recd"] = (RES_DIR / f"recd_rr_{rec}.json").exists()
    # weighted
    status["weighted"] = (RES_DIR / f"recd_weighted_rr_{rec}.json").exists()
    # figures dir has the late ones
    status["figures"] = (BASE / "figures" / rec / "19_recd_weighted_box.png").exists()
    return status

def run_stage(rec: str, stage: str, script_name: str, force: bool, dry_run: bool, extra_args=None):
    out_status = outputs_exist(rec)
    if out_status.get(stage, False) and not force:
        print(f"  [{rec}] {stage}: outputs exist, skipping (use --force to rerun)")
        return True
    cmd = [sys.executable, str(CODE_DIR / script_name), "--record", rec]
    if extra_args:
        cmd.extend(extra_args)
    print(f"  [{rec}] {'[DRY] would run' if dry_run else 'running'} {stage} → {script_name} {' '.join(extra_args or [])}")
    if dry_run:
        return True
    try:
        res = subprocess.run(cmd, cwd=BASE, capture_output=True, text=True, timeout=600)
        if res.returncode != 0:
            print(f"    ERROR (rc={res.returncode}): {res.stderr[-500:] if res.stderr else ''}")
            return False
        print(f"    OK")
        return True
    except subprocess.TimeoutExpired:
        print(f"    TIMEOUT for {rec} {stage}")
        return False
    except Exception as e:
        print(f"    EXCEPTION: {e}")
        return False

def load_key_metrics(rec: str) -> dict:
    """Aggregate key numbers from the authoritative JSONs (structures from actual files) + quality from npz."""
    m = {"record": rec}
    # Pilot summary
    pfile = RES_DIR / f"cctp_pilot_summary_{rec}.json"
    if pfile.exists():
        with open(pfile) as f:
            p = json.load(f)
        m["duration_h"] = p.get("total_hours") or p.get("duration_h")
        m["event_hr"] = p.get("event_hr")
        mets = p.get("metrics", {})
        ts = mets.get("tau_s", {})
        m["delta_tau"] = ts.get("delta")
        m["p_tau"] = ts.get("p_welch")
        m["delta_var"] = mets.get("var", {}).get("delta")
        m["delta_ar1"] = mets.get("ar1", {}).get("delta")
    # surrogates
    sfile = RES_DIR / f"surrogate_cctp_{rec}.json"
    if sfile.exists():
        with open(sfile) as f:
            s = json.load(f)
        # support both legacy p_value and current p_direction_specific / p_two_sided_abs
        m["p_tau_surrogate"] = s.get("p_value") or s.get("p_direction_specific") or s.get("p_two_sided_abs")
        m["n_surrogates"] = s.get("n_surr") or s.get("n_surrogates")
    # RECD levels (unweighted)
    rfile = RES_DIR / f"recd_rr_{rec}.json"
    if rfile.exists():
        with open(rfile) as f:
            r = json.load(f)
        ex = r.get("mean_excess3") or {}
        m["mean_excess3_basal"] = ex.get("basal")
        m["mean_excess3_approach"] = ex.get("approach")
        m["delta_excess3_unw"] = ex.get("delta")
        m["p_excess3_unw"] = ex.get("p_welch")
        hl = r.get("high_level3_rate") or {}
        m["high_level3_basal"] = hl.get("basal")
        m["high_level3_approach"] = hl.get("approach")
        m["delta_high_level3"] = hl.get("delta")
    # weighted (preferred for primary columns when available - full pipeline)
    wfile = RES_DIR / f"recd_weighted_rr_{rec}.json"
    if wfile.exists():
        with open(wfile) as f:
            w = json.load(f)
        exw = w.get("mean_excess3") or {}
        m["delta_excess3_w"] = exw.get("delta")
        m["p_excess3_w"] = exw.get("p_welch")
        fc = w.get("frac_contrib3") or {}
        m["delta_frac_contrib3"] = fc.get("delta")
        m["p_frac3"] = fc.get("p_welch")
        lam = w.get("mean_lambda") or {}
        m["mean_lambda"] = lam.get("approach") if isinstance(lam, dict) else lam
        hlw = w.get("high_level3_rate") or {}
        m["high_level3_basal_w"] = hlw.get("basal")
        m["high_level3_approach_w"] = hlw.get("approach")
        m["delta_high_level3_w"] = hlw.get("delta")
        # promote weighted to primary view
        m["delta_excess3"] = exw.get("delta")
        m["p_excess3"] = exw.get("p_welch")
    else:
        # fallback to unweighted for primary columns
        m["delta_excess3"] = m.get("delta_excess3_unw")
        m["p_excess3"] = m.get("p_excess3_unw")

    # === Quality + event timing from clean npz (authoritative for event_hr) ===
    npz_path = DATA_DIR / f"rr_{rec}_clean.npz"
    if npz_path.exists():
        try:
            npz = np.load(npz_path, allow_pickle=True)
            m["n_beats"] = int(npz["n_beats"]) if "n_beats" in npz else None
            m["interp_frac"] = float(npz["interp_frac"]) if "interp_frac" in npz else None
            m["pacing_detected"] = bool(npz["pacing_detected"]) if "pacing_detected" in npz else None
            m["known_pacing_type"] = str(npz["known_pacing_type"]) if "known_pacing_type" in npz else "none"
            m["cv_rr"] = float(npz["cv_rr"]) if "cv_rr" in npz else None
            # Always resolve event_hr from RR series + vfon (pilot JSON often omits it)
            try:
                from cctp_metrics_core import resolve_event_timing_from_npz

                timing = resolve_event_timing_from_npz(rec, str(npz_path))
                m["event_hr"] = timing["event_hr"]
                if m.get("duration_h") is None:
                    m["duration_h"] = timing["duration_h"]
                m["event_type"] = timing["event_type"]
            except Exception as te:
                print(f"  [{rec}] WARNING: could not resolve event_hr from npz: {te}")
        except Exception as e:
            print(f"  [{rec}] WARNING: could not load quality from npz: {e}")

    # Explicit flag for whether full weighted RECD was completed for this record
    m["has_weighted"] = (RES_DIR / f"recd_weighted_rr_{rec}.json").exists()
    return m

def write_batch_summary(metrics_list: list):
    if not metrics_list:
        return
    keys = set()
    for m in metrics_list:
        keys.update(m.keys())
    # order preferred (richer with quality + re-cal high_level3)
    pref = ["record", "duration_h", "event_hr", "delta_tau", "p_tau", "p_tau_surrogate",
            "delta_var", "delta_ar1", "delta_excess3", "p_excess3", "delta_excess3_w",
            "p_excess3_w", "delta_frac_contrib3", "p_frac3", "mean_lambda", "n_surrogates",
            "n_beats", "interp_frac", "pacing_detected", "known_pacing_type", "cv_rr",
            "has_weighted", "delta_high_level3", "high_level3_basal", "high_level3_approach",
            "delta_high_level3_w"]
    fieldnames = [k for k in pref if k in keys] + sorted(k for k in keys if k not in pref)

    out_csv = RES_DIR / "cctp_batch_summary.csv"
    with open(out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for m in metrics_list:
            row = {k: m.get(k) for k in fieldnames}
            w.writerow(row)
    print(f"\nWrote consolidated summary: {out_csv}")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--records", default="", help="Comma list e.g. 30,31,35")
    ap.add_argument("--list", default="", help="selected_records.txt path")
    ap.add_argument("--force", action="store_true", help="Re-run stages even if outputs exist")
    ap.add_argument("--skip-extract", action="store_true")
    ap.add_argument("--dry-run", action="store_true", help="Only check what would run, do not execute heavy stages")
    # Re-calibration light overrides (passed to recd/weighted stages)
    ap.add_argument("--theta3", type=float, default=None, help="Re-cal theta3 e.g. 0.08")
    ap.add_argument("--high-thresh", type=float, default=None, help="Re-cal high_thresh e.g. 0.65")
    ap.add_argument("--lambda-theta", type=float, default=None, help="Re-cal lambda theta_chaos e.g. 0.09")
    ap.add_argument("--lambda-relative", action="store_true", help="Use record-relative lambda scaling")
    args = ap.parse_args()

    recs = load_records(args.records, args.list)
    print(f"CCTP batch for records: {recs}")
    print(f"Force={args.force}  Dry-run={args.dry_run}")

    metrics = []

    for rec in recs:
        print(f"\n{'='*60}\nRECORD {rec}\n{'='*60}")
        if not args.skip_extract:
            ensure_clean_npz(rec)

        # Run the four stages
        # theta3/high go to both recd and weighted
        common_recal = []
        if args.theta3 is not None:
            common_recal.extend(["--theta3", str(args.theta3)])
        if args.high_thresh is not None:
            common_recal.extend(["--high-thresh", str(args.high_thresh)])

        # lambda args only for weighted
        weighted_only = []
        if args.lambda_theta is not None:
            weighted_only.extend(["--lambda-theta", str(args.lambda_theta)])
        if args.lambda_relative:
            weighted_only.append("--lambda-relative")

        all_ok = True
        for stage_name, script in SCRIPTS:
            if stage_name == "recd":
                stage_extra = common_recal if common_recal else None
            elif stage_name == "weighted":
                stage_extra = (common_recal + weighted_only) if (common_recal or weighted_only) else None
            else:
                stage_extra = None
            ok = run_stage(rec, stage_name, script, args.force, args.dry_run, extra_args=stage_extra)
            all_ok = all_ok and ok

        if all_ok:
            m = load_key_metrics(rec)
            metrics.append(m)
            print(f"  [{rec}] metrics collected: delta_tau={m.get('delta_tau')}, delta_excess3={m.get('delta_excess3')}")
        else:
            print(f"  [{rec}] WARNING: some stages failed or skipped, partial metrics")

    write_batch_summary(metrics)

    # Minimal batch comparative note
    print("\nBatch complete.")
    print("Next: inspect results/cctp_batch_summary.csv + figures/*/ + figures/batch/ (add comparative plots here if needed).")
    print("To generate simple cross-record comparison, you can extend this script or use pandas/plot later.")

    # Generate richer batch summary report + basic comparative figures
    if metrics:
        print_quality_report(metrics)
        if HAS_PANDAS_MPL:
            generate_batch_figures(metrics)
        else:
            print("  (pandas/matplotlib not available for batch figures; install if needed)")

def print_quality_report(metrics_list):
    print("\n" + "="*60)
    print("BATCH QUALITY REPORT (per record)")
    print("="*60)
    for m in metrics_list:
        rec = m.get("record")
        nb = m.get("n_beats")
        ifrac = m.get("interp_frac")
        pdet = m.get("pacing_detected")
        ptype = m.get("known_pacing_type", "none")
        cv = m.get("cv_rr")
        ifrac_str = f"{ifrac*100:.1f}%" if ifrac is not None else "n/a"
        cv_str = f"{cv:.4f}" if cv is not None else "n/a"
        print(f"  {rec}: n_beats={nb}  interp={ifrac_str}  pacing_detected={pdet} ({ptype})  cv={cv_str}")
    print("Use this for exclusion decisions and methods section.")

def generate_batch_figures(metrics_list):
    """Minimal but useful comparative plots for preprint (saved to figures/batch/)."""
    import pandas as pd
    import matplotlib.pyplot as plt
    df = pd.DataFrame(metrics_list)
    FIG_BATCH.mkdir(parents=True, exist_ok=True)

    # 1. Delta tau and delta_excess3 box-ish (actually bar per record)
    fig, ax = plt.subplots(figsize=(10,4))
    recs = df["record"].astype(str).tolist()
    x = np.arange(len(recs))
    width=0.35
    if "delta_tau" in df:
        ax.bar(x - width/2, df["delta_tau"].fillna(0), width, label="Δτ_s")
    if "delta_excess3" in df:
        ax.bar(x + width/2, df["delta_excess3"].fillna(0), width, label="Δexcess3")
    ax.set_xticks(x)
    ax.set_xticklabels(recs)
    ax.axhline(0, color="k", lw=0.5)
    ax.set_ylabel("Delta (approach - basal)")
    ax.set_title("CCTP Batch: Δτ_s and Δexcess3 across records")
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIG_BATCH / "batch_delta_tau_excess3.png", dpi=150)
    plt.close(fig)

    # 2. Significance heatmap-ish (simple p-values)
    fig, ax = plt.subplots(figsize=(8,3))
    pcols = []
    labels = []
    for col, lab in [("p_tau", "p_Δτ"), ("p_excess3", "p_excess3")]:
        if col in df:
            pcols.append( (-np.log10( df[col].replace(0, 1e-300).astype(float) + 1e-300 )) )
            labels.append(lab)
    if pcols:
        mat = np.array(pcols)
        im = ax.imshow(mat, aspect="auto", cmap="viridis")
        ax.set_yticks(range(len(labels)))
        ax.set_yticklabels(labels)
        ax.set_xticks(range(len(recs)))
        ax.set_xticklabels(recs)
        ax.set_title("-log10(p) for key deltas (higher = stronger)")
        fig.colorbar(im, ax=ax)
        fig.tight_layout()
        fig.savefig(FIG_BATCH / "batch_significance.png", dpi=150)
        plt.close(fig)

    # 3. Quality: interp % and pacing flag
    if "interp_frac" in df:
        fig, ax = plt.subplots(figsize=(9,3))
        ax.bar(recs, (df["interp_frac"].fillna(0)*100))
        ax.set_ylabel("% interpolated RR")
        ax.set_title("Batch data quality: interpolation fraction")
        for i, (r, p) in enumerate(zip(recs, df.get("pacing_detected", [False]*len(recs)))):
            if p:
                ax.text(i, (df["interp_frac"].iloc[i] or 0)*100 + 1, "P", ha="center", color="red")
        fig.tight_layout()
        fig.savefig(FIG_BATCH / "batch_quality_interp.png", dpi=150)
        plt.close(fig)

    print(f"  Batch figures written to {FIG_BATCH} (delta, significance, quality)")

if __name__ == "__main__":
    main()
