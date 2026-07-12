#!/usr/bin/env python3
"""
External Validation Phase 1 — frozen pipeline on VFDB / CUDB + NSRDB controls.

Frozen discovery parameters (Jul-12 2026) — DO NOT RETUNE:
  θ₃=0.08, high-threshold=0.65, W_TAU=101, stride=5,
  detector abs-z ≥ 2.0 sustained ≥ 3 consecutive windows.

Outputs under results/:
  external_phase1_inventory.csv
  external_phase1_per_record.csv
  external_phase1_sensitivity.json
  external_phase1_controls.csv
  external_phase1_far.json
  external_phase1_summary.json

Does NOT claim clinical deployability.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _bootstrap import import_systemictau_core
from cctp_metrics_core import (
    FROZEN_HIGH_THRESHOLD,
    FROZEN_MIN_CONSECUTIVE,
    FROZEN_THETA3,
    FROZEN_Z_THRESHOLD,
    STRIDE,
    W_EWS,
    W_TAU,
    build_bivariate_proxy,
    count_alarm_episodes,
    detect_lead_time,
    detector_performance,
    false_alarm_rate,
    regime_delta,
    rolling_ar1,
    rolling_var,
    short_db_windows,
    sign_concordance,
)
from recd_ordinal_levels import compute_phi3, generate_multivariate_symbols

BASE = Path(__file__).resolve().parent.parent
DATA = BASE / "data"
RR_EXT = DATA / "rr_external"
RES = BASE / "results"
RES.mkdir(parents=True, exist_ok=True)

# Inclusion floors (hours)
MIN_PRE_EVENT_PREFERRED_H = 1.0  # 60 min preferred
MIN_PRE_EVENT_ABSOLUTE_H = 0.25  # 15 min absolute (short-DB stratum)
MAX_INTERP_FRAC = 0.15


def compute_tau_series(rr: np.ndarray) -> np.ndarray:
    compute_taus, _, has = import_systemictau_core()
    X = build_bivariate_proxy(rr)
    if has and compute_taus is not None:
        taus_global, _ = compute_taus(X, window_size=W_TAU, stride=STRIDE)
        return np.asarray(taus_global, dtype=float)
    n = len(rr)
    out = np.full(n, np.nan)
    for i in range(W_TAU - 1, n, STRIDE):
        win = X[i - W_TAU + 1 : i + 1]
        r0 = win[:, 0].argsort().argsort().astype(float)
        r1 = win[:, 1].argsort().argsort().astype(float)
        if np.std(r0) > 0 and np.std(r1) > 0:
            out[i] = np.corrcoef(r0, r1)[0, 1]
    return out


def compute_excess3_series(rr: np.ndarray, theta3: float = FROZEN_THETA3) -> Tuple[int, np.ndarray]:
    X = build_bivariate_proxy(rr)
    S = generate_multivariate_symbols(X, m=3, delay=1)
    offset = (3 - 1) * 1
    phi3, excess3 = compute_phi3(S, window=W_TAU, theta=theta3, stride=STRIDE)
    return offset, np.asarray(excess3, dtype=float)


def load_npz(path: Path) -> dict:
    d = np.load(path, allow_pickle=True)
    out = {}
    for k in d.files:
        v = d[k]
        if isinstance(v, np.ndarray) and v.shape == () and v.dtype == object:
            out[k] = str(v.item())
        elif isinstance(v, np.ndarray) and v.dtype.kind in ("U", "S") and v.shape == ():
            out[k] = str(v.item())
        else:
            out[k] = v
    return out


def inventory_row_from_npz(path: Path) -> dict:
    d = load_npz(path)
    db = str(d.get("source_db", "unknown"))
    rec = str(d.get("record_id", path.stem))
    total_h = float(d.get("total_hours", 0.0))
    pre_h = float(d.get("pre_event_hours", 0.0))
    interp = float(d.get("interp_frac", 0.0))
    n_beats = int(d.get("n_beats", 0))
    is_control = bool(d.get("is_control", False)) or db == "nsrdb"
    event_label = str(d.get("event_label", ""))

    reasons = []
    include = True
    if n_beats < 200:
        include = False
        reasons.append("n_beats<200")
    if interp > MAX_INTERP_FRAC:
        include = False
        reasons.append(f"interp_frac>{MAX_INTERP_FRAC}")
    if not is_control:
        if pre_h < MIN_PRE_EVENT_ABSOLUTE_H:
            include = False
            reasons.append(f"pre_event<{MIN_PRE_EVENT_ABSOLUTE_H}h")
        # event time uncertain / missing
        if not np.isfinite(float(d.get("vfon_sec", np.nan))):
            include = False
            reasons.append("no_event_time")
    if is_control and total_h < 1.0:
        include = False
        reasons.append("control_total_h<1")

    if is_control:
        stratum = "control_nsrdb" if db == "nsrdb" else "control"
    elif pre_h >= 6.0:
        stratum = "holter_ge6h"
    elif pre_h >= 1.0:
        stratum = "short_ge60min"
    elif pre_h >= MIN_PRE_EVENT_ABSOLUTE_H:
        stratum = "short_15_60min"
    else:
        stratum = "too_short"

    return {
        "source_db": db,
        "record_id": rec,
        "npz_path": str(path.relative_to(BASE)) if path.is_relative_to(BASE) else str(path),
        "is_control": int(is_control),
        "event_label": event_label,
        "total_hours": total_h,
        "pre_event_hours": pre_h,
        "n_beats": n_beats,
        "interp_frac": interp,
        "rr_method": str(d.get("rr_method", "")),
        "duration_stratum": stratum,
        "include": int(include),
        "exclusion_reason": ";".join(reasons) if reasons else "",
        "independence": "independent" if db in ("vfdb", "cudb", "nsrdb") else "internal_extension",
    }


def analyze_event_record(
    d: dict,
    *,
    z_threshold: float,
    min_consecutive: int,
    theta3: float,
) -> List[dict]:
    rr = np.asarray(d["rr_ms"], dtype=float)
    t_sec = np.asarray(d["t_sec"], dtype=float)
    t_hr = t_sec / 3600.0
    event_sec = float(d["vfon_sec"])
    event_hr = event_sec / 3600.0
    total_h = float(d.get("total_hours", np.nanmax(t_hr) if len(t_hr) else 0.0))
    db = str(d.get("source_db", "unknown"))
    rec = str(d.get("record_id", "?"))

    # use only pre-event samples for metrics (avoid VF-contaminated RR)
    pre_mask = t_hr < event_hr
    if pre_mask.sum() < W_TAU + 10:
        return []
    rr_pre = rr[pre_mask]
    t_pre = t_hr[pre_mask]
    event_hr_use, basal, approach, dur_stratum = short_db_windows(event_hr, event_hr)

    series: Dict[str, Tuple[np.ndarray, np.ndarray]] = {}
    tau = compute_tau_series(rr_pre)
    series["tau_s"] = (t_pre, tau)
    series["var"] = (t_pre, rolling_var(rr_pre, W_EWS if len(rr_pre) > W_EWS else max(51, len(rr_pre) // 5), STRIDE))
    series["ar1"] = (t_pre, rolling_ar1(rr_pre, W_EWS if len(rr_pre) > W_EWS else max(51, len(rr_pre) // 5), STRIDE))
    off, ex3 = compute_excess3_series(rr_pre, theta3=theta3)
    t_ex = t_pre[off : off + len(ex3)] if off < len(t_pre) else t_pre[:0]
    if len(t_ex) < len(ex3):
        ex3 = ex3[: len(t_ex)]
    elif len(t_ex) > len(ex3):
        t_ex = t_ex[: len(ex3)]
    series["excess3"] = (t_ex, ex3)

    rows = []
    for name, (t, m) in series.items():
        # align lengths
        n = min(len(t), len(m))
        t, m = t[:n], m[:n]
        det = detect_lead_time(
            m,
            t,
            event_hr_use,
            basal,
            z_threshold=z_threshold,
            min_consecutive=min_consecutive,
            use_abs=True,
        )
        rd = regime_delta(m, t, basal, approach)
        rows.append(
            {
                "source_db": db,
                "record_id": rec,
                "metric": name,
                "event_hr": event_hr_use,
                "event_label": str(d.get("event_label", "")),
                "duration_stratum": dur_stratum,
                "basal_start": basal[0],
                "basal_end": basal[1],
                "approach_start": approach[0],
                "approach_end": approach[1],
                "pre_event_hours": float(d.get("pre_event_hours", event_hr)),
                "total_hours": total_h,
                "z_threshold": z_threshold,
                "min_consecutive": min_consecutive,
                "theta3": theta3,
                "high_threshold": FROZEN_HIGH_THRESHOLD,
                "W_TAU": W_TAU,
                "stride": STRIDE,
                "delta": rd["delta"],
                "basal_mean_metric": rd["basal_mean"],
                "approach_mean_metric": rd["approach_mean"],
                "independence": "independent",
                **det,
            }
        )
    return rows


def analyze_control_record(
    d: dict,
    *,
    z_threshold: float,
    min_consecutive: int,
    theta3: float,
    basal_hours: float = 2.0,
    max_hours: Optional[float] = 12.0,
) -> List[dict]:
    """Run abs-z detector on control Holter; count FAR episodes for τ_s and excess3."""
    rr = np.asarray(d["rr_ms"], dtype=float)
    t_sec = np.asarray(d["t_sec"], dtype=float)
    t_hr = t_sec / 3600.0
    db = str(d.get("source_db", "unknown"))
    rec = str(d.get("record_id", "?"))
    total_h = float(t_hr[-1]) if len(t_hr) else 0.0
    if max_hours is not None and total_h > max_hours:
        # cap compute: first max_hours only (still multi-hour FAR estimate)
        keep = t_hr <= max_hours
        rr = rr[keep]
        t_hr = t_hr[keep]
        total_h = float(t_hr[-1]) if len(t_hr) else 0.0

    if len(rr) < W_TAU + 50 or total_h < basal_hours + 0.5:
        return []

    basal = (0.25, min(basal_hours, total_h * 0.25))
    if basal[1] <= basal[0]:
        basal = (0.0, max(total_h * 0.2, 0.1))

    series = {}
    series["tau_s"] = (t_hr, compute_tau_series(rr))
    off, ex3 = compute_excess3_series(rr, theta3=theta3)
    t_ex = t_hr[off : off + len(ex3)] if off < len(t_hr) else t_hr[:0]
    if len(t_ex) < len(ex3):
        ex3 = ex3[: len(t_ex)]
    series["excess3"] = (t_ex, ex3)

    rows = []
    for name, (t, m) in series.items():
        n = min(len(t), len(m))
        t, m = t[:n], m[:n]
        ep = count_alarm_episodes(
            m,
            t,
            basal,
            z_threshold=z_threshold,
            min_consecutive=min_consecutive,
            use_abs=True,
            refractory_h=0.5,
        )
        rows.append(
            {
                "source_db": db,
                "record_id": rec,
                "metric": name,
                "is_control": 1,
                "total_hours_used": total_h,
                "basal_start": basal[0],
                "basal_end": basal[1],
                "z_threshold": z_threshold,
                "min_consecutive": min_consecutive,
                "theta3": theta3,
                "high_threshold": FROZEN_HIGH_THRESHOLD,
                "W_TAU": W_TAU,
                "stride": STRIDE,
                **ep,
            }
        )
    return rows


def write_csv(path: Path, rows: List[dict]) -> None:
    if not rows:
        path.write_text("")
        return
    # union of keys
    keys: List[str] = []
    seen = set()
    for r in rows:
        for k in r.keys():
            if k not in seen:
                seen.add(k)
                keys.append(k)
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            # jsonify non-jsonable
            out = {}
            for k, v in r.items():
                if isinstance(v, (np.floating, float)):
                    out[k] = float(v)
                elif isinstance(v, (np.integer, int)):
                    out[k] = int(v)
                else:
                    out[k] = v
            w.writerow(out)


def main() -> int:
    ap = argparse.ArgumentParser(description="External Validation Phase 1 (frozen params)")
    ap.add_argument("--z-threshold", type=float, default=FROZEN_Z_THRESHOLD)
    ap.add_argument("--min-consecutive", type=int, default=FROZEN_MIN_CONSECUTIVE)
    ap.add_argument("--theta3", type=float, default=FROZEN_THETA3)
    ap.add_argument("--control-max-hours", type=float, default=12.0)
    ap.add_argument("--include-sddb-extension", action="store_true", help="also run SDDB rec 44 as labeled internal extension")
    args = ap.parse_args()

    # Guard: refuse silent retuning away from frozen defaults without explicit CLI
    frozen_params = {
        "theta3": FROZEN_THETA3,
        "high_threshold": FROZEN_HIGH_THRESHOLD,
        "W_TAU": W_TAU,
        "W_EWS": W_EWS,
        "stride": STRIDE,
        "z_threshold": FROZEN_Z_THRESHOLD,
        "min_consecutive": FROZEN_MIN_CONSECUTIVE,
        "alarm_rule": "abs_z_from_basal >= z_threshold sustained min_consecutive",
        "rr_clean_ms": [250.0, 2000.0],
        "note": "Defaults match Jul-12 2026 discovery; CLI may log overrides but Phase-1 report uses CLI values as run params.",
    }
    # Use CLI values (defaults = frozen)
    z_th = float(args.z_threshold)
    min_c = int(args.min_consecutive)
    theta3 = float(args.theta3)

    print("=== External Validation Phase 1 ===")
    print(f"params: theta3={theta3} z={z_th} min_c={min_c} W_TAU={W_TAU} stride={STRIDE}")
    print(f"frozen defaults: {frozen_params}")

    # Inventory all external npz
    inv_rows = []
    if RR_EXT.exists():
        for p in sorted(RR_EXT.glob("*_clean.npz")):
            inv_rows.append(inventory_row_from_npz(p))

    # Optional SDDB extension (record 44)
    if args.include_sddb_extension:
        p44 = DATA / "rr_44_clean.npz"
        if p44.exists():
            d = np.load(p44)
            inv_rows.append(
                {
                    "source_db": "sddb",
                    "record_id": "44",
                    "npz_path": "data/rr_44_clean.npz",
                    "is_control": 0,
                    "event_label": "VF",
                    "total_hours": float(d["total_hours"]),
                    "pre_event_hours": float(d["vfon_sec"]) / 3600.0,
                    "n_beats": int(d["n_beats"]),
                    "interp_frac": float(d["interp_frac"]) if "interp_frac" in d.files else float("nan"),
                    "rr_method": "beat_ann",
                    "duration_stratum": "holter_ge6h",
                    "include": 1,
                    "exclusion_reason": "",
                    "independence": "internal_extension",
                }
            )

    inv_path = RES / "external_phase1_inventory.csv"
    write_csv(inv_path, inv_rows)
    print(f"Wrote {inv_path} ({len(inv_rows)} rows)")

    # Process events
    event_rows: List[dict] = []
    for inv in inv_rows:
        if inv["is_control"] or not inv["include"]:
            continue
        path = BASE / inv["npz_path"] if not Path(inv["npz_path"]).is_absolute() else Path(inv["npz_path"])
        if inv["source_db"] == "sddb":
            # load SDDB schema npz
            raw = np.load(path)
            d = {
                "rr_ms": raw["rr_ms"],
                "t_sec": raw["t_sec"],
                "vfon_sec": float(raw["vfon_sec"]),
                "total_hours": float(raw["total_hours"]),
                "source_db": "sddb",
                "record_id": inv["record_id"],
                "event_label": inv["event_label"],
                "pre_event_hours": inv["pre_event_hours"],
            }
        else:
            d = load_npz(path)
        try:
            rows = analyze_event_record(d, z_threshold=z_th, min_consecutive=min_c, theta3=theta3)
            for r in rows:
                r["independence"] = inv["independence"]
            event_rows.extend(rows)
            print(f"  event {inv['source_db']}/{inv['record_id']}: {len(rows)} metric rows")
        except Exception as e:
            print(f"  event {inv['source_db']}/{inv['record_id']}: FAIL {e}")

    per_path = RES / "external_phase1_per_record.csv"
    write_csv(per_path, event_rows)
    print(f"Wrote {per_path} ({len(event_rows)} rows)")

    # Sensitivity by metric (independent only)
    indep = [r for r in event_rows if r.get("independence") == "independent"]
    metrics = sorted({r["metric"] for r in indep}) or ["tau_s", "excess3"]
    sens = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "params": {
            "theta3": theta3,
            "high_threshold": FROZEN_HIGH_THRESHOLD,
            "W_TAU": W_TAU,
            "stride": STRIDE,
            "z_threshold": z_th,
            "min_consecutive": min_c,
            "alarm_rule": frozen_params["alarm_rule"],
            "rr_clean_ms": [250.0, 2000.0],
        },
        "n_inventory": len(inv_rows),
        "n_included_events_indep": len({(r["source_db"], r["record_id"]) for r in indep}),
        "by_metric": {},
        "by_stratum_tau_s": {},
        "concordance_tau_excess3": {},
        "limitations": [
            "Phase-1 public DBs (VFDB/CU) are short vs SDDB Holters; lead times in minutes not multi-hour.",
            "Sensitivity on event-only sets; clinical PPV not estimated.",
            "No threshold retuning; frozen Jul-12 discovery params.",
            "Not a multi-center clinical validation; no FDA/deployability claim.",
        ],
    }
    for met in metrics:
        subset = [r for r in indep if r["metric"] == met]
        perf = detector_performance(subset)
        sens["by_metric"][met] = perf
        print(
            f"  sens {met}: n={perf['n_records']} det={perf['n_detected']} "
            f"sens={perf['sensitivity']} med_lead_h={perf['median_lead_time_h']}"
        )

    # stratum for tau_s
    for stratum in sorted({r.get("duration_stratum", "?") for r in indep}):
        sub = [r for r in indep if r["metric"] == "tau_s" and r.get("duration_stratum") == stratum]
        if sub:
            sens["by_stratum_tau_s"][stratum] = detector_performance(sub)

    # concordance of deltas
    tau_d = { (r["source_db"], r["record_id"]): r["delta"] for r in indep if r["metric"] == "tau_s" }
    ex_d = { (r["source_db"], r["record_id"]): r["delta"] for r in indep if r["metric"] == "excess3" }
    keys = sorted(set(tau_d) & set(ex_d))
    if keys:
        sens["concordance_tau_excess3"] = sign_concordance(
            [tau_d[k] for k in keys], [ex_d[k] for k in keys]
        )

    sens_path = RES / "external_phase1_sensitivity.json"
    with sens_path.open("w") as f:
        json.dump(sens, f, indent=2, default=float)
    print(f"Wrote {sens_path}")

    # Controls / FAR
    control_rows: List[dict] = []
    for inv in inv_rows:
        if not inv["is_control"] or not inv["include"]:
            continue
        path = BASE / inv["npz_path"] if not Path(inv["npz_path"]).is_absolute() else Path(inv["npz_path"])
        d = load_npz(path)
        try:
            rows = analyze_control_record(
                d,
                z_threshold=z_th,
                min_consecutive=min_c,
                theta3=theta3,
                max_hours=args.control_max_hours,
            )
            control_rows.extend(rows)
            print(f"  control {inv['source_db']}/{inv['record_id']}: {len(rows)} metric rows")
        except Exception as e:
            print(f"  control {inv['source_db']}/{inv['record_id']}: FAIL {e}")

    ctrl_path = RES / "external_phase1_controls.csv"
    write_csv(ctrl_path, control_rows)
    print(f"Wrote {ctrl_path} ({len(control_rows)} rows)")

    far_report = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "params": sens["params"],
        "control_source": "PhysioNet NSRDB (negative-control Holters; no VF)",
        "by_metric": {},
        "notes": [
            "FAR = total alarm episodes / total search hours * 24.",
            "Basal = early control window; search = remainder (capped).",
            "NSRDB is rhythm-healthy Holter — not device-matched to VFDB/CU telemetry.",
            "If no controls processed, far_per_24h is nan (honest fallback).",
        ],
    }
    for met in ("tau_s", "excess3"):
        sub = [r for r in control_rows if r["metric"] == met]
        far_report["by_metric"][met] = false_alarm_rate(sub)
        fr = far_report["by_metric"][met]
        print(
            f"  FAR {met}: n_ctrl={fr['n_controls']} hours={fr['total_search_hours']:.2f} "
            f"episodes={fr['total_episodes']} FAR/24h={fr['far_per_24h']} reason={fr['reason']}"
        )

    far_path = RES / "external_phase1_far.json"
    with far_path.open("w") as f:
        json.dump(far_report, f, indent=2, default=float)
    print(f"Wrote {far_path}")

    # Master summary
    summary = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "phase": "external_validation_phase1",
        "status": "started_and_reported",
        "params_frozen": sens["params"],
        "inventory_csv": str(inv_path.relative_to(BASE)),
        "per_record_csv": str(per_path.relative_to(BASE)),
        "sensitivity_json": str(sens_path.relative_to(BASE)),
        "controls_csv": str(ctrl_path.relative_to(BASE)),
        "far_json": str(far_path.relative_to(BASE)),
        "n_independent_events_processable": sens["n_included_events_indep"],
        "sensitivity_by_metric": sens["by_metric"],
        "far_by_metric": far_report["by_metric"],
        "clinical_claim": False,
        "deployability_claim": False,
        "next_step": (
            "Expand processable independent n (more CU/VFDB with pre≥15–60 min); "
            "add institutional matched non-event Holters for device-matched FAR (Phase 2); "
            "do not retune thresholds."
        ),
        "success_criteria_note": (
            "S1–S6 from docs/EXTERNAL_VALIDATION_PLAN.md are NOT claimed fully met; "
            "Phase 1 provides preliminary external sensitivity + control FAR estimate only."
        ),
    }
    sum_path = RES / "external_phase1_summary.json"
    with sum_path.open("w") as f:
        json.dump(summary, f, indent=2, default=float)
    print(f"Wrote {sum_path}")
    print("=== Phase 1 run complete (no clinical deployability claim) ===")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
