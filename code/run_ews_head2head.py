#!/usr/bin/env python3
"""
Head-to-head comparison: Systemic Tau + excess3 vs classic var / AR(1).

Uses:
  - results/cctp_batch_summary.csv (basal→approach deltas)
  - results/leadtime_per_record.csv + leadtime_detector_summary.json (detector/lead-time)

Outputs:
  results/ews_head2head.csv
  results/ews_head2head_report.json
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cctp_metrics_core import detector_performance, sign_concordance

BASE = Path(__file__).resolve().parent.parent
RES = BASE / "results"


def load_csv(path: Path) -> list:
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def fnum(v):
    if v in (None, ""):
        return float("nan")
    try:
        return float(v)
    except ValueError:
        return float("nan")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--batch", default=str(RES / "cctp_batch_summary.csv"))
    ap.add_argument("--leadtime", default=str(RES / "leadtime_per_record.csv"))
    ap.add_argument("--lead-summary", default=str(RES / "leadtime_detector_summary.json"))
    ap.add_argument("--out-csv", default=str(RES / "ews_head2head.csv"))
    ap.add_argument("--out-json", default=str(RES / "ews_head2head_report.json"))
    args = ap.parse_args()

    batch = load_csv(Path(args.batch))
    lead_path = Path(args.leadtime)
    lead_rows = load_csv(lead_path) if lead_path.exists() else []
    lead_sum = {}
    if Path(args.lead_summary).exists():
        with open(args.lead_summary) as f:
            lead_sum = json.load(f)

    # per-record comparative table
    out_rows = []
    for m in batch:
        rec = str(m["record"])
        dt = fnum(m.get("delta_tau"))
        de = fnum(m.get("delta_excess3"))
        dv = fnum(m.get("delta_var"))
        da = fnum(m.get("delta_ar1"))

        def lead_for(metric):
            matches = [r for r in lead_rows if r.get("record") == rec and r.get("metric") == metric]
            if not matches:
                return float("nan"), 0.0
            r0 = matches[0]
            return fnum(r0.get("lead_time_h")), fnum(r0.get("alarmed"))

        lt_tau, al_tau = lead_for("tau_s")
        lt_ex, al_ex = lead_for("excess3")
        lt_var, al_var = lead_for("var")
        lt_ar1, al_ar1 = lead_for("ar1")

        out_rows.append(
            {
                "record": rec,
                "delta_tau": dt,
                "delta_excess3": de,
                "delta_var": dv,
                "delta_ar1": da,
                "sign_tau": float(np.sign(dt)) if np.isfinite(dt) else float("nan"),
                "sign_excess3": float(np.sign(de)) if np.isfinite(de) else float("nan"),
                "sign_var": float(np.sign(dv)) if np.isfinite(dv) else float("nan"),
                "sign_ar1": float(np.sign(da)) if np.isfinite(da) else float("nan"),
                "tau_excess3_concordant": int(
                    np.isfinite(dt) and np.isfinite(de) and (np.sign(dt) == np.sign(de)) and np.sign(dt) != 0
                ),
                "tau_var_concordant": int(
                    np.isfinite(dt) and np.isfinite(dv) and (np.sign(dt) == np.sign(dv)) and np.sign(dt) != 0
                ),
                "lead_time_tau_h": lt_tau,
                "lead_time_excess3_h": lt_ex,
                "lead_time_var_h": lt_var,
                "lead_time_ar1_h": lt_ar1,
                "alarmed_tau": al_tau,
                "alarmed_excess3": al_ex,
                "alarmed_var": al_var,
                "alarmed_ar1": al_ar1,
                "p_tau": fnum(m.get("p_tau")),
                "p_excess3": fnum(m.get("p_excess3")),
                "p_tau_surrogate": fnum(m.get("p_tau_surrogate")),
                "interp_frac": fnum(m.get("interp_frac")),
                "pacing_detected": m.get("pacing_detected"),
            }
        )

    with open(args.out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(out_rows[0].keys()))
        w.writeheader()
        w.writerows(out_rows)
    print(f"Wrote {args.out_csv} ({len(out_rows)} rows)")

    dt_all = [r["delta_tau"] for r in out_rows]
    de_all = [r["delta_excess3"] for r in out_rows]
    dv_all = [r["delta_var"] for r in out_rows]
    da_all = [r["delta_ar1"] for r in out_rows]

    def metric_effect(name, vals):
        a = np.array(vals, dtype=float)
        a = a[np.isfinite(a)]
        return {
            "n": int(len(a)),
            "mean_delta": float(np.mean(a)) if len(a) else float("nan"),
            "median_delta": float(np.median(a)) if len(a) else float("nan"),
            "n_positive": int(np.sum(a > 0)) if len(a) else 0,
            "n_negative": int(np.sum(a < 0)) if len(a) else 0,
            "frac_positive": float(np.mean(a > 0)) if len(a) else float("nan"),
        }

    # detector side-by-side from lead summary or recompute
    by_metric = lead_sum.get("by_metric") or {}
    if not by_metric and lead_rows:
        for met in ("tau_s", "excess3", "var", "ar1"):
            by_metric[met] = detector_performance([r for r in lead_rows if r.get("metric") == met])

    report = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "n_records": len(out_rows),
        "effect_sizes_deltas": {
            "tau_s": metric_effect("tau_s", dt_all),
            "excess3": metric_effect("excess3", de_all),
            "var": metric_effect("var", dv_all),
            "ar1": metric_effect("ar1", da_all),
        },
        "direction_concordance": {
            "tau_vs_excess3": sign_concordance(dt_all, de_all),
            "tau_vs_var": sign_concordance(dt_all, dv_all),
            "tau_vs_ar1": sign_concordance(dt_all, da_all),
            "excess3_vs_var": sign_concordance(de_all, dv_all),
            "var_vs_ar1": sign_concordance(dv_all, da_all),
        },
        "detector_performance": by_metric,
        "interpretation_notes": [
            "All metrics share the same basal/approach windows as the manuscript pipeline.",
            "Direction concordance is sign(Δ) agreement; τ_s and excess3 are expected to align (relational).",
            "Classic var often rises (critical slowing-like) while AR(1) often falls on these Holters — relational metrics need not match univariate sign.",
            "Detector sensitivity uses abs-z departure from basal; FAR requires independent control Holters (reported nan until external validation).",
        ],
    }
    with open(args.out_json, "w") as f:
        json.dump(report, f, indent=2)
    print(f"Wrote {args.out_json}")

    c = report["direction_concordance"]["tau_vs_excess3"]
    print(
        f"τ_s vs excess3 concordance: {c['n_concordant']}/{c['n']} = {c['concordance']}"
    )
    for met, perf in by_metric.items():
        print(
            f"  detector {met}: sens={perf.get('sensitivity')} "
            f"median_lead_h={perf.get('median_lead_time_h')}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
