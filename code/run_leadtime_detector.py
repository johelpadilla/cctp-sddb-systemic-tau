#!/usr/bin/env python3
"""
Lead-time and detector performance for τ_s, excess3, var, and AR(1).

Uses manuscript windows (basal vs pre-event search) and z-threshold alarms on
rolling metric series computed from cleaned RR npz files.

Outputs:
  results/leadtime_per_record.csv
  results/leadtime_detector_summary.json
  results/leadtime_cumulative.csv
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
from _bootstrap import import_systemictau_core
from cctp_metrics_core import (
    STRIDE,
    W_EWS,
    W_TAU,
    build_bivariate_proxy,
    cumulative_detection_curve,
    detect_lead_time,
    detector_performance,
    get_event_and_windows,
    rolling_ar1,
    rolling_var,
)
from recd_ordinal_levels import compute_phi3, generate_multivariate_symbols

BASE = Path(__file__).resolve().parent.parent
DATA = BASE / "data"
RES = BASE / "results"
RES.mkdir(parents=True, exist_ok=True)

ANALYTIC = ["30", "31", "32", "35", "36", "38", "45", "47", "50", "51"]


def compute_tau_series(rr: np.ndarray) -> np.ndarray:
    compute_taus, _, has = import_systemictau_core()
    X = build_bivariate_proxy(rr)
    if has and compute_taus is not None:
        taus_global, _ = compute_taus(X, window_size=W_TAU, stride=STRIDE)
        return np.asarray(taus_global, dtype=float)
    # fallback: rolling Spearman of the two proxy columns
    n = len(rr)
    out = np.full(n, np.nan)
    for i in range(W_TAU - 1, n, STRIDE):
        win = X[i - W_TAU + 1 : i + 1]
        r0 = win[:, 0].argsort().argsort().astype(float)
        r1 = win[:, 1].argsort().argsort().astype(float)
        if np.std(r0) > 0 and np.std(r1) > 0:
            out[i] = np.corrcoef(r0, r1)[0, 1]
    return out


def compute_excess3_series(rr: np.ndarray, theta3: float = 0.08) -> tuple:
    X = build_bivariate_proxy(rr)
    S = generate_multivariate_symbols(X, m=3, delay=1)
    offset = (3 - 1) * 1
    phi3, excess3 = compute_phi3(S, window=W_TAU, theta=theta3, stride=STRIDE)
    return offset, np.asarray(excess3, dtype=float)


def analyze_record(rec: str, z_threshold: float, min_consecutive: int, theta3: float) -> list:
    npz = DATA / f"rr_{rec}_clean.npz"
    if not npz.exists():
        print(f"  skip {rec}: missing {npz.name}")
        return []
    d = np.load(npz)
    rr = d["rr_ms"].astype(float)
    t_sec = d["t_sec"].astype(float)
    vfon_sec = float(d["vfon_sec"])
    t_hr = t_sec / 3600.0
    vfon_hr = vfon_sec / 3600.0
    event_hr, basal, approach = get_event_and_windows(rec, t_hr, vfon_hr)

    series = {}
    series["tau_s"] = (t_hr, compute_tau_series(rr))
    series["var"] = (t_hr, rolling_var(rr, W_EWS, STRIDE))
    series["ar1"] = (t_hr, rolling_ar1(rr, W_EWS, STRIDE))
    off, ex3 = compute_excess3_series(rr, theta3=theta3)
    t_ex = t_hr[off : off + len(ex3)]
    if len(t_ex) < len(ex3):
        ex3 = ex3[: len(t_ex)]
    elif len(t_ex) > len(ex3):
        t_ex = t_ex[: len(ex3)]
    series["excess3"] = (t_ex, ex3)

    rows = []
    for name, (t, m) in series.items():
        det = detect_lead_time(
            m,
            t,
            event_hr,
            basal,
            z_threshold=z_threshold,
            min_consecutive=min_consecutive,
            use_abs=True,
        )
        rows.append(
            {
                "record": rec,
                "metric": name,
                "event_hr": event_hr,
                "basal_start": basal[0],
                "basal_end": basal[1],
                "approach_start": approach[0],
                "approach_end": approach[1],
                "z_threshold": z_threshold,
                "min_consecutive": min_consecutive,
                **det,
            }
        )
        print(
            f"  {rec} {name}: alarmed={int(det['alarmed'])} "
            f"lead_h={det['lead_time_h']:.3f}" if np.isfinite(det["lead_time_h"]) else
            f"  {rec} {name}: alarmed={int(det['alarmed'])} lead_h=nan"
        )
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--records", default=",".join(ANALYTIC))
    ap.add_argument("--z-threshold", type=float, default=2.0)
    ap.add_argument("--min-consecutive", type=int, default=3)
    ap.add_argument("--theta3", type=float, default=0.08)
    ap.add_argument("--out-csv", default=str(RES / "leadtime_per_record.csv"))
    ap.add_argument("--out-summary", default=str(RES / "leadtime_detector_summary.json"))
    ap.add_argument("--out-cum", default=str(RES / "leadtime_cumulative.csv"))
    args = ap.parse_args()

    recs = [r.strip() for r in args.records.split(",") if r.strip()]
    print(f"Lead-time / detector analysis: records={recs}")
    print(f"z_threshold={args.z_threshold} min_consecutive={args.min_consecutive} theta3={args.theta3}")

    all_rows = []
    for rec in recs:
        all_rows.extend(
            analyze_record(rec, args.z_threshold, args.min_consecutive, args.theta3)
        )

    if not all_rows:
        raise SystemExit("No lead-time rows produced")

    fieldnames = list(all_rows[0].keys())
    with open(args.out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(all_rows)
    print(f"Wrote {args.out_csv} ({len(all_rows)} rows)")

    # per-metric detector performance + cumulative curves
    metrics = sorted({r["metric"] for r in all_rows})
    summary = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "params": {
            "z_threshold": args.z_threshold,
            "min_consecutive": args.min_consecutive,
            "theta3": args.theta3,
            "records": recs,
            "W_TAU": W_TAU,
            "W_EWS": W_EWS,
            "stride": STRIDE,
            "alarm_rule": "abs_z_from_basal >= z_threshold sustained min_consecutive",
        },
        "by_metric": {},
        "cumulative_by_metric": {},
    }
    cum_rows = []
    for met in metrics:
        subset = [r for r in all_rows if r["metric"] == met]
        perf = detector_performance(subset)
        summary["by_metric"][met] = perf
        curve = cumulative_detection_curve(subset)
        summary["cumulative_by_metric"][met] = curve
        for c in curve:
            cum_rows.append({"metric": met, **c})
        print(
            f"  {met}: sens={perf['sensitivity']:.2f} "
            f"median_lead={perf['median_lead_time_h']} n_det={perf['n_detected']}/{perf['n_records']}"
        )

    with open(args.out_summary, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"Wrote {args.out_summary}")

    with open(args.out_cum, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["metric", "horizon_h", "detection_rate", "n_detected", "n_records"])
        w.writeheader()
        w.writerows(cum_rows)
    print(f"Wrote {args.out_cum}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
