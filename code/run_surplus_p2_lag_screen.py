#!/usr/bin/env python3
"""
Phase β mini — P2 lag screen (diagnostic only).

Compares three frozen lags of the P2 proxy under identical Bandt–Pompe /
TV surplus settings as Phase α:

  P2-5   ℓ=5
  P2-10  ℓ=10  (Phase-α reference)
  P2-20  ℓ=20

Not a detector, not a continuous lag grid, not I0.

Outputs under results/:
  surplus_p2lag_nsrdb_per_record.csv
  surplus_p2lag_event_per_record.csv
  surplus_p2lag_summary.csv
  surplus_p2lag_summary.json

Reproduce:
  PYTHONPATH=code python3 code/run_surplus_p2_lag_screen.py
  PYTHONPATH=code python3 code/run_surplus_p2_lag_screen.py --smoke

Observational only — no clinical claims.
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from proxy_variants import (  # noqa: E402
    P2_LAG_DESCRIPTIONS,
    P2_LAG_IDS,
    P2_LAG_SCREEN,
)
from run_surplus_proxy_bakeoff import (  # noqa: E402
    CONTROL_MAX_HOURS,
    DELAY,
    HIGH_S_ABS,
    K_JOINT,
    L_SURPLUS,
    M_EMB,
    MERGE_GAP,
    MIN_EPISODE_LEN,
    RES,
    aggregate_proxy,
    list_event_records,
    list_nsrdb,
    process_event_record,
    process_nsrdb_record,
    write_csv,
)

REF_ID = "P2-10"  # Phase-α lag reference within this screen


def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser(description="P2 lag screen diagnostic (ℓ=5,10,20)")
    ap.add_argument("--smoke", action="store_true", help="2 NSRDB + 2 events")
    ap.add_argument("--L", type=int, default=L_SURPLUS)
    ap.add_argument("--max-hours", type=float, default=CONTROL_MAX_HOURS)
    ap.add_argument("--high-abs", type=float, default=HIGH_S_ABS)
    ap.add_argument(
        "--lags",
        type=str,
        default=",".join(str(x) for x in P2_LAG_SCREEN),
        help="Comma-separated lag values (default: 5,10,20)",
    )
    args = ap.parse_args(argv)

    lags = [int(x.strip()) for x in args.lags.split(",") if x.strip()]
    proxy_ids = [f"P2-{lag}" for lag in lags]
    for pid in proxy_ids:
        if pid not in P2_LAG_IDS and not pid.startswith("P2-"):
            print(f"Unexpected proxy id {pid}", file=sys.stderr)
            return 2
        # Allow only lags registered in P2_LAG_SCREEN (or dynamic build via build_proxy_p2_lagged)
        # build_proxy only knows registered P2_LAG_IDS; reject others
        from proxy_variants import ALL_PROXY_BUILDERS

        if pid not in ALL_PROXY_BUILDERS:
            print(
                f"Unknown lag proxy {pid}; registered: {list(P2_LAG_IDS)}",
                file=sys.stderr,
            )
            return 2

    nsrdb_paths = list_nsrdb()
    events = list_event_records()
    if args.smoke:
        nsrdb_paths = nsrdb_paths[:2]
        events = events[:2]

    print(f"[p2lag] lags={lags} proxies={proxy_ids}")
    print(f"[p2lag] NSRDB n={len(nsrdb_paths)} events n={len(events)}")
    print(
        f"[p2lag] L={args.L} high_abs={args.high_abs} max_h={args.max_hours} "
        f"ref={REF_ID}"
    )

    nsrdb_rows: List[dict] = []
    event_rows: List[dict] = []

    for pid in proxy_ids:
        desc = P2_LAG_DESCRIPTIONS.get(pid, pid)
        print(f"[p2lag] --- {pid}: {desc}")
        for path in nsrdb_paths:
            row = process_nsrdb_record(
                path, pid, L=args.L, max_hours=args.max_hours, high_abs=args.high_abs
            )
            nsrdb_rows.append(row)
            if not row.get("skipped"):
                print(
                    f"  nsrdb/{row['record_id']}: meanS={row['mean_S']:.4f} "
                    f"p90={row['p90_S']:.4f} f40={row['frac_S_ge_0p40']:.3f} "
                    f"eps24={row['episode_rate_per_24h']:.1f}"
                )
            else:
                print(f"  nsrdb/{row['record_id']}: SKIP {row.get('skip_reason')}")

        for source, rec, path in events:
            row = process_event_record(source, rec, path, pid, L=args.L)
            event_rows.append(row)
            if not row.get("skipped"):
                print(
                    f"  {source}/{rec}: basal={row['mean_S_basal']:.4f} "
                    f"app={row['mean_S_approach']:.4f} "
                    f"Δ={row['delta_mean_approach_minus_basal']:.4f}"
                )
            else:
                print(f"  {source}/{rec}: SKIP {row.get('skip_reason')}")

    summaries = [aggregate_proxy(pid, nsrdb_rows, event_rows) for pid in proxy_ids]
    by_proxy = {s["proxy_id"]: s for s in summaries}
    ref = by_proxy.get(REF_ID)

    for s in summaries:
        lag = int(str(s["proxy_id"]).split("-")[-1]) if "-" in s["proxy_id"] else -1
        s["lag"] = lag
        if ref and s["proxy_id"] != REF_ID:
            rm = float(ref["nsrdb_mean_of_mean_S"])
            rf = float(ref["nsrdb_mean_frac_S_ge_0p40"])
            rd = float(ref["cohens_d_event_basal_vs_nsrdb"])
            m = float(s["nsrdb_mean_of_mean_S"])
            f = float(s["nsrdb_mean_frac_S_ge_0p40"])
            d = float(s["cohens_d_event_basal_vs_nsrdb"])
            s["rel_mean_S_vs_p2_10"] = m / rm if rm > 1e-12 else float("nan")
            s["mean_S_drop_frac_vs_p2_10"] = (rm - m) / rm if rm > 1e-12 else float("nan")
            s["frac40_drop_frac_vs_p2_10"] = (rf - f) / rf if rf > 1e-12 else float("nan")
            s["delta_cohens_d_vs_p2_10"] = d - rd if np.isfinite(d) and np.isfinite(rd) else float("nan")
        else:
            s["rel_mean_S_vs_p2_10"] = 1.0
            s["mean_S_drop_frac_vs_p2_10"] = 0.0
            s["frac40_drop_frac_vs_p2_10"] = 0.0
            s["delta_cohens_d_vs_p2_10"] = 0.0
        # Balance notes (descriptive only)
        s["event_basal_ge_nsrdb"] = bool(s.get("event_basal_gt_nsrdb")) or (
            np.isfinite(float(s.get("diff_event_basal_minus_nsrdb_mean", np.nan)))
            and float(s["diff_event_basal_minus_nsrdb_mean"]) >= -1e-6
        )

    p_nsr = RES / "surplus_p2lag_nsrdb_per_record.csv"
    p_ev = RES / "surplus_p2lag_event_per_record.csv"
    p_sum = RES / "surplus_p2lag_summary.csv"
    p_json = RES / "surplus_p2lag_summary.json"

    write_csv(p_nsr, nsrdb_rows)
    write_csv(p_ev, event_rows)
    write_csv(p_sum, summaries)

    payload: Dict[str, Any] = {
        "exercise": "surplus_p2_lag_screen",
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "smoke": bool(args.smoke),
        "methods": {
            "m": M_EMB,
            "delay": DELAY,
            "L": args.L,
            "K": K_JOINT,
            "statistic": "TV(P_joint, P1⊗P2)",
            "high_abs": args.high_abs,
            "min_episode_len": MIN_EPISODE_LEN,
            "merge_gap": MERGE_GAP,
            "control_max_hours": args.max_hours,
            "lags": lags,
            "reference_within_screen": REF_ID,
            "proxies": {
                pid: P2_LAG_DESCRIPTIONS.get(pid, pid) for pid in proxy_ids
            },
            "note": (
                "Discrete lag screen only. Same BP/TV as Phase α. "
                "Not a detector; not continuous lag optimization."
            ),
        },
        "lag_summaries": summaries,
        "artifacts": {
            "nsrdb_per_record": str(p_nsr.name),
            "event_per_record": str(p_ev.name),
            "summary": str(p_sum.name),
        },
        "disclaimer": "Diagnostic only. No clinical claims. Not a detector.",
    }
    p_json.write_text(json.dumps(payload, indent=2, default=str))

    print("\n[p2lag] === LAG SUMMARY ===")
    hdr = (
        f"{'id':6s} {'lag':>3s} {'meanS':>7s} {'p90':>7s} {'f≥0.40':>7s} "
        f"{'eps/24h':>8s} {'evBas':>7s} {'Δapp':>7s} {'d_ev':>6s}"
    )
    print(hdr)
    for s in summaries:
        print(
            f"{s['proxy_id']:6s} "
            f"{int(s.get('lag', -1)):3d} "
            f"{s['nsrdb_mean_of_mean_S']:7.4f} "
            f"{s['nsrdb_mean_of_p90_S']:7.4f} "
            f"{s['nsrdb_mean_frac_S_ge_0p40']:7.3f} "
            f"{s['nsrdb_mean_episode_rate_per_24h']:8.1f} "
            f"{s['event_mean_of_mean_S_basal']:7.4f} "
            f"{s['event_mean_delta_approach_minus_basal']:7.4f} "
            f"{s['cohens_d_event_basal_vs_nsrdb']:6.2f}"
        )
    print(f"[p2lag] wrote {p_nsr.name}, {p_ev.name}, {p_sum.name}, {p_json.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
