#!/usr/bin/env python3
"""
Structural I0 surplus arms (Track S1) — fixed 6 arms, no θ-grid.

R0 I0-ref:     L=50, mean+0.08, hop=1, θ_R=5
R1 hop:        L=50, mean+0.08, hop=25, θ_R=3
R2 basal-q90:  L=50, thr=q90(basal), hop=1, θ_R=5
R3 basal-MAD:  L=50, med+2.5·MAD, hop=1, θ_R=5
R4 combo:      L=50, thr=q90, hop=25, θ_R=3
R5 combo+L_S:  L=100, thr=q90, hop=50, θ_R=3

Protocol identical to OPSP Holter / I0 grid (SDDB+VFDB sens, NSRDB FAR).

Pre-registered structural success (relaxed vs original I0 grid clear-advance):
  FAR/24h ≤ 2 × OPC (≤ ~7.467) AND sens_all ≥ 0.55
Original clear-advance (sens ≥ 0.65 + FAR ≤ 2×OPC) still scored as jackpot.

Usage:
  PYTHONPATH=code python3 code/run_i0_structural_arms.py
  PYTHONPATH=code python3 code/run_i0_structural_arms.py --smoke
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from run_opsp_integrated_holter_eval import (  # noqa: E402
    BASELINE_ANCHORS,
    list_event_records,
    list_nsrdb,
    process_control,
    process_event_record,
    sens_block,
    far_block,
    write_csv,
    BASAL_HOURS,
    CONTROL_MAX_HOURS,
    REFRACTORY_H,
    RES,
)

OPC_SENS_ALL = float(BASELINE_ANCHORS["opc_L50"]["sens_all"])
OPC_FAR = float(BASELINE_ANCHORS["opc_L50"]["far_per_24h"])
ABSZ_FAR = float(BASELINE_ANCHORS["absz_tau_s"]["far_per_24h"])
FAR_CAP_2X_OPC = 2.0 * OPC_FAR
SENS_STRUCTURAL = 0.55
SENS_CLEAR = 0.65

PRE_REGISTERED = {
    "structural_win": {
        "far_le": FAR_CAP_2X_OPC,
        "sens_ge": SENS_STRUCTURAL,
        "description": "FAR ≤ 2×OPC and sens_all ≥ 0.55",
    },
    "clear_advance_original": {
        "far_le": FAR_CAP_2X_OPC,
        "sens_ge": SENS_CLEAR,
        "description": "FAR ≤ 2×OPC and sens_all ≥ 0.65 (original I0 grid jackpot)",
    },
    "secondary_interesting": {
        "far_le": 0.5 * ABSZ_FAR,
        "sens_ge": SENS_CLEAR,
        "description": "sens ≥ 0.65 and FAR ≤ 0.5×abs-z",
    },
    "stop_rule": (
        "If no arm reaches structural_win nor approaches "
        "(FAR ≤ 12 and sens ≥ 0.55), stop surplus-primary competitive track."
    ),
}


def build_structural_arms() -> List[dict]:
    """Fixed R0–R5 configs (no Cartesian product)."""
    common = dict(
        collapse_role="none",
        theta_D=0.35,
        confirm_window=5,
        modulate_delta_R=2,
        modulate_delta_S=0.02,
    )
    return [
        {
            **common,
            "name": "R0_I0_ref",
            "L": 50,
            "theta_delta_S": 0.08,
            "theta_R": 5,
            "basal_mode": "mean_delta",
            "hop": 1,
        },
        {
            **common,
            "name": "R1_hop",
            "L": 50,
            "theta_delta_S": 0.08,
            "theta_R": 3,
            "basal_mode": "mean_delta",
            "hop": 25,
        },
        {
            **common,
            "name": "R2_basal_q90",
            "L": 50,
            "theta_delta_S": None,
            "theta_R": 5,
            "basal_mode": "percentile",
            "basal_q": 90.0,
            "hop": 1,
        },
        {
            **common,
            "name": "R3_basal_mad",
            "L": 50,
            "theta_delta_S": None,
            "theta_R": 5,
            "basal_mode": "mad",
            "mad_kappa": 2.5,
            "hop": 1,
        },
        {
            **common,
            "name": "R4_combo",
            "L": 50,
            "theta_delta_S": None,
            "theta_R": 3,
            "basal_mode": "percentile",
            "basal_q": 90.0,
            "hop": 25,
        },
        {
            **common,
            "name": "R5_combo_Ls100",
            "L": 100,
            "theta_delta_S": None,
            "theta_R": 3,
            "basal_mode": "percentile",
            "basal_q": 90.0,
            "hop": 50,
        },
    ]


def score_arm(sens_all: float, far: float) -> Dict[str, Any]:
    """Score one arm against pre-registered structural rules."""
    s = float(sens_all)
    f = float(far)
    structural = bool(
        np.isfinite(s)
        and np.isfinite(f)
        and f <= FAR_CAP_2X_OPC + 1e-12
        and s >= SENS_STRUCTURAL - 1e-12
    )
    clear = bool(
        np.isfinite(s)
        and np.isfinite(f)
        and f <= FAR_CAP_2X_OPC + 1e-12
        and s >= SENS_CLEAR - 1e-12
    )
    secondary = bool(
        np.isfinite(s)
        and np.isfinite(f)
        and f <= 0.5 * ABSZ_FAR + 1e-12
        and s >= SENS_CLEAR - 1e-12
    )
    approaches = bool(
        np.isfinite(s)
        and np.isfinite(f)
        and f <= 12.0 + 1e-12
        and s >= SENS_STRUCTURAL - 1e-12
    )
    reasons: List[str] = []
    if clear:
        reasons.append("clear_advance_original (jackpot)")
    elif structural:
        reasons.append("structural_win")
    elif secondary:
        reasons.append("secondary_interesting")
    elif approaches:
        reasons.append("approaches (FAR≤12, sens≥0.55) but not structural_win")
    else:
        if not np.isfinite(f) or f > FAR_CAP_2X_OPC:
            reasons.append(f"FAR={f:.3f} > 2×OPC={FAR_CAP_2X_OPC:.3f}")
        if not np.isfinite(s) or s < SENS_STRUCTURAL:
            reasons.append(f"sens={s:.3f} < structural bar {SENS_STRUCTURAL}")
    return {
        "structural_win": structural,
        "clear_advance_original": clear,
        "secondary_interesting": secondary,
        "approaches": approaches,
        "far_over_opc": float(f / OPC_FAR) if np.isfinite(f) and OPC_FAR > 0 else float("nan"),
        "reasons": reasons,
    }


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(description="I0 structural surplus arms Holter bake-off")
    ap.add_argument("--control-max-hours", type=float, default=CONTROL_MAX_HOURS)
    ap.add_argument("--basal-hours", type=float, default=BASAL_HOURS)
    ap.add_argument("--refractory-h", type=float, default=REFRACTORY_H)
    ap.add_argument("--smoke", action="store_true")
    args = ap.parse_args(argv)

    configs = build_structural_arms()
    events = list_event_records()
    controls = list_nsrdb()
    if args.smoke:
        sddb = [e for e in events if e[0] == "sddb"][:2]
        vfdb = [e for e in events if e[0] == "vfdb"][:2]
        events = sddb + vfdb
        controls = controls[:2]
        print(f"SMOKE: {len(events)} events, {len(controls)} controls")

    print(f"Event records: {len(events)}  NSRDB: {len(controls)}")
    print(f"Arms: {[c['name'] for c in configs]}")

    sens_rows: List[dict] = []
    for source, rec, path in events:
        print(f"  event {source}/{rec} ...", flush=True)
        rows = process_event_record(source, rec, path, configs)
        sens_rows.extend(rows)
        for r in rows:
            if r.get("error"):
                print(f"    {r['detector']}: ERROR {r['error']}")
            else:
                print(
                    f"    {r['detector']}: alarmed={r['alarmed']} "
                    f"lead={r.get('lead_time_h', float('nan'))}"
                )

    far_rows: List[dict] = []
    for path in controls:
        print(f"  control {path.name} ...", flush=True)
        rows = process_control(
            path,
            configs,
            basal_hours=args.basal_hours,
            max_hours=args.control_max_hours,
            refractory_h=args.refractory_h,
        )
        far_rows.extend(rows)
        for r in rows:
            if r.get("skipped"):
                print(f"    {r['detector']}: skipped {r.get('skip_reason')}")
            else:
                print(
                    f"    {r['detector']}: FAR={r['far_per_24h']:.2f} "
                    f"ep={r['n_episodes']}"
                )

    det_names = [c["name"] for c in configs]
    comparison: List[dict] = []
    n_structural = 0
    n_clear = 0
    n_approach = 0

    for cfg in configs:
        name = cfg["name"]
        s_all = sens_block(sens_rows, name, None)
        s_sddb = sens_block(sens_rows, name, "sddb")
        s_vfdb = sens_block(sens_rows, name, "vfdb")
        f = far_block(far_rows, name)
        score = score_arm(s_all["sensitivity"], f["far_per_24h"])
        if score["structural_win"]:
            n_structural += 1
        if score["clear_advance_original"]:
            n_clear += 1
        if score["approaches"] or score["structural_win"]:
            n_approach += 1
        comparison.append(
            {
                "arm": name,
                "L": cfg["L"],
                "basal_mode": cfg.get("basal_mode"),
                "basal_q": cfg.get("basal_q", ""),
                "mad_kappa": cfg.get("mad_kappa", ""),
                "theta_delta_S": cfg.get("theta_delta_S", ""),
                "theta_R": cfg["theta_R"],
                "hop": cfg.get("hop", 1),
                "sens_sddb": s_sddb["sensitivity"],
                "n_detected_sddb": s_sddb["n_detected"],
                "n_sddb": s_sddb["n"],
                "sens_vfdb": s_vfdb["sensitivity"],
                "n_detected_vfdb": s_vfdb["n_detected"],
                "n_vfdb": s_vfdb["n"],
                "sens_all": s_all["sensitivity"],
                "n_detected_all": s_all["n_detected"],
                "n_all": s_all["n"],
                "median_lead_h": s_all["median_lead_h"],
                "far_per_24h": f["far_per_24h"],
                "total_episodes": f["total_episodes"],
                "total_search_hours": f["total_search_hours"],
                "n_controls": f["n_controls"],
                "far_over_opc": score["far_over_opc"],
                "structural_win": int(score["structural_win"]),
                "clear_advance_original": int(score["clear_advance_original"]),
                "secondary_interesting": int(score["secondary_interesting"]),
                "approaches": int(score["approaches"]),
                "score_reasons": "; ".join(score["reasons"]),
            }
        )

    stop_recommend = n_approach == 0 and not args.smoke

    summary: Dict[str, Any] = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "phase": "i0_structural_arms",
        "smoke": bool(args.smoke),
        "clinical_claim": False,
        "pre_registered": PRE_REGISTERED,
        "opc_anchors": {"sens_all": OPC_SENS_ALL, "far_per_24h": OPC_FAR},
        "arms": [{k: v for k, v in c.items() if k != "name"} | {"name": c["name"]}
                 for c in configs],
        "n_structural_wins": n_structural,
        "n_clear_advance_original": n_clear,
        "n_approaches": n_approach,
        "stop_surplus_primary_recommended": stop_recommend,
        "comparison": comparison,
        "sensitivity": {
            name: {
                "sddb": sens_block(sens_rows, name, "sddb"),
                "vfdb": sens_block(sens_rows, name, "vfdb"),
                "all_events": sens_block(sens_rows, name, None),
            }
            for name in det_names
        },
        "far": {name: far_block(far_rows, name) for name in det_names},
    }

    sens_path = RES / "i0_structural_arms_sens_per_record.csv"
    far_path = RES / "i0_structural_arms_nsrdb_far_per_record.csv"
    cmp_path = RES / "i0_structural_arms_comparison.csv"
    sum_path = RES / "i0_structural_arms_summary.json"

    write_csv(sens_path, sens_rows)
    write_csv(far_path, far_rows)
    write_csv(cmp_path, comparison)
    sum_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print("\n=== STRUCTURAL ARMS SUMMARY ===")
    print(
        f"OPC anchors: sens={OPC_SENS_ALL:.3f} FAR={OPC_FAR:.3f}  "
        f"structural bar: sens≥{SENS_STRUCTURAL} FAR≤{FAR_CAP_2X_OPC:.3f}"
    )
    for row in comparison:
        print(
            f"  {row['arm']}: sens={row['sens_all']:.3f} "
            f"({row['n_detected_all']}/{row['n_all']})  "
            f"FAR={row['far_per_24h']:.3f} ({row['far_over_opc']:.2f}×OPC)  "
            f"struct={row['structural_win']} clear={row['clear_advance_original']}  "
            f"| {row['score_reasons']}"
        )
    print(
        f"\nstructural_wins={n_structural}/{len(configs)}  "
        f"clear_advance={n_clear}  approaches={n_approach}  "
        f"stop_recommend={stop_recommend}"
    )
    print(f"\nWrote:\n  {sens_path}\n  {far_path}\n  {cmp_path}\n  {sum_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
