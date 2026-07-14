#!/usr/bin/env python3
"""
Pre-registered I0 surplus-persist parameter grid vs OPC companion.

Grid (fixed L=50, collapse_role=none — surplus-primary only):
  θ_ΔS ∈ {0.08, 0.10, 0.12, 0.15}
  θ_R^S ∈ {5, 8, 10, 12}

Protocol (identical family to OPSP Holter eval / ordinal bake-off):
  - Sensitivity: first alarm in post-basal pre-event window (SDDB + VFDB, n=33)
  - FAR: NSRDB Phase-2 controls, episode refractory 0.5 h, cap 12 h

Pre-registered clear-advance rule (written before any grid judgment):
  A cell is a **clear advance** vs OPC only if BOTH hold:
    (a) FAR/24h ≤ 2 × FAR_OPC   (cap ≈ 7.466… with FAR_OPC ≈ 3.733)
    (b) sens_all ≥ 0.65
  One-sided improvements, FAR still > 2×OPC with high sens, or marginal moves
  inside the high-FAR band are **not** success. I-confirm / free OR do not count.

OPC anchors (same protocol family):
  sens_all ≈ 0.424, FAR ≈ 3.733/24h

Usage:
  PYTHONPATH=code python3 code/run_i0_surplus_param_grid.py
  PYTHONPATH=code python3 code/run_i0_surplus_param_grid.py --smoke
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

# Reuse Holter event/control processing from the integrated eval runner
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
    BASE,
    RES,
)

# ---------------------------------------------------------------------------
# Pre-registered constants (locked before scoring any cell)
# ---------------------------------------------------------------------------

DEFAULT_THETA_DELTA_S: Tuple[float, ...] = (0.08, 0.10, 0.12, 0.15)
DEFAULT_THETA_R: Tuple[int, ...] = (5, 8, 10, 12)
FIXED_L = 50
FIXED_COLLAPSE_ROLE = "none"

OPC_SENS_ALL = float(BASELINE_ANCHORS["opc_L50"]["sens_all"])  # ≈ 0.424
OPC_FAR = float(BASELINE_ANCHORS["opc_L50"]["far_per_24h"])  # ≈ 3.733
FAR_MULTIPLIER_CAP = 2.0
SENS_MIN_CLEAR = 0.65

PRE_REGISTERED_RULES = {
    "name": "clear_advance_vs_opc",
    "opc_sens_all": OPC_SENS_ALL,
    "opc_far_per_24h": OPC_FAR,
    "require_far_le_multiplier_times_opc": FAR_MULTIPLIER_CAP,
    "far_cap_per_24h": FAR_MULTIPLIER_CAP * OPC_FAR,
    "require_sens_all_ge": SENS_MIN_CLEAR,
    "both_required": True,
    "non_success": [
        "FAR ≤ 2×OPC but sens_all in [OPC, 0.65)",
        "sens_all ≥ 0.65 but FAR > 2×OPC",
        "any marginal move inside the high-FAR band without both bars",
        "I-confirm / free OR / collapse-filter substituted for I0 primary",
    ],
    "judged_arm": "I0 surplus-persist only (collapse_role=none)",
}


# ---------------------------------------------------------------------------
# Pure helpers (unit-tested; no Holter I/O)
# ---------------------------------------------------------------------------


def build_i0_grid(
    theta_delta_S_values: Sequence[float] = DEFAULT_THETA_DELTA_S,
    theta_R_values: Sequence[int] = DEFAULT_THETA_R,
    *,
    L: int = FIXED_L,
) -> List[dict]:
    """Cartesian product θ_ΔS × θ_R with fixed L and collapse_role=none."""
    configs: List[dict] = []
    for th_ds in theta_delta_S_values:
        for th_r in theta_R_values:
            name = f"I0_dS{float(th_ds):.2f}_R{int(th_r)}"
            configs.append(
                {
                    "name": name,
                    "collapse_role": FIXED_COLLAPSE_ROLE,
                    "L": int(L),
                    "theta_delta_S": float(th_ds),
                    "theta_R": int(th_r),
                    "theta_D": 0.35,
                    "confirm_window": 5,
                    "modulate_delta_R": 2,
                    "modulate_delta_S": 0.02,
                }
            )
    return configs


def score_clear_advance(
    sens_all: float,
    far_per_24h: float,
    *,
    opc_sens_all: float = OPC_SENS_ALL,
    opc_far: float = OPC_FAR,
    far_multiplier: float = FAR_MULTIPLIER_CAP,
    sens_min: float = SENS_MIN_CLEAR,
) -> Dict[str, Any]:
    """
    Pre-registered clear-advance scorer.

    Pass only if FAR ≤ far_multiplier × opc_far AND sens_all ≥ sens_min.
    """
    far_cap = float(far_multiplier) * float(opc_far)
    sens = float(sens_all)
    far = float(far_per_24h)

    far_ok = bool(np.isfinite(far) and far <= far_cap + 1e-12)
    sens_ok = bool(np.isfinite(sens) and sens >= float(sens_min) - 1e-12)
    clear = bool(far_ok and sens_ok)

    reasons: List[str] = []
    if clear:
        reasons.append("meets both FAR ≤ 2×OPC and sens_all ≥ 0.65")
    else:
        if not far_ok:
            reasons.append(
                f"FAR={far:.4f} exceeds cap 2×OPC={far_cap:.4f}"
                if np.isfinite(far)
                else "FAR non-finite"
            )
        if not sens_ok:
            reasons.append(
                f"sens_all={sens:.4f} below min {sens_min}"
                if np.isfinite(sens)
                else "sens_all non-finite"
            )
        if far_ok and not sens_ok and np.isfinite(sens) and sens >= opc_sens_all:
            reasons.append(
                "one-sided: FAR ok but sens only ≥ OPC, not ≥ 0.65 (not success)"
            )
        if sens_ok and not far_ok:
            reasons.append(
                "one-sided: sens meets 0.65 but FAR still > 2×OPC (not success)"
            )

    return {
        "clear_advance": clear,
        "far_ok": far_ok,
        "sens_ok": sens_ok,
        "far_cap": far_cap,
        "sens_min": float(sens_min),
        "opc_sens_all": float(opc_sens_all),
        "opc_far": float(opc_far),
        "sens_all": sens,
        "far_per_24h": far,
        "delta_sens_vs_opc": sens - float(opc_sens_all) if np.isfinite(sens) else float("nan"),
        "far_ratio_vs_opc": far / float(opc_far) if opc_far > 0 and np.isfinite(far) else float("nan"),
        "reason": "; ".join(reasons),
    }


def score_grid_rows(
    rows: Sequence[Dict[str, Any]],
    *,
    opc_sens_all: float = OPC_SENS_ALL,
    opc_far: float = OPC_FAR,
) -> List[Dict[str, Any]]:
    """Attach clear_advance scoring to summary rows with sens_all and far_per_24h."""
    out: List[Dict[str, Any]] = []
    for r in rows:
        sc = score_clear_advance(
            float(r["sens_all"]),
            float(r["far_per_24h"]),
            opc_sens_all=opc_sens_all,
            opc_far=opc_far,
        )
        merged = dict(r)
        merged.update(
            {
                "clear_advance": sc["clear_advance"],
                "far_ok": sc["far_ok"],
                "sens_ok": sc["sens_ok"],
                "far_cap": sc["far_cap"],
                "delta_sens_vs_opc": sc["delta_sens_vs_opc"],
                "far_ratio_vs_opc": sc["far_ratio_vs_opc"],
                "score_reason": sc["reason"],
            }
        )
        out.append(merged)
    return out


# ---------------------------------------------------------------------------
# Main grid run
# ---------------------------------------------------------------------------


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(
        description="Pre-registered I0 surplus-persist 4×4 grid vs OPC"
    )
    ap.add_argument("--control-max-hours", type=float, default=CONTROL_MAX_HOURS)
    ap.add_argument("--basal-hours", type=float, default=BASAL_HOURS)
    ap.add_argument("--refractory-h", type=float, default=REFRACTORY_H)
    ap.add_argument(
        "--smoke",
        action="store_true",
        help="First 2 SDDB + 2 VFDB + 2 NSRDB only (wiring; not gating)",
    )
    ap.add_argument(
        "--scratch-dir",
        type=str,
        default="",
        help="Copy logs/artifacts here",
    )
    args = ap.parse_args(argv)

    configs = build_i0_grid()
    events = list_event_records()
    controls = list_nsrdb()
    if args.smoke:
        sddb = [e for e in events if e[0] == "sddb"][:2]
        vfdb = [e for e in events if e[0] == "vfdb"][:2]
        events = sddb + vfdb
        controls = controls[:2]
        print(f"SMOKE mode: {len(events)} events, {len(controls)} controls")

    print(f"Event records: {len(events)}  NSRDB controls: {len(controls)}")
    print(f"Grid cells: {len(configs)}  (L={FIXED_L}, collapse_role={FIXED_COLLAPSE_ROLE})")
    print("Pre-registered clear-advance:")
    print(
        f"  FAR ≤ {FAR_MULTIPLIER_CAP}×OPC ({FAR_MULTIPLIER_CAP * OPC_FAR:.4f}/24h) "
        f"AND sens_all ≥ {SENS_MIN_CLEAR}"
    )
    print(f"  OPC anchors: sens_all={OPC_SENS_ALL:.4f}, FAR={OPC_FAR:.4f}/24h")

    sens_rows: List[dict] = []
    for source, rec, path in events:
        print(f"  event {source}/{rec} ...", flush=True)
        rows = process_event_record(source, rec, path, configs)
        sens_rows.extend(rows)
        n_hit = sum(int(r.get("alarmed") or 0) for r in rows if not r.get("error"))
        print(f"    alarmed cells: {n_hit}/{len(rows)}", flush=True)

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

    # Per-cell aggregation
    grid_raw: List[dict] = []
    for cfg in configs:
        name = cfg["name"]
        s_all = sens_block(sens_rows, name, None)
        s_sddb = sens_block(sens_rows, name, "sddb")
        s_vfdb = sens_block(sens_rows, name, "vfdb")
        f = far_block(far_rows, name)
        grid_raw.append(
            {
                "detector": name,
                "L": cfg["L"],
                "theta_delta_S": cfg["theta_delta_S"],
                "theta_R": cfg["theta_R"],
                "collapse_role": cfg["collapse_role"],
                "sens_sddb": s_sddb["sensitivity"],
                "n_sddb": s_sddb["n"],
                "n_detected_sddb": s_sddb["n_detected"],
                "sens_vfdb": s_vfdb["sensitivity"],
                "n_vfdb": s_vfdb["n"],
                "n_detected_vfdb": s_vfdb["n_detected"],
                "sens_all": s_all["sensitivity"],
                "n_all": s_all["n"],
                "n_detected_all": s_all["n_detected"],
                "median_lead_h_all": s_all["median_lead_h"],
                "far_per_24h": f["far_per_24h"],
                "total_episodes": f["total_episodes"],
                "total_search_hours": f["total_search_hours"],
                "n_controls": f["n_controls"],
                "fraction_controls_alarmed": f["fraction_alarmed"],
            }
        )

    scored = score_grid_rows(grid_raw)
    n_clear = sum(1 for r in scored if r["clear_advance"])

    summary: Dict[str, Any] = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "phase": "i0_surplus_param_grid",
        "smoke": bool(args.smoke),
        "clinical_claim": False,
        "pre_registered_rules": PRE_REGISTERED_RULES,
        "grid_spec": {
            "theta_delta_S": list(DEFAULT_THETA_DELTA_S),
            "theta_R": list(DEFAULT_THETA_R),
            "L": FIXED_L,
            "collapse_role": FIXED_COLLAPSE_ROLE,
            "n_cells": len(configs),
        },
        "params": {
            "control_max_hours": args.control_max_hours,
            "basal_hours": args.basal_hours,
            "refractory_h": args.refractory_h,
            "n_events": len(events),
            "n_controls": len(controls),
        },
        "baseline_anchors": {
            "opc_L50": BASELINE_ANCHORS["opc_L50"],
            "absz_tau_s": BASELINE_ANCHORS["absz_tau_s"],
        },
        "n_clear_advance": n_clear,
        "any_clear_advance": n_clear > 0,
        "grid": scored,
        "strategic_default_if_zero": (
            "no clear advance under this proxy/protocol; "
            "prefer OPC as primary ordinal specificity arm or structural surplus change"
        ),
    }

    # Paths
    sens_path = RES / "i0_surplus_param_grid_sens_per_record.csv"
    far_path = RES / "i0_surplus_param_grid_nsrdb_far_per_record.csv"
    grid_path = RES / "i0_surplus_param_grid.csv"
    sum_path = RES / "i0_surplus_param_grid_summary.json"

    write_csv(sens_path, sens_rows)
    write_csv(far_path, far_rows)
    write_csv(grid_path, scored)
    sum_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print("\n=== PRE-REGISTERED RULES (locked before judgment) ===")
    print(json.dumps(PRE_REGISTERED_RULES, indent=2))
    print("\n=== GRID RESULTS ===")
    print(
        f"{'cell':<22} {'sens_all':>8} {'det/n':>8} {'FAR/24h':>10} "
        f"{'FAR/OPC':>8} {'clear?':>7}  reason"
    )
    for r in scored:
        print(
            f"{r['detector']:<22} {r['sens_all']:8.3f} "
            f"{r['n_detected_all']}/{r['n_all']:<5} {r['far_per_24h']:10.3f} "
            f"{r['far_ratio_vs_opc']:8.2f} "
            f"{'YES' if r['clear_advance'] else 'no':>7}  {r['score_reason']}"
        )
    print(f"\nClear-advance count: {n_clear} / {len(scored)}")
    print(f"\nWrote:\n  {sens_path}\n  {far_path}\n  {grid_path}\n  {sum_path}")

    if args.scratch_dir:
        scratch = Path(args.scratch_dir)
        scratch.mkdir(parents=True, exist_ok=True)
        for p in (sens_path, far_path, grid_path, sum_path):
            (scratch / p.name).write_bytes(p.read_bytes())
        log = scratch / "i0_grid_run.log"
        lines = [
            f"generated_at={summary['generated_at']}",
            f"smoke={args.smoke}",
            f"n_events={len(events)} n_controls={len(controls)} n_cells={len(configs)}",
            f"pre_registered: FAR≤{FAR_MULTIPLIER_CAP}×OPC={FAR_MULTIPLIER_CAP * OPC_FAR:.4f} "
            f"AND sens_all≥{SENS_MIN_CLEAR}",
            f"opc anchors: sens={OPC_SENS_ALL:.4f} far={OPC_FAR:.4f}",
            f"n_clear_advance={n_clear}",
            "",
        ]
        for r in scored:
            lines.append(
                f"{r['detector']}: sens_all={r['sens_all']:.4f} "
                f"({r['n_detected_all']}/{r['n_all']}) "
                f"FAR={r['far_per_24h']:.4f} ratio={r['far_ratio_vs_opc']:.3f} "
                f"clear={r['clear_advance']} | {r['score_reason']}"
            )
        log.write_text("\n".join(lines) + "\n", encoding="utf-8")
        print(f"Scratch copies → {scratch}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
