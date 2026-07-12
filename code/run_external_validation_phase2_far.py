#!/usr/bin/env python3
"""
External Validation Phase 2 — interim NSRDB control FAR expansion.

Reuses Phase-1 frozen primary path (analyze_control_record, false_alarm_rate,
FROZEN_* constants). Does NOT retune thresholds. Does NOT claim S5 met.

Outputs under results/:
  external_phase2_controls.csv
  external_phase2_far.json
  external_phase2_summary.json

Device-mismatch caveat is mandatory in every report artifact.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cctp_metrics_core import (
    FROZEN_HIGH_THRESHOLD,
    FROZEN_MIN_CONSECUTIVE,
    FROZEN_THETA3,
    FROZEN_Z_THRESHOLD,
    STRIDE,
    W_TAU,
    false_alarm_rate,
)
from run_external_validation_phase1 import (
    analyze_control_record,
    inventory_row_from_npz,
    load_npz,
    write_csv,
)

BASE = Path(__file__).resolve().parent.parent
DATA = BASE / "data"
RR_EXT = DATA / "rr_external"
RES = BASE / "results"
RES.mkdir(parents=True, exist_ok=True)

FROZEN_PARAMS = {
    "theta3": FROZEN_THETA3,
    "high_threshold": FROZEN_HIGH_THRESHOLD,
    "W_TAU": W_TAU,
    "stride": STRIDE,
    "z_threshold": FROZEN_Z_THRESHOLD,
    "min_consecutive": FROZEN_MIN_CONSECUTIVE,
    "alarm_rule": "abs_z_from_basal >= z_threshold sustained min_consecutive",
    "rr_clean_ms": [250.0, 2000.0],
}


def list_nsrdb_npz() -> List[Path]:
    return sorted(RR_EXT.glob("nsrdb_*_clean.npz"))


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Phase 2 interim FAR on all local NSRDB (frozen primary rule)"
    )
    ap.add_argument("--control-max-hours", type=float, default=12.0)
    ap.add_argument(
        "--z-threshold",
        type=float,
        default=FROZEN_Z_THRESHOLD,
        help="Must match frozen primary (default FROZEN_Z_THRESHOLD)",
    )
    ap.add_argument(
        "--min-consecutive",
        type=int,
        default=FROZEN_MIN_CONSECUTIVE,
        help="Must match frozen primary (default FROZEN_MIN_CONSECUTIVE)",
    )
    ap.add_argument(
        "--theta3",
        type=float,
        default=FROZEN_THETA3,
        help="Must match frozen primary (default FROZEN_THETA3)",
    )
    args = ap.parse_args()

    # Refuse silent retune away from frozen defaults
    if abs(args.z_threshold - FROZEN_Z_THRESHOLD) > 1e-12:
        print(f"ERROR: z_threshold={args.z_threshold} != frozen {FROZEN_Z_THRESHOLD}")
        return 2
    if args.min_consecutive != FROZEN_MIN_CONSECUTIVE:
        print(
            f"ERROR: min_consecutive={args.min_consecutive} != frozen {FROZEN_MIN_CONSECUTIVE}"
        )
        return 2
    if abs(args.theta3 - FROZEN_THETA3) > 1e-12:
        print(f"ERROR: theta3={args.theta3} != frozen {FROZEN_THETA3}")
        return 2

    npz_paths = list_nsrdb_npz()
    print(f"=== Phase 2 FAR (NSRDB controls): {len(npz_paths)} npz ===")
    print(f"Frozen params: {FROZEN_PARAMS}")

    inv_rows: List[dict] = []
    control_rows: List[dict] = []
    per_record: List[dict] = []

    for path in npz_paths:
        inv = inventory_row_from_npz(path)
        inv_rows.append(inv)
        if not inv["is_control"] or not inv["include"]:
            print(
                f"  skip {inv['record_id']}: include={inv['include']} "
                f"reason={inv.get('exclusion_reason')}"
            )
            continue
        d = load_npz(path)
        try:
            rows = analyze_control_record(
                d,
                z_threshold=args.z_threshold,
                min_consecutive=args.min_consecutive,
                theta3=args.theta3,
                max_hours=args.control_max_hours,
            )
            control_rows.extend(rows)
            # compact per-record summary for inventory updates
            by_m: Dict[str, dict] = {r["metric"]: r for r in rows}
            per_record.append(
                {
                    "record_id": inv["record_id"],
                    "total_hours_npz": inv["total_hours"],
                    "n_beats": inv["n_beats"],
                    "interp_frac": inv["interp_frac"],
                    "include": inv["include"],
                    "tau_s_episodes": int(by_m.get("tau_s", {}).get("n_episodes", 0) or 0),
                    "excess3_episodes": int(
                        by_m.get("excess3", {}).get("n_episodes", 0) or 0
                    ),
                    "total_hours_used": float(
                        by_m.get("tau_s", {}).get("total_hours_used", 0) or 0
                    ),
                    "search_hours_tau_s": float(
                        by_m.get("tau_s", {}).get("search_hours", 0) or 0
                    ),
                }
            )
            print(
                f"  control nsrdb/{inv['record_id']}: {len(rows)} metric rows "
                f"(tau_ep={per_record[-1]['tau_s_episodes']} "
                f"ex3_ep={per_record[-1]['excess3_episodes']} "
                f"h_used={per_record[-1]['total_hours_used']:.2f})"
            )
        except Exception as e:
            print(f"  control nsrdb/{inv['record_id']}: FAIL {type(e).__name__}: {e}")
            per_record.append(
                {
                    "record_id": inv["record_id"],
                    "total_hours_npz": inv["total_hours"],
                    "n_beats": inv["n_beats"],
                    "interp_frac": inv["interp_frac"],
                    "include": inv["include"],
                    "error": str(e),
                }
            )

    ctrl_path = RES / "external_phase2_controls.csv"
    write_csv(ctrl_path, control_rows)
    print(f"Wrote {ctrl_path} ({len(control_rows)} rows)")

    per_path = RES / "external_phase2_controls_per_record.csv"
    write_csv(per_path, per_record)
    print(f"Wrote {per_path} ({len(per_record)} rows)")

    far_report = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "phase": "external_validation_phase2_interim_far",
        "params": dict(FROZEN_PARAMS),
        "control_source": "PhysioNet NSRDB full local set (negative-control Holters; no VF)",
        "control_max_hours": args.control_max_hours,
        "n_npz_found": len(npz_paths),
        "n_controls_processed": len(
            {r["record_id"] for r in control_rows if r.get("metric") == "tau_s"}
        ),
        "by_metric": {},
        "device_mismatch": True,
        "device_mismatch_note": (
            "NSRDB is rhythm-healthy Holter ECG — NOT device-matched to VFDB/CU "
            "telemetry. Expanded n does not remove this mismatch; FAR remains an "
            "upper-bound / interim public estimate only."
        ),
        "s5_claim": False,
        "clinical_claim": False,
        "deployability_claim": False,
        "notes": [
            "FAR = total alarm episodes / total search hours * 24.",
            "Basal = early control window; search = remainder (capped at control_max_hours).",
            "Primary rule frozen: abs-z ≥ 2 sustained ≥ 3; θ₃=0.08; high=0.65; W_TAU=101.",
            "NSRDB is rhythm-healthy Holter — not device-matched to VFDB/CU telemetry.",
            "S5 (FAR ≤ 2/24h) is NOT claimed met solely from more healthy Holter controls.",
            "Phase 1 baseline (n=6): τ_s FAR ~34.4/24h, excess3 ~28.8/24h.",
        ],
        "phase1_baseline": {
            "n_controls": 6,
            "tau_s_far_per_24h": 34.400425526097,
            "excess3_far_per_24h": 28.800356254406793,
            "total_search_hours": 59.9992578125,
        },
    }
    for met in ("tau_s", "excess3"):
        sub = [r for r in control_rows if r["metric"] == met]
        far_report["by_metric"][met] = false_alarm_rate(sub)
        fr = far_report["by_metric"][met]
        print(
            f"  FAR {met}: n_ctrl={fr['n_controls']} hours={fr['total_search_hours']:.2f} "
            f"episodes={fr['total_episodes']} FAR/24h={fr['far_per_24h']} reason={fr['reason']}"
        )

    far_path = RES / "external_phase2_far.json"
    with far_path.open("w") as f:
        json.dump(far_report, f, indent=2, default=float)
    print(f"Wrote {far_path}")

    tau_fr = far_report["by_metric"].get("tau_s", {})
    ex_fr = far_report["by_metric"].get("excess3", {})
    summary = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "phase": "external_validation_phase2_interim_far",
        "status": "nsrdb_full_set_far_reported",
        "params_frozen": dict(FROZEN_PARAMS),
        "controls_csv": str(ctrl_path.relative_to(BASE)),
        "per_record_csv": str(per_path.relative_to(BASE)),
        "far_json": str(far_path.relative_to(BASE)),
        "n_controls": int(tau_fr.get("n_controls") or 0),
        "total_search_hours": float(tau_fr.get("total_search_hours") or 0.0),
        "far_by_metric": far_report["by_metric"],
        "phase1_vs_phase2": {
            "phase1_n": 6,
            "phase2_n": int(tau_fr.get("n_controls") or 0),
            "phase1_tau_s_far_per_24h": 34.400425526097,
            "phase2_tau_s_far_per_24h": float(tau_fr.get("far_per_24h") or float("nan")),
            "phase1_excess3_far_per_24h": 28.800356254406793,
            "phase2_excess3_far_per_24h": float(ex_fr.get("far_per_24h") or float("nan")),
        },
        "device_mismatch": True,
        "device_mismatch_note": far_report["device_mismatch_note"],
        "s5_claim": False,
        "s5_note": (
            "S5 (FAR ≤ 2/24h) is NOT met and is NOT claimed; expanding healthy NSRDB "
            "cannot alone establish clinical specificity or pass S5."
        ),
        "clinical_claim": False,
        "deployability_claim": False,
        "recommendation": (
            "Stop expanding public healthy Holter as the main specificity path. "
            "Next highest-value step: prepare institutional / device-matched non-event "
            "controls via docs/PHASE2_IRB_DATA_CHECKLIST.md (Tier A). Optional: keep "
            "full NSRDB FAR as public interim reference only."
        ),
        "next_step": (
            "Institutional pathway: complete partner identification + de-ID schema "
            "from PHASE2_IRB_DATA_CHECKLIST; do not retune frozen primary rule."
        ),
    }
    sum_path = RES / "external_phase2_summary.json"
    with sum_path.open("w") as f:
        json.dump(summary, f, indent=2, default=float)
    print(f"Wrote {sum_path}")
    print("=== Phase 2 interim FAR complete (no S5 / clinical claim) ===")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
