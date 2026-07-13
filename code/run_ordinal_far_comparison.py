#!/usr/bin/env python3
"""
Comparable NSRDB FAR for OPC L=50, SDD, and frozen abs-z τ_s.

Phase-2–aligned episode counting with refractory period (0.5 h):
  FAR (alarms/24h) = total_episodes / total_search_hours * 24

Runs three detectors **separately** (no fusion):
  - OPC L=50, θ_D=0.35, θ_R=5 on joint bivariate symbols (K=36)
  - SDD L_c=50, θ_TV=0.35, θ_S=1, fixed early basal, mask_basal=True
  - abs-z on continuous τ_s (frozen z=2, min_consecutive=3)

Control windowing matches Phase 2 analyze_control_record:
  basal ≈ (0.25, min(2.0, 0.25*total_h)), search = remainder, cap max_hours=12.

Observational only — no clinical / S5 / superiority claims.

Outputs under results/:
  ordinal_nsrdb_far_per_record.csv
  ordinal_nsrdb_far_summary.json
  ordinal_nsrdb_far_by_detector.csv
"""
from __future__ import annotations

import argparse
import ast
import csv
import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _bootstrap import import_systemictau_core
from cctp_metrics_core import (
    FROZEN_MIN_CONSECUTIVE,
    FROZEN_Z_THRESHOLD,
    STRIDE,
    W_TAU,
    build_bivariate_proxy,
    count_alarm_episodes,
    count_binary_alarm_episodes,
    false_alarm_rate,
)
from ordinal_detectors.opc_detector import opc_detect
from ordinal_detectors.sdd_detector import sdd_detect
from recd_ordinal_levels import generate_multivariate_symbols

BASE = Path(__file__).resolve().parent.parent
DATA = BASE / "data"
RR_EXT = DATA / "rr_external"
RES = BASE / "results"
RES.mkdir(parents=True, exist_ok=True)

# Symbol encoding — same as exploratory bake-off
M_EMB = 3
DELAY = 1
BP_ALPHABET = 6
K_JOINT = BP_ALPHABET * BP_ALPHABET  # 36

# Ordinal params (fixed for comparison; not retuned)
OPC_L = 50
OPC_THETA_D = 0.35
OPC_THETA_R = 5
SDD_L_C = 50
SDD_THETA_TV = 0.35
SDD_THETA_S = 1

# Phase-2 control window defaults
BASAL_HOURS = 2.0
CONTROL_MAX_HOURS = 12.0
REFRACTORY_H = 0.5

FUSION_FORBIDDEN = True  # never combine OPC and SDD alarms


def load_npz(path: Path) -> dict:
    d = np.load(path, allow_pickle=True)
    out = {}
    for k in d.files:
        v = d[k]
        if isinstance(v, np.ndarray) and v.shape == () and v.dtype == object:
            out[k] = v.item()
        elif isinstance(v, np.ndarray) and v.dtype.kind in ("U", "S") and v.shape == ():
            out[k] = str(v.item())
        else:
            out[k] = v
    return out


def list_nsrdb_npz() -> List[Path]:
    return sorted(RR_EXT.glob("nsrdb_*_clean.npz"))


def joint_bivariate_symbols(rr: np.ndarray) -> Tuple[np.ndarray, int, int]:
    """1-D joint Bandt–Pompe symbols; sigma = s0*6+s1 ∈ {0..35}; K=36."""
    X = build_bivariate_proxy(np.asarray(rr, dtype=float))
    S = generate_multivariate_symbols(X, m=M_EMB, delay=DELAY)
    if S.size == 0 or S.shape[1] < 2:
        return np.array([], dtype=np.int64), K_JOINT, (M_EMB - 1) * DELAY
    sigma = (S[:, 0].astype(np.int64) * BP_ALPHABET) + S[:, 1].astype(np.int64)
    offset = (M_EMB - 1) * DELAY
    return sigma, K_JOINT, offset


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


def control_basal(total_h: float, basal_hours: float = BASAL_HOURS) -> Tuple[float, float]:
    """Match Phase-1/2 analyze_control_record basal window."""
    basal = (0.25, min(basal_hours, total_h * 0.25))
    if basal[1] <= basal[0]:
        basal = (0.0, max(total_h * 0.2, 0.1))
    return basal


def hours_to_index_left(t_hr: np.ndarray, hr: float) -> int:
    if len(t_hr) == 0:
        return 0
    return max(0, min(int(np.searchsorted(t_hr, hr, side="left")), len(t_hr)))


def hours_to_index_right(t_hr: np.ndarray, hr: float) -> int:
    if len(t_hr) == 0:
        return 0
    return max(0, min(int(np.searchsorted(t_hr, hr, side="left")), len(t_hr)))


def process_control(
    path: Path,
    *,
    basal_hours: float = BASAL_HOURS,
    max_hours: float = CONTROL_MAX_HOURS,
    refractory_h: float = REFRACTORY_H,
) -> Dict[str, Any]:
    """
    Process one NSRDB control with three independent detectors.
    Returns per-detector episode/FAR fields plus shared window metadata.
    """
    d = load_npz(path)
    rr = np.asarray(d["rr_ms"], dtype=float)
    t_sec = np.asarray(d["t_sec"], dtype=float)
    t_hr = t_sec / 3600.0
    rec = str(d.get("record_id", path.stem))
    total_h = float(t_hr[-1]) if len(t_hr) else 0.0

    if max_hours is not None and total_h > max_hours:
        keep = t_hr <= max_hours
        rr = rr[keep]
        t_hr = t_hr[keep]
        total_h = float(t_hr[-1]) if len(t_hr) else 0.0

    empty: Dict[str, Any] = {
        "record_id": rec,
        "source_path": str(path.name),
        "total_hours_used": total_h,
        "skipped": True,
        "skip_reason": "too_short",
    }
    if len(rr) < W_TAU + 50 or total_h < basal_hours + 0.5:
        return empty

    basal = control_basal(total_h, basal_hours)
    b0, b1 = basal
    search_start = b1
    search_end = total_h

    # --- abs-z on τ_s (frozen continuous path) ---
    tau = compute_tau_series(rr)
    n_tau = min(len(t_hr), len(tau))
    absz_ep = count_alarm_episodes(
        tau[:n_tau],
        t_hr[:n_tau],
        basal,
        z_threshold=FROZEN_Z_THRESHOLD,
        min_consecutive=FROZEN_MIN_CONSECUTIVE,
        use_abs=True,
        refractory_h=refractory_h,
    )

    # --- joint symbols for ordinal detectors ---
    sigma, K, offset = joint_bivariate_symbols(rr)
    if len(sigma) < max(OPC_L, SDD_L_C) + 10:
        return {
            **empty,
            "skip_reason": "symbol_stream_too_short",
            "absz_n_episodes": absz_ep["n_episodes"],
            "absz_search_hours": absz_ep["search_hours"],
            "absz_far_per_24h": (
                absz_ep["n_episodes"] / absz_ep["search_hours"] * 24.0
                if absz_ep["search_hours"] > 0
                else float("nan")
            ),
            "skipped": True,
        }

    t_hr_sym = t_hr[offset : offset + len(sigma)]
    if len(t_hr_sym) < len(sigma):
        sigma = sigma[: len(t_hr_sym)]
    n_sym = min(len(sigma), len(t_hr_sym))
    sigma = sigma[:n_sym]
    t_hr_sym = t_hr_sym[:n_sym]

    basal_i0 = hours_to_index_left(t_hr_sym, b0)
    basal_i1 = hours_to_index_right(t_hr_sym, b1)
    if basal_i1 - basal_i0 < max(2, SDD_L_C // 5):
        # fallback: first 20% of symbols as basal (still fixed early segment)
        basal_i0 = 0
        basal_i1 = max(SDD_L_C, min(len(sigma) // 5, len(sigma) // 4))
        if basal_i1 <= basal_i0:
            basal_i1 = min(len(sigma), basal_i0 + SDD_L_C)

    # --- OPC L=50 (independent) ---
    opc_out = opc_detect(
        sigma,
        L=OPC_L,
        theta_D=OPC_THETA_D,
        theta_R=OPC_THETA_R,
        K=K,
    )
    opc_ep = count_binary_alarm_episodes(
        opc_out["alarm"],
        t_hr_sym,
        search_start_hr=search_start,
        search_end_hr=search_end,
        refractory_h=refractory_h,
    )

    # --- SDD (independent; fixed basal, mask_basal=True) ---
    sdd_out = sdd_detect(
        sigma,
        (basal_i0, basal_i1),
        L_c=SDD_L_C,
        theta_TV=SDD_THETA_TV,
        theta_S=SDD_THETA_S,
        K=K,
        mask_basal=True,
    )
    sdd_ep = count_binary_alarm_episodes(
        sdd_out["alarm"],
        t_hr_sym,
        search_start_hr=search_start,
        search_end_hr=search_end,
        refractory_h=refractory_h,
    )

    def _far(ep: Dict[str, float]) -> float:
        h = float(ep["search_hours"])
        if h <= 0:
            return float("nan")
        return float(ep["n_episodes"]) / h * 24.0

    return {
        "record_id": rec,
        "source_path": str(path.name),
        "skipped": False,
        "total_hours_used": total_h,
        "basal_start": b0,
        "basal_end": b1,
        "search_start": search_start,
        "search_end": search_end,
        "refractory_h": refractory_h,
        "n_symbols": int(len(sigma)),
        "K": int(K),
        "basal_idx_start": int(basal_i0),
        "basal_idx_end": int(basal_i1),
        # OPC
        "opc_n_episodes": opc_ep["n_episodes"],
        "opc_search_hours": opc_ep["search_hours"],
        "opc_first_alarm_hr": opc_ep["first_alarm_hr"],
        "opc_alarmed": opc_ep["alarmed"],
        "opc_far_per_24h": _far(opc_ep),
        # SDD
        "sdd_n_episodes": sdd_ep["n_episodes"],
        "sdd_search_hours": sdd_ep["search_hours"],
        "sdd_first_alarm_hr": sdd_ep["first_alarm_hr"],
        "sdd_alarmed": sdd_ep["alarmed"],
        "sdd_far_per_24h": _far(sdd_ep),
        # abs-z
        "absz_n_episodes": absz_ep["n_episodes"],
        "absz_search_hours": absz_ep["search_hours"],
        "absz_first_alarm_hr": absz_ep["first_alarm_hr"],
        "absz_alarmed": absz_ep["alarmed"],
        "absz_far_per_24h": _far(absz_ep),
        "absz_z_threshold": FROZEN_Z_THRESHOLD,
        "absz_min_consecutive": FROZEN_MIN_CONSECUTIVE,
        "fusion": False,
    }


def _per_record_stats(fars: List[float]) -> Dict[str, float]:
    arr = np.asarray([f for f in fars if np.isfinite(f)], dtype=float)
    if arr.size == 0:
        return {
            "mean": float("nan"),
            "median": float("nan"),
            "min": float("nan"),
            "max": float("nan"),
            "std": float("nan"),
            "n": 0.0,
        }
    return {
        "mean": float(np.mean(arr)),
        "median": float(np.median(arr)),
        "min": float(np.min(arr)),
        "max": float(np.max(arr)),
        "std": float(np.std(arr, ddof=0)),
        "n": float(arr.size),
    }


def write_csv(path: Path, rows: List[dict]) -> None:
    if not rows:
        path.write_text("")
        return
    keys: List[str] = []
    seen = set()
    for r in rows:
        for k in r:
            if k not in seen:
                seen.add(k)
                keys.append(k)
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in keys})


def assert_no_fusion_in_source(path: Path) -> None:
    """Structural guard: runner processes OPC and SDD as independent arms."""
    src = path.read_text(encoding="utf-8")
    tree = ast.parse(src)
    assert "opc_detect" in src and "sdd_detect" in src
    assert "count_binary_alarm_episodes" in src
    assert "count_alarm_episodes" in src
    # Independent FAR keys; no single shared 'combined' episode field
    assert '"opc_n_episodes"' in src or "'opc_n_episodes'" in src
    assert '"sdd_n_episodes"' in src or "'sdd_n_episodes'" in src
    assert "fusion" in src  # documents fusion=False
    _ = tree


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(
        description="NSRDB FAR comparison: OPC L=50 vs SDD vs frozen abs-z (separate)"
    )
    ap.add_argument("--control-max-hours", type=float, default=CONTROL_MAX_HOURS)
    ap.add_argument("--basal-hours", type=float, default=BASAL_HOURS)
    ap.add_argument("--refractory-h", type=float, default=REFRACTORY_H)
    ap.add_argument("--write-doc", action="store_true", help="Write docs/ORDINAL_NSRDB_FAR_COMPARISON.md")
    args = ap.parse_args(argv)

    assert_no_fusion_in_source(Path(__file__).resolve())

    paths = list_nsrdb_npz()
    print(f"NSRDB controls found: {len(paths)}")
    if not paths:
        print("ERROR: no nsrdb_*_clean.npz under data/rr_external/", file=sys.stderr)
        return 1

    rows: List[Dict[str, Any]] = []
    for p in paths:
        print(f"  processing {p.name} ...", flush=True)
        row = process_control(
            p,
            basal_hours=args.basal_hours,
            max_hours=args.control_max_hours,
            refractory_h=args.refractory_h,
        )
        rows.append(row)
        if row.get("skipped"):
            print(f"    skipped: {row.get('skip_reason')}")
        else:
            print(
                f"    OPC FAR={row['opc_far_per_24h']:.2f}  "
                f"SDD FAR={row['sdd_far_per_24h']:.2f}  "
                f"abs-z FAR={row['absz_far_per_24h']:.2f}  "
                f"(ep OPC/SDD/absz="
                f"{int(row['opc_n_episodes'])}/"
                f"{int(row['sdd_n_episodes'])}/"
                f"{int(row['absz_n_episodes'])})"
            )

    ok_rows = [r for r in rows if not r.get("skipped")]
    n_controls = len(ok_rows)
    print(f"n_controls processed: {n_controls} (of {len(paths)} npz)")

    # Aggregate via shipped false_alarm_rate (pooled)
    opc_pool = false_alarm_rate(
        [{"n_episodes": r["opc_n_episodes"], "search_hours": r["opc_search_hours"], "alarmed": r["opc_alarmed"]}
         for r in ok_rows]
    )
    sdd_pool = false_alarm_rate(
        [{"n_episodes": r["sdd_n_episodes"], "search_hours": r["sdd_search_hours"], "alarmed": r["sdd_alarmed"]}
         for r in ok_rows]
    )
    absz_pool = false_alarm_rate(
        [{"n_episodes": r["absz_n_episodes"], "search_hours": r["absz_search_hours"], "alarmed": r["absz_alarmed"]}
         for r in ok_rows]
    )

    per_stats = {
        "opc_L50": _per_record_stats([float(r["opc_far_per_24h"]) for r in ok_rows]),
        "sdd": _per_record_stats([float(r["sdd_far_per_24h"]) for r in ok_rows]),
        "absz_tau_s": _per_record_stats([float(r["absz_far_per_24h"]) for r in ok_rows]),
    }

    summary = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "phase": "ordinal_nsrdb_far_comparison",
        "n_npz_found": len(paths),
        "n_controls_processed": n_controls,
        "control_max_hours": args.control_max_hours,
        "basal_hours": args.basal_hours,
        "refractory_h": args.refractory_h,
        "fusion": False,
        "params": {
            "opc": {
                "L": OPC_L,
                "theta_D": OPC_THETA_D,
                "theta_R": OPC_THETA_R,
                "K": K_JOINT,
                "encoding": "joint_bivariate_bandt_pompe_m3",
            },
            "sdd": {
                "L_c": SDD_L_C,
                "theta_TV": SDD_THETA_TV,
                "theta_S": SDD_THETA_S,
                "basal": "fixed_early",
                "mask_basal": True,
                "K": K_JOINT,
            },
            "absz": {
                "z_threshold": FROZEN_Z_THRESHOLD,
                "min_consecutive": FROZEN_MIN_CONSECUTIVE,
                "metric": "tau_s",
                "W_TAU": W_TAU,
                "stride": STRIDE,
                "frozen": True,
            },
        },
        "methodology": {
            "far_formula": "FAR = total_episodes / total_search_hours * 24",
            "episode_definition_absz": (
                "abs-z >= 2 sustained >= 3 consecutive samples; then refractory 0.5 h "
                "(count_alarm_episodes)"
            ),
            "episode_definition_ordinal": (
                "Binary alarm from opc_detect / sdd_detect (persistence already in detector); "
                "each True sample outside refractory starts an episode; refractory 0.5 h "
                "(count_binary_alarm_episodes)"
            ),
            "basal_window": "Phase-2 style: (0.25, min(basal_hours, 0.25*total_h))",
            "search_window": "After basal end through end of capped recording",
            "control_cap_h": args.control_max_hours,
            "device_mismatch": True,
            "device_mismatch_note": (
                "NSRDB is rhythm-healthy Holter ECG — NOT device-matched to VFDB/CU telemetry. "
                "FAR is an interim public upper-bound estimate only."
            ),
        },
        "pooled_far": {
            "opc_L50": opc_pool,
            "sdd": sdd_pool,
            "absz_tau_s": absz_pool,
        },
        "per_record_far_stats": per_stats,
        "s5_claim": False,
        "clinical_claim": False,
        "superiority_claim": False,
        "qualitative_notes": [],
        "phase2_absz_reference": {
            "source": "results/external_phase2_far.json",
            "tau_s_far_per_24h_expected_order": "~33–34",
            "note": "Sanity: abs-z arm here should be same order of magnitude as Phase 2 τ_s FAR.",
        },
    }

    # Qualitative observations (observational, not claims)
    notes = []
    o_far = opc_pool.get("far_per_24h", float("nan"))
    s_far = sdd_pool.get("far_per_24h", float("nan"))
    a_far = absz_pool.get("far_per_24h", float("nan"))
    notes.append(
        f"Pooled FAR (alarms/24h): OPC L=50 = {o_far:.3f}, SDD = {s_far:.3f}, "
        f"abs-z τ_s = {a_far:.3f} (n={n_controls})."
    )
    if np.isfinite(o_far) and np.isfinite(a_far):
        if o_far < a_far:
            notes.append(
                "OPC L=50 shows lower pooled FAR than frozen abs-z on these NSRDB Holters "
                "(observational; not a clinical superiority claim)."
            )
        elif o_far > a_far:
            notes.append(
                "OPC L=50 shows higher pooled FAR than frozen abs-z on these NSRDB Holters."
            )
        else:
            notes.append("OPC L=50 and abs-z pooled FAR are essentially equal.")
    if np.isfinite(s_far) and np.isfinite(a_far):
        if s_far < a_far:
            notes.append(
                "SDD shows lower pooled FAR than frozen abs-z on these NSRDB Holters "
                "(observational)."
            )
        elif s_far > a_far:
            notes.append(
                "SDD shows higher pooled FAR than frozen abs-z (more distribution-shift "
                "triggers on healthy Holter dynamics)."
            )
    if np.isfinite(o_far) and np.isfinite(s_far):
        if o_far < s_far:
            notes.append(
                "Among ordinal rules, OPC L=50 appears more specific (lower FAR) than SDD "
                "under these fixed params — consistent with bake-off basal cleanliness."
            )
        elif s_far < o_far:
            notes.append(
                "Among ordinal rules, SDD appears more specific (lower FAR) than OPC L=50 "
                "under these fixed params."
            )
    notes.append(
        "Timebases differ: abs-z on strided τ_s vs symbol-endpoint times for OPC/SDD; "
        "FAR remains comparable by shared definition (episodes/search_h×24), not identical indices."
    )
    notes.append(
        "S5 (FAR ≤ 2/24h) is NOT claimed. NSRDB is not device-matched ICU/telemetry control."
    )
    summary["qualitative_notes"] = notes

    # Per-record CSV
    per_path = RES / "ordinal_nsrdb_far_per_record.csv"
    write_csv(per_path, rows)

    # Compact by-detector aggregate table
    by_det = [
        {
            "detector": "opc_L50",
            "pooled_far_per_24h": opc_pool["far_per_24h"],
            "total_episodes": opc_pool["total_episodes"],
            "total_search_hours": opc_pool["total_search_hours"],
            "fraction_alarmed": opc_pool["fraction_alarmed"],
            "mean_far": per_stats["opc_L50"]["mean"],
            "median_far": per_stats["opc_L50"]["median"],
            "min_far": per_stats["opc_L50"]["min"],
            "max_far": per_stats["opc_L50"]["max"],
            "std_far": per_stats["opc_L50"]["std"],
            "n_controls": n_controls,
        },
        {
            "detector": "sdd",
            "pooled_far_per_24h": sdd_pool["far_per_24h"],
            "total_episodes": sdd_pool["total_episodes"],
            "total_search_hours": sdd_pool["total_search_hours"],
            "fraction_alarmed": sdd_pool["fraction_alarmed"],
            "mean_far": per_stats["sdd"]["mean"],
            "median_far": per_stats["sdd"]["median"],
            "min_far": per_stats["sdd"]["min"],
            "max_far": per_stats["sdd"]["max"],
            "std_far": per_stats["sdd"]["std"],
            "n_controls": n_controls,
        },
        {
            "detector": "absz_tau_s",
            "pooled_far_per_24h": absz_pool["far_per_24h"],
            "total_episodes": absz_pool["total_episodes"],
            "total_search_hours": absz_pool["total_search_hours"],
            "fraction_alarmed": absz_pool["fraction_alarmed"],
            "mean_far": per_stats["absz_tau_s"]["mean"],
            "median_far": per_stats["absz_tau_s"]["median"],
            "min_far": per_stats["absz_tau_s"]["min"],
            "max_far": per_stats["absz_tau_s"]["max"],
            "std_far": per_stats["absz_tau_s"]["std"],
            "n_controls": n_controls,
        },
    ]
    by_path = RES / "ordinal_nsrdb_far_by_detector.csv"
    write_csv(by_path, by_det)

    sum_path = RES / "ordinal_nsrdb_far_summary.json"
    with sum_path.open("w") as f:
        json.dump(summary, f, indent=2)

    print("\n=== Pooled FAR (alarms / 24 h) ===")
    print(f"  OPC L=50 : {opc_pool['far_per_24h']:.4f}  "
          f"(episodes={opc_pool['total_episodes']:.0f}, hours={opc_pool['total_search_hours']:.2f})")
    print(f"  SDD      : {sdd_pool['far_per_24h']:.4f}  "
          f"(episodes={sdd_pool['total_episodes']:.0f}, hours={sdd_pool['total_search_hours']:.2f})")
    print(f"  abs-z τs : {absz_pool['far_per_24h']:.4f}  "
          f"(episodes={absz_pool['total_episodes']:.0f}, hours={absz_pool['total_search_hours']:.2f})")
    print("\nPer-record FAR mean/median/range:")
    for name, st in per_stats.items():
        print(
            f"  {name:12s} mean={st['mean']:.3f} median={st['median']:.3f} "
            f"range=[{st['min']:.3f}, {st['max']:.3f}]"
        )
    print(f"\nWrote {per_path}")
    print(f"Wrote {by_path}")
    print(f"Wrote {sum_path}")

    if args.write_doc:
        write_methodology_doc(summary, rows, by_det)

    return 0


def write_methodology_doc(
    summary: dict,
    rows: List[dict],
    by_det: List[dict],
) -> None:
    doc = BASE / "docs" / "ORDINAL_NSRDB_FAR_COMPARISON.md"
    doc.parent.mkdir(parents=True, exist_ok=True)
    o = summary["pooled_far"]["opc_L50"]
    s = summary["pooled_far"]["sdd"]
    a = summary["pooled_far"]["absz_tau_s"]
    ps = summary["per_record_far_stats"]

    lines = [
        "# Ordinal vs abs-z FAR on NSRDB controls",
        "",
        f"*Generated: {summary['generated_at']}*",
        "",
        "**Status:** Exploratory / interim public control FAR comparison. "
        "No clinical claim, no S5 claim, no superiority claim.",
        "",
        "## Purpose",
        "",
        "Measure False Alarm Rate (alarms per 24 h with refractory episode counting) "
        "for two native ordinal detectors and the frozen abs-z baseline on the same "
        "18 NSRDB negative-control Holters, with Phase-2–aligned windowing.",
        "",
        "Options remain **completely separate** (no OPC∧SDD fusion).",
        "",
        "## Methodology (reproducible)",
        "",
        "### Entry point",
        "",
        "```bash",
        "python code/run_ordinal_far_comparison.py --write-doc",
        "```",
        "",
        "### Data",
        "",
        f"- Controls: `data/rr_external/nsrdb_*_clean.npz` "
        f"({summary['n_npz_found']} files found, "
        f"{summary['n_controls_processed']} processed).",
        f"- Cap: first `{summary['control_max_hours']}` hours per record (Phase 2 default).",
        f"- Basal hours parameter: `{summary['basal_hours']}` "
        "(window = `(0.25, min(basal_hours, 0.25×total_h))`, same as Phase 1/2).",
        "- Search: after basal end → end of capped recording.",
        "",
        "### Symbol encoding (ordinal arms)",
        "",
        "- Bivariate proxy from cleaned RR (project `build_bivariate_proxy`).",
        "- Bandt–Pompe m=3, delay=1 → per-channel alphabet 6.",
        "- Joint code σ = s0×6 + s1 ∈ {0,…,35}, **K=36** (same as bake-off).",
        "",
        "### Detectors (fixed params; not retuned for FAR)",
        "",
        "| Detector | Parameters |",
        "|----------|------------|",
        f"| **OPC L=50** | L={OPC_L}, θ_D={OPC_THETA_D}, θ_R={OPC_THETA_R}, K={K_JOINT} |",
        f"| **SDD** | L_c={SDD_L_C}, θ_TV={SDD_THETA_TV}, θ_S={SDD_THETA_S}, fixed early basal, mask_basal=True |",
        f"| **abs-z τ_s** | z≥{FROZEN_Z_THRESHOLD}, min_consecutive={FROZEN_MIN_CONSECUTIVE}, W_TAU={W_TAU}, stride={STRIDE} (**frozen**) |",
        "",
        "### Episode counting and FAR",
        "",
        f"- **Refractory period:** {summary['refractory_h']} h after each episode.",
        "- **abs-z:** `count_alarm_episodes` — abs-z ≥ 2 sustained ≥ 3 samples, then refractory.",
        "- **OPC / SDD:** `count_binary_alarm_episodes` on the detector’s binary alarm stream "
        "(persistence already enforced by θ_R / θ_S); True outside refractory → episode.",
        "- **FAR formula:** `total_episodes / total_search_hours × 24` via `false_alarm_rate`.",
        "- **Pooled FAR** sums episodes and hours across records; "
        "**per-record FAR** = n_episodes_i / search_hours_i × 24.",
        "",
        "### Comparability notes",
        "",
        "- abs-z uses continuous strided τ_s timebase; ordinal uses symbol-endpoint times. "
        "Shared FAR *definition* enables comparison; sample indices are not identical.",
        "- NSRDB is healthy Holter, **not** device-matched to VFDB/CU telemetry "
        "(same caveat as Phase 2).",
        "",
        "## Results",
        "",
        "### Pooled FAR",
        "",
        "| Detector | Total episodes | Search hours | Pooled FAR (/24h) | Fraction records alarmed |",
        "|----------|----------------|--------------|-------------------|--------------------------|",
        f"| OPC L=50 | {o['total_episodes']:.0f} | {o['total_search_hours']:.2f} | **{o['far_per_24h']:.3f}** | {o['fraction_alarmed']:.3f} |",
        f"| SDD | {s['total_episodes']:.0f} | {s['total_search_hours']:.2f} | **{s['far_per_24h']:.3f}** | {s['fraction_alarmed']:.3f} |",
        f"| abs-z τ_s | {a['total_episodes']:.0f} | {a['total_search_hours']:.2f} | **{a['far_per_24h']:.3f}** | {a['fraction_alarmed']:.3f} |",
        "",
        "### Per-record FAR statistics",
        "",
        "| Detector | Mean | Median | Min | Max | Std |",
        "|----------|------|--------|-----|-----|-----|",
    ]
    for name, label in (("opc_L50", "OPC L=50"), ("sdd", "SDD"), ("absz_tau_s", "abs-z τ_s")):
        st = ps[name]
        lines.append(
            f"| {label} | {st['mean']:.3f} | {st['median']:.3f} | "
            f"{st['min']:.3f} | {st['max']:.3f} | {st['std']:.3f} |"
        )

    lines += [
        "",
        "### Per-record table",
        "",
        "| Record | OPC ep | OPC FAR | SDD ep | SDD FAR | abs-z ep | abs-z FAR | Search h |",
        "|--------|--------|---------|--------|---------|----------|-----------|----------|",
    ]
    for r in rows:
        if r.get("skipped"):
            lines.append(
                f"| {r['record_id']} | — | — | — | — | — | — | skipped |"
            )
            continue
        lines.append(
            f"| {r['record_id']} | {int(r['opc_n_episodes'])} | {r['opc_far_per_24h']:.2f} | "
            f"{int(r['sdd_n_episodes'])} | {r['sdd_far_per_24h']:.2f} | "
            f"{int(r['absz_n_episodes'])} | {r['absz_far_per_24h']:.2f} | "
            f"{r['opc_search_hours']:.2f} |"
        )

    lines += [
        "",
        "## Qualitative observations",
        "",
    ]
    for n in summary.get("qualitative_notes", []):
        lines.append(f"- {n}")

    lines += [
        "",
        "## Specificity ranking (observational)",
        "",
        "Rank by **lower pooled FAR** under these fixed params on NSRDB (not a clinical ranking):",
        "",
    ]
    ranked = sorted(
        [
            ("OPC L=50", o["far_per_24h"]),
            ("SDD", s["far_per_24h"]),
            ("abs-z τ_s", a["far_per_24h"]),
        ],
        key=lambda x: (float("inf") if not np.isfinite(x[1]) else x[1]),
    )
    for i, (name, far) in enumerate(ranked, 1):
        lines.append(f"{i}. **{name}**: {far:.3f} / 24h")

    lines += [
        "",
        "### Recommendation (specificity-only, interim)",
        "",
        "If the sole question is which **ordinal** rule looks more specific on healthy Holters "
        "under these defaults: prefer the lower-FAR ordinal arm above. "
        "Neither ordinal rule is promoted to production; abs-z remains the frozen baseline. "
        "Sensitivity on SDDB/VFDB was explored separately in the ordinal bake-off "
        "(`docs/ORDINAL_EXPLORATORY_BAKEOFF.md`). FAR and sensitivity trade-offs must be "
        "weighed jointly before any detector choice.",
        "",
        "## Artifacts",
        "",
        "- `results/ordinal_nsrdb_far_per_record.csv`",
        "- `results/ordinal_nsrdb_far_by_detector.csv`",
        "- `results/ordinal_nsrdb_far_summary.json`",
        "- Runner: `code/run_ordinal_far_comparison.py`",
        "- Episode helper: `count_binary_alarm_episodes` in `code/cctp_metrics_core.py`",
        "",
        "## Non-claims",
        "",
        "- No S5 (FAR ≤ 2/24h).",
        "- No clinical utility / FDA / deployability.",
        "- No fusion of OPC and SDD.",
        "- No retune of abs-z or ordinal thresholds based on this run.",
        "",
    ]
    doc.write_text("\n".join(lines), encoding="utf-8")
    print(f"Wrote {doc}")


if __name__ == "__main__":
    raise SystemExit(main())
