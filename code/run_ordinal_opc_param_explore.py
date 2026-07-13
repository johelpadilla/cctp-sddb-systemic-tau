#!/usr/bin/env python3
"""
Modest OPC (L, θ_D, θ_R) parameter exploration vs baseline L=50, θ_D=0.35, θ_R=5.

Exploratory only — not clinical optimization, not abs-z retune, not fusion.

Grid (fixed, modest, justified):
  L      ∈ {40, 50, 60, 70}   around L=50 companion scale for K=36
  θ_D    ∈ {0.30, 0.35, 0.40} around 0.35 diversity threshold
  θ_R    ∈ {4, 5, 6}          around 5 min consecutive low-div windows

Methodology (same as bake-off + NSRDB FAR):
  - Event sensitivity: first OPC alarm in post-basal pre-event window (SDDB n=11, VFDB n=22)
  - Control FAR: binary-alarm episodes, refractory 0.5 h, Phase-2 basal/search, cap 12 h (NSRDB n=18)
  - FAR = total_episodes / total_search_hours * 24

Symbol streams are built once per record and reused across all 36 cells.

Outputs under results/:
  ordinal_opc_param_explore_grid.csv
  ordinal_opc_param_explore_summary.json
  ordinal_opc_param_explore_report.md  (optional via --write-report)

Usage:
  python code/run_ordinal_opc_param_explore.py
  python code/run_ordinal_opc_param_explore.py --smoke  # tiny subset for tests
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cctp_metrics_core import (  # noqa: E402
    FROZEN_MIN_CONSECUTIVE,
    FROZEN_Z_THRESHOLD,
    build_bivariate_proxy,
    count_binary_alarm_episodes,
    false_alarm_rate,
    get_event_and_windows,
    short_db_windows,
)
from ordinal_detectors.opc_detector import opc_detect, opc_first_alarm_index  # noqa: E402
from recd_ordinal_levels import generate_multivariate_symbols  # noqa: E402

BASE = Path(__file__).resolve().parent.parent
DATA = BASE / "data"
RR_EXT = DATA / "rr_external"
RES = BASE / "results"
RES.mkdir(parents=True, exist_ok=True)

# Encoding — frozen with bake-off / FAR runners
M_EMB = 3
DELAY = 1
BP_ALPHABET = 6
K_JOINT = BP_ALPHABET * BP_ALPHABET  # 36

# Modest justified grid (plan / objective)
DEFAULT_L_VALUES: Tuple[int, ...] = (40, 50, 60, 70)
DEFAULT_THETA_D_VALUES: Tuple[float, ...] = (0.30, 0.35, 0.40)
DEFAULT_THETA_R_VALUES: Tuple[int, ...] = (4, 5, 6)

BASELINE_L = 50
BASELINE_THETA_D = 0.35
BASELINE_THETA_R = 5

ANALYTIC_SDDB = ["30", "31", "32", "35", "36", "38", "45", "47", "50", "51"]
EXTRA_SDDB = ["44"]  # full n=11 with 44

BASAL_HOURS = 2.0
CONTROL_MAX_HOURS = 12.0
REFRACTORY_H = 0.5

# Exploratory “not much worse FAR” slack for candidate selection
FAR_SLACK_FACTOR = 1.5  # allow up to 1.5× baseline FAR
FAR_SLACK_ABS = 2.0  # or +2 alarms/24h, whichever is larger slack

GRID_JUSTIFICATION = {
    "L": (
        "Around L=50 companion used for K=36 so min(L,K)/K can exceed θ_D "
        "(need L > θ_D·K ≈ 12.6). Explore ±10–20: {40,50,60,70}."
    ),
    "theta_D": (
        "Around 0.35 baseline diversity collapse threshold. Lower (0.30) is "
        "stricter collapse → fewer alarms; higher (0.40) is looser → more sensitive."
    ),
    "theta_R": (
        "Around 5 consecutive low-diversity windows. Lower (4) eases alarm "
        "(more sens / more FAR); higher (6) requires longer collapse."
    ),
    "scope": (
        "Modest 4×3×3=36 cell product — not exhaustive grid search or nested CV. "
        "Exploratory only; abs-z remains frozen; no fusion with SDD."
    ),
}


# ---------------------------------------------------------------------------
# Pure helpers (unit-tested)
# ---------------------------------------------------------------------------


def build_param_grid(
    L_values: Sequence[int] = DEFAULT_L_VALUES,
    theta_D_values: Sequence[float] = DEFAULT_THETA_D_VALUES,
    theta_R_values: Sequence[int] = DEFAULT_THETA_R_VALUES,
) -> List[Tuple[int, float, int]]:
    """Cartesian product of (L, θ_D, θ_R), stable nested order L → θ_D → θ_R."""
    grid: List[Tuple[int, float, int]] = []
    for L in L_values:
        for theta_D in theta_D_values:
            for theta_R in theta_R_values:
                grid.append((int(L), float(theta_D), int(theta_R)))
    return grid


def is_baseline(
    L: int,
    theta_D: float,
    theta_R: int,
    *,
    base_L: int = BASELINE_L,
    base_theta_D: float = BASELINE_THETA_D,
    base_theta_R: int = BASELINE_THETA_R,
    atol: float = 1e-12,
) -> bool:
    return (
        int(L) == int(base_L)
        and abs(float(theta_D) - float(base_theta_D)) <= atol
        and int(theta_R) == int(base_theta_R)
    )


def param_key(L: int, theta_D: float, theta_R: int) -> str:
    return f"L{int(L)}_D{float(theta_D):.2f}_R{int(theta_R)}"


def aggregate_sensitivity_from_flags(
    flags: Sequence[Dict[str, Any]],
) -> Dict[str, Any]:
    """
    Aggregate event sensitivity from per-record dicts with keys:
      source ('sddb'|'vfdb'), alarmed (0/1).
    """
    def _block(source: Optional[str] = None) -> Dict[str, Any]:
        sub = list(flags)
        if source is not None:
            sub = [r for r in sub if r.get("source") == source]
        n = len(sub)
        det = sum(1 for r in sub if int(r.get("alarmed", 0)) == 1)
        return {
            "n": n,
            "n_detected": det,
            "sensitivity": float(det / n) if n else float("nan"),
        }

    return {
        "sddb": _block("sddb"),
        "vfdb": _block("vfdb"),
        "all_events": _block(None),
    }


def aggregate_far_from_episodes(
    episode_rows: Sequence[Dict[str, float]],
) -> Dict[str, Any]:
    """
    Pool FAR from rows with n_episodes, search_hours, alarmed (shipped false_alarm_rate).
    """
    fr = false_alarm_rate(list(episode_rows))
    return {
        "n_controls": int(fr.get("n_controls", 0)),
        "total_episodes": float(fr.get("total_episodes", 0.0)),
        "total_search_hours": float(fr.get("total_search_hours", 0.0)),
        "far_per_24h": float(fr.get("far_per_24h", float("nan"))),
        "fraction_alarmed": float(fr.get("fraction_alarmed", float("nan")))
        if fr.get("fraction_alarmed") is not None
        else float("nan"),
        "reason": str(fr.get("reason", "")),
    }


def baseline_deltas(
    sens_all: float,
    far_per_24h: float,
    base_sens_all: float,
    base_far: float,
) -> Dict[str, float]:
    return {
        "delta_sens_all": float(sens_all - base_sens_all),
        "delta_far_per_24h": float(far_per_24h - base_far),
        "sens_all_ratio": float(sens_all / base_sens_all)
        if base_sens_all > 0
        else float("nan"),
        "far_ratio": float(far_per_24h / base_far) if base_far > 0 else float("nan"),
    }


def far_within_slack(
    far: float,
    base_far: float,
    *,
    factor: float = FAR_SLACK_FACTOR,
    abs_slack: float = FAR_SLACK_ABS,
) -> bool:
    """True if FAR is not much worse than baseline (exploratory slack)."""
    if not np.isfinite(far) or not np.isfinite(base_far):
        return False
    limit = max(base_far * factor, base_far + abs_slack)
    return far <= limit + 1e-12


def classify_cell(
    sens_all: float,
    far: float,
    base_sens: float,
    base_far: float,
    *,
    factor: float = FAR_SLACK_FACTOR,
    abs_slack: float = FAR_SLACK_ABS,
) -> str:
    """
    Exploratory labels (not clinical gates):
      baseline | higher_sens_far_ok | higher_sens_far_worse |
      lower_sens_lower_far | lower_sens_higher_far | similar | other
    """
    if abs(sens_all - base_sens) < 1e-12 and abs(far - base_far) < 1e-9:
        return "baseline"
    sens_up = sens_all > base_sens + 1e-12
    sens_down = sens_all < base_sens - 1e-12
    far_ok = far_within_slack(far, base_far, factor=factor, abs_slack=abs_slack)
    far_better = far < base_far - 1e-9
    far_worse = far > base_far + 1e-9
    if sens_up and far_ok:
        return "higher_sens_far_ok"
    if sens_up and not far_ok:
        return "higher_sens_far_worse"
    if sens_down and far_better:
        return "lower_sens_lower_far"
    if sens_down and far_worse:
        return "lower_sens_higher_far"
    if abs(sens_all - base_sens) < 1e-12 and far_ok:
        return "same_sens_far_ok"
    return "other"


def build_recommendation(
    cell_rows: Sequence[Dict[str, Any]],
    *,
    factor: float = FAR_SLACK_FACTOR,
    abs_slack: float = FAR_SLACK_ABS,
) -> Dict[str, Any]:
    """
    Honest adopt / keep-baseline decision from labeled cells.

    Prefer cells with higher_sens_far_ok ranked by (delta_sens desc, delta_far asc).
    If none, recommend keep_baseline.
    """
    baseline = next((r for r in cell_rows if r.get("is_baseline")), None)
    if baseline is None:
        return {
            "decision": "keep_baseline",
            "reason": "baseline row missing from results",
            "candidates": [],
            "adopt_params": None,
            "exploratory_only": True,
            "clinical_claim": False,
        }

    candidates = [r for r in cell_rows if r.get("class") == "higher_sens_far_ok"]
    candidates = sorted(
        candidates,
        key=lambda r: (
            -float(r.get("delta_sens_all", 0.0)),
            float(r.get("delta_far_per_24h", 0.0)),
            int(r.get("L", 0)),
            float(r.get("theta_D", 0.0)),
            int(r.get("theta_R", 0)),
        ),
    )
    top = candidates[:5]
    if not candidates:
        # Check if any higher sens exists at all
        higher = [r for r in cell_rows if float(r.get("delta_sens_all", 0)) > 1e-12]
        if higher:
            reason = (
                "No cell increased all-event sensitivity while keeping FAR within "
                f"exploratory slack (≤ max({factor}× baseline, baseline+{abs_slack})). "
                "Some cells raised sensitivity only with substantially higher FAR."
            )
        else:
            reason = (
                "No cell in the modest grid increased all-event sensitivity above "
                "baseline OPC L=50, θ_D=0.35, θ_R=5."
            )
        return {
            "decision": "keep_baseline",
            "reason": reason,
            "candidates": [],
            "adopt_params": {
                "L": int(baseline["L"]),
                "theta_D": float(baseline["theta_D"]),
                "theta_R": int(baseline["theta_R"]),
            },
            "baseline_sens_all": float(baseline["sens_all"]),
            "baseline_far_per_24h": float(baseline["far_per_24h"]),
            "far_slack_factor": factor,
            "far_slack_abs": abs_slack,
            "exploratory_only": True,
            "clinical_claim": False,
            "superiority_claim": False,
        }

    best = candidates[0]
    return {
        "decision": "consider_adopt_exploratory",
        "reason": (
            "At least one grid cell raised sensitivity with FAR within exploratory "
            "slack vs baseline. Still exploratory — not a clinical retune; small n."
        ),
        "candidates": [
            {
                "param_key": c["param_key"],
                "L": int(c["L"]),
                "theta_D": float(c["theta_D"]),
                "theta_R": int(c["theta_R"]),
                "sens_all": float(c["sens_all"]),
                "far_per_24h": float(c["far_per_24h"]),
                "delta_sens_all": float(c["delta_sens_all"]),
                "delta_far_per_24h": float(c["delta_far_per_24h"]),
            }
            for c in top
        ],
        "adopt_params": {
            "L": int(best["L"]),
            "theta_D": float(best["theta_D"]),
            "theta_R": int(best["theta_R"]),
        },
        "baseline_sens_all": float(baseline["sens_all"]),
        "baseline_far_per_24h": float(baseline["far_per_24h"]),
        "far_slack_factor": factor,
        "far_slack_abs": abs_slack,
        "exploratory_only": True,
        "clinical_claim": False,
        "superiority_claim": False,
    }


def qualitative_parameter_effects(
    cell_rows: Sequence[Dict[str, Any]],
) -> List[str]:
    """
    Qualitative notes from grid margins (mean sens/FAR holding other axes when possible).
    Uses simple one-way means — descriptive only.
    """
    notes: List[str] = []
    if not cell_rows:
        return ["No cells to summarize."]

    def _mean(rows: List[Dict[str, Any]], key: str) -> float:
        vals = [float(r[key]) for r in rows if np.isfinite(float(r.get(key, float("nan"))))]
        return float(np.mean(vals)) if vals else float("nan")

    for axis, key, fmt in (
        ("L", "L", "d"),
        ("θ_D", "theta_D", ".2f"),
        ("θ_R", "theta_R", "d"),
    ):
        levels = sorted({r[key] for r in cell_rows})
        parts = []
        for lv in levels:
            sub = [r for r in cell_rows if r[key] == lv]
            parts.append(
                f"{axis}={format(lv, fmt)}: mean sens_all={_mean(sub, 'sens_all'):.3f}, "
                f"mean FAR={_mean(sub, 'far_per_24h'):.2f}"
            )
        notes.append(" | ".join(parts))

    # Directional heuristics from sorted levels
    L_levels = sorted({int(r["L"]) for r in cell_rows})
    if len(L_levels) >= 2:
        low_L = [r for r in cell_rows if int(r["L"]) == L_levels[0]]
        high_L = [r for r in cell_rows if int(r["L"]) == L_levels[-1]]
        ds = _mean(high_L, "sens_all") - _mean(low_L, "sens_all")
        df = _mean(high_L, "far_per_24h") - _mean(low_L, "far_per_24h")
        notes.append(
            f"L effect (high vs low): Δmean_sens={ds:+.3f}, Δmean_FAR={df:+.2f} "
            "(larger L averages diversity over more symbols — tends to dampen brief collapses)."
        )

    D_levels = sorted({float(r["theta_D"]) for r in cell_rows})
    if len(D_levels) >= 2:
        low_D = [r for r in cell_rows if float(r["theta_D"]) == D_levels[0]]
        high_D = [r for r in cell_rows if float(r["theta_D"]) == D_levels[-1]]
        ds = _mean(high_D, "sens_all") - _mean(low_D, "sens_all")
        df = _mean(high_D, "far_per_24h") - _mean(low_D, "far_per_24h")
        notes.append(
            f"θ_D effect (high vs low): Δmean_sens={ds:+.3f}, Δmean_FAR={df:+.2f} "
            "(higher θ_D is a looser diversity threshold → easier low-div declaration)."
        )

    R_levels = sorted({int(r["theta_R"]) for r in cell_rows})
    if len(R_levels) >= 2:
        low_R = [r for r in cell_rows if int(r["theta_R"]) == R_levels[0]]
        high_R = [r for r in cell_rows if int(r["theta_R"]) == R_levels[-1]]
        ds = _mean(high_R, "sens_all") - _mean(low_R, "sens_all")
        df = _mean(high_R, "far_per_24h") - _mean(low_R, "far_per_24h")
        notes.append(
            f"θ_R effect (high vs low): Δmean_sens={ds:+.3f}, Δmean_FAR={df:+.2f} "
            "(higher θ_R requires longer consecutive collapse before alarm)."
        )

    notes.append(
        "Framing: descriptive grid margins only — not causal clinical effects; n is small."
    )
    return notes


# ---------------------------------------------------------------------------
# Data I/O / symbolization (shared with bake-off / FAR)
# ---------------------------------------------------------------------------


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


def joint_bivariate_symbols(rr: np.ndarray) -> Tuple[np.ndarray, int, int]:
    X = build_bivariate_proxy(np.asarray(rr, dtype=float))
    S = generate_multivariate_symbols(X, m=M_EMB, delay=DELAY)
    if S.size == 0 or S.shape[1] < 2:
        return np.array([], dtype=np.int64), K_JOINT, (M_EMB - 1) * DELAY
    sigma = (S[:, 0].astype(np.int64) * BP_ALPHABET) + S[:, 1].astype(np.int64)
    offset = (M_EMB - 1) * DELAY
    return sigma, K_JOINT, offset


def hours_to_symbol_index(t_hr_sym: np.ndarray, hr: float) -> int:
    if len(t_hr_sym) == 0:
        return 0
    idx = int(np.searchsorted(t_hr_sym, hr, side="left"))
    return max(0, min(idx, len(t_hr_sym)))


def hours_to_symbol_index_right(t_hr_sym: np.ndarray, hr: float) -> int:
    if len(t_hr_sym) == 0:
        return 0
    idx = int(np.searchsorted(t_hr_sym, hr, side="left"))
    return max(0, min(idx, len(t_hr_sym)))


def control_basal(total_h: float, basal_hours: float = BASAL_HOURS) -> Tuple[float, float]:
    basal = (0.25, min(basal_hours, total_h * 0.25))
    if basal[1] <= basal[0]:
        basal = (0.0, max(total_h * 0.2, 0.1))
    return basal


def resolve_windows(
    source: str,
    record: str,
    t_hr: np.ndarray,
    vfon_hr: float,
    total_h: float,
) -> Tuple[float, Tuple[float, float], Tuple[float, float], str]:
    if source == "sddb":
        event_hr, basal, approach = get_event_and_windows(record, t_hr, vfon_hr)
        return float(event_hr), basal, approach, "holter_analytic"
    event_hr, basal, approach, stratum = short_db_windows(vfon_hr, total_h)
    return float(event_hr), basal, approach, stratum


# ---------------------------------------------------------------------------
# Cached record packs
# ---------------------------------------------------------------------------


def load_event_pack(source: str, record: str, path: Path) -> Optional[Dict[str, Any]]:
    d = load_npz(path)
    rr = np.asarray(d["rr_ms"], dtype=float)
    t_sec = np.asarray(d["t_sec"], dtype=float)
    t_hr = t_sec / 3600.0
    vfon_hr = float(d["vfon_sec"]) / 3600.0
    total_h = float(d.get("total_hours", float(np.nanmax(t_hr))))
    event_hr, basal, approach, stratum = resolve_windows(
        source, record, t_hr, vfon_hr, total_h
    )
    sigma, K, offset = joint_bivariate_symbols(rr)
    if len(sigma) == 0:
        return None
    t_hr_sym = t_hr[offset : offset + len(sigma)]
    if len(t_hr_sym) < len(sigma):
        sigma = sigma[: len(t_hr_sym)]
    elif len(t_hr_sym) > len(sigma):
        t_hr_sym = t_hr_sym[: len(sigma)]
    b0, b1 = basal
    search_start = hours_to_symbol_index(t_hr_sym, b1)
    search_end = hours_to_symbol_index_right(t_hr_sym, event_hr)
    return {
        "source": source,
        "record": record,
        "path": path.name,
        "sigma": sigma,
        "K": K,
        "t_hr_sym": t_hr_sym,
        "event_hr": event_hr,
        "basal": basal,
        "approach": approach,
        "stratum": stratum,
        "search_start": search_start,
        "search_end": search_end,
    }


def load_control_pack(
    path: Path,
    *,
    basal_hours: float = BASAL_HOURS,
    max_hours: float = CONTROL_MAX_HOURS,
) -> Optional[Dict[str, Any]]:
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
    if len(rr) < 150 or total_h < basal_hours + 0.5:
        return None
    basal = control_basal(total_h, basal_hours)
    b0, b1 = basal
    search_start = b1
    search_end = total_h
    sigma, K, offset = joint_bivariate_symbols(rr)
    if len(sigma) < 60:
        return None
    t_hr_sym = t_hr[offset : offset + len(sigma)]
    if len(t_hr_sym) < len(sigma):
        sigma = sigma[: len(t_hr_sym)]
    n_sym = min(len(sigma), len(t_hr_sym))
    sigma = sigma[:n_sym]
    t_hr_sym = t_hr_sym[:n_sym]
    return {
        "record_id": rec,
        "path": path.name,
        "sigma": sigma,
        "K": K,
        "t_hr_sym": t_hr_sym,
        "basal": basal,
        "search_start": search_start,
        "search_end": search_end,
        "total_hours_used": total_h,
    }


def evaluate_opc_on_event(
    pack: Dict[str, Any],
    *,
    L: int,
    theta_D: float,
    theta_R: int,
) -> Dict[str, Any]:
    idx, _ = opc_first_alarm_index(
        pack["sigma"],
        L=L,
        theta_D=theta_D,
        theta_R=theta_R,
        K=pack["K"],
        search_start=pack["search_start"],
        search_end=pack["search_end"],
    )
    alarmed = 1 if idx is not None else 0
    lead = float("nan")
    det_hr = float("nan")
    if idx is not None:
        det_hr = float(pack["t_hr_sym"][idx])
        lead = float(pack["event_hr"] - det_hr)
    return {
        "source": pack["source"],
        "record": pack["record"],
        "alarmed": alarmed,
        "lead_time_h": lead,
        "detection_hr": det_hr,
        "first_alarm_idx": idx if idx is not None else "",
    }


def evaluate_opc_on_control(
    pack: Dict[str, Any],
    *,
    L: int,
    theta_D: float,
    theta_R: int,
    refractory_h: float = REFRACTORY_H,
) -> Dict[str, float]:
    out = opc_detect(
        pack["sigma"],
        L=L,
        theta_D=theta_D,
        theta_R=theta_R,
        K=pack["K"],
    )
    ep = count_binary_alarm_episodes(
        out["alarm"],
        pack["t_hr_sym"],
        search_start_hr=pack["search_start"],
        search_end_hr=pack["search_end"],
        refractory_h=refractory_h,
    )
    return {
        "n_episodes": float(ep["n_episodes"]),
        "search_hours": float(ep["search_hours"]),
        "alarmed": float(ep["alarmed"]),
        "first_alarm_hr": float(ep["first_alarm_hr"]),
    }


def evaluate_grid_cell(
    event_packs: Sequence[Dict[str, Any]],
    control_packs: Sequence[Dict[str, Any]],
    *,
    L: int,
    theta_D: float,
    theta_R: int,
    refractory_h: float = REFRACTORY_H,
) -> Dict[str, Any]:
    """Evaluate one (L, θ_D, θ_R) on cached packs; return aggregated row fields."""
    event_flags = [
        evaluate_opc_on_event(p, L=L, theta_D=theta_D, theta_R=theta_R)
        for p in event_packs
    ]
    sens = aggregate_sensitivity_from_flags(event_flags)
    ep_rows = [
        evaluate_opc_on_control(
            p, L=L, theta_D=theta_D, theta_R=theta_R, refractory_h=refractory_h
        )
        for p in control_packs
    ]
    far = aggregate_far_from_episodes(ep_rows)
    return {
        "L": int(L),
        "theta_D": float(theta_D),
        "theta_R": int(theta_R),
        "param_key": param_key(L, theta_D, theta_R),
        "is_baseline": is_baseline(L, theta_D, theta_R),
        "n_sddb": int(sens["sddb"]["n"]),
        "n_detected_sddb": int(sens["sddb"]["n_detected"]),
        "sens_sddb": float(sens["sddb"]["sensitivity"]),
        "n_vfdb": int(sens["vfdb"]["n"]),
        "n_detected_vfdb": int(sens["vfdb"]["n_detected"]),
        "sens_vfdb": float(sens["vfdb"]["sensitivity"]),
        "n_all": int(sens["all_events"]["n"]),
        "n_detected_all": int(sens["all_events"]["n_detected"]),
        "sens_all": float(sens["all_events"]["sensitivity"]),
        "n_controls": int(far["n_controls"]),
        "total_episodes": float(far["total_episodes"]),
        "total_search_hours": float(far["total_search_hours"]),
        "far_per_24h": float(far["far_per_24h"]),
        "fraction_alarmed": float(far["fraction_alarmed"]),
        "event_flags": event_flags,
        "control_episode_rows": ep_rows,
    }


def finalize_cell_rows(raw_rows: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Attach baseline deltas, class labels; drop bulky nested fields for table."""
    baseline = next((r for r in raw_rows if r.get("is_baseline")), None)
    if baseline is None:
        # synthetic baseline metrics from known pilot if missing (should not happen)
        base_sens = 14 / 33
        base_far = 3.7333769938233505
    else:
        base_sens = float(baseline["sens_all"])
        base_far = float(baseline["far_per_24h"])

    out: List[Dict[str, Any]] = []
    for r in raw_rows:
        deltas = baseline_deltas(
            float(r["sens_all"]),
            float(r["far_per_24h"]),
            base_sens,
            base_far,
        )
        cls = classify_cell(
            float(r["sens_all"]),
            float(r["far_per_24h"]),
            base_sens,
            base_far,
        )
        if r.get("is_baseline"):
            cls = "baseline"
        row = {
            "param_key": r["param_key"],
            "L": r["L"],
            "theta_D": r["theta_D"],
            "theta_R": r["theta_R"],
            "is_baseline": bool(r["is_baseline"]),
            "n_sddb": r["n_sddb"],
            "n_detected_sddb": r["n_detected_sddb"],
            "sens_sddb": r["sens_sddb"],
            "n_vfdb": r["n_vfdb"],
            "n_detected_vfdb": r["n_detected_vfdb"],
            "sens_vfdb": r["sens_vfdb"],
            "n_all": r["n_all"],
            "n_detected_all": r["n_detected_all"],
            "sens_all": r["sens_all"],
            "n_controls": r["n_controls"],
            "total_episodes": r["total_episodes"],
            "total_search_hours": r["total_search_hours"],
            "far_per_24h": r["far_per_24h"],
            "fraction_alarmed": r["fraction_alarmed"],
            "delta_sens_all": deltas["delta_sens_all"],
            "delta_far_per_24h": deltas["delta_far_per_24h"],
            "sens_all_ratio": deltas["sens_all_ratio"],
            "far_ratio": deltas["far_ratio"],
            "class": cls,
            "far_within_slack": far_within_slack(
                float(r["far_per_24h"]), base_far
            ),
        }
        out.append(row)
    return out


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


def write_report_md(
    path: Path,
    cell_rows: List[Dict[str, Any]],
    recommendation: Dict[str, Any],
    notes: List[str],
) -> None:
    baseline = next((r for r in cell_rows if r.get("is_baseline")), None)
    lines = [
        "# OPC parameter exploration (exploratory)",
        "",
        f"Generated: {datetime.now(timezone.utc).isoformat()}",
        "",
        "## Framing",
        "",
        "- Exploratory only — not clinical optimization, not S5, not FDA, not superiority.",
        "- Abs-z remains frozen (z=2, min_run=3); not retuned here.",
        "- No fusion with SDD.",
        "",
        "## Parameter grid (justified modest)",
        "",
        f"- **L** ∈ {list(DEFAULT_L_VALUES)} — {GRID_JUSTIFICATION['L']}",
        f"- **θ_D** ∈ {list(DEFAULT_THETA_D_VALUES)} — {GRID_JUSTIFICATION['theta_D']}",
        f"- **θ_R** ∈ {list(DEFAULT_THETA_R_VALUES)} — {GRID_JUSTIFICATION['theta_R']}",
        f"- {GRID_JUSTIFICATION['scope']}",
        "",
        f"- Cells evaluated: **{len(cell_rows)}**",
        "",
        "## Methodology",
        "",
        "- Sensitivity: first OPC alarm in post-basal pre-event window (SDDB n=11, VFDB n=22).",
        "- FAR: binary alarm episodes, refractory 0.5 h, Phase-2 basal/search, control cap 12 h, NSRDB n=18.",
        "- FAR formula: total_episodes / total_search_hours × 24.",
        "- Encoding: joint bivariate Bandt–Pompe m=3, K=36 (unchanged).",
        "",
        "## Baseline (L=50, θ_D=0.35, θ_R=5)",
        "",
    ]
    if baseline:
        lines.extend(
            [
                f"- sens_sddb = {baseline['sens_sddb']:.4f} "
                f"({baseline['n_detected_sddb']}/{baseline['n_sddb']})",
                f"- sens_vfdb = {baseline['sens_vfdb']:.4f} "
                f"({baseline['n_detected_vfdb']}/{baseline['n_vfdb']})",
                f"- sens_all = {baseline['sens_all']:.4f} "
                f"({baseline['n_detected_all']}/{baseline['n_all']})",
                f"- FAR = {baseline['far_per_24h']:.4f} /24h "
                f"({int(baseline['total_episodes'])} ep / "
                f"{baseline['total_search_hours']:.2f} h search)",
                "",
            ]
        )
    lines.extend(
        [
            "## Recommendation",
            "",
            f"- **Decision:** `{recommendation.get('decision')}`",
            f"- **Reason:** {recommendation.get('reason')}",
            f"- **Adopt params (if any):** `{recommendation.get('adopt_params')}`",
            "",
            "## Best candidates (higher sens, FAR within slack)",
            "",
        ]
    )
    cands = recommendation.get("candidates") or []
    if not cands:
        lines.append("_None._")
    else:
        lines.append(
            "| param_key | L | θ_D | θ_R | sens_all | FAR | Δsens | ΔFAR |"
        )
        lines.append("|---|---:|---:|---:|---:|---:|---:|---:|")
        for c in cands:
            lines.append(
                f"| {c['param_key']} | {c['L']} | {c['theta_D']:.2f} | {c['theta_R']} | "
                f"{c['sens_all']:.4f} | {c['far_per_24h']:.3f} | "
                f"{c['delta_sens_all']:+.4f} | {c['delta_far_per_24h']:+.3f} |"
            )
    lines.extend(["", "## Qualitative parameter effects", ""])
    for n in notes:
        lines.append(f"- {n}")
    lines.extend(
        [
            "",
            "## Full grid (sorted by sens_all desc, then FAR asc)",
            "",
            "| is_baseline | L | θ_D | θ_R | sens_sddb | sens_vfdb | sens_all | FAR | class |",
            "|---:|---:|---:|---:|---:|---:|---:|---:|---|",
        ]
    )
    sorted_rows = sorted(
        cell_rows,
        key=lambda r: (-float(r["sens_all"]), float(r["far_per_24h"])),
    )
    for r in sorted_rows:
        lines.append(
            f"| {int(r['is_baseline'])} | {r['L']} | {r['theta_D']:.2f} | {r['theta_R']} | "
            f"{r['sens_sddb']:.3f} | {r['sens_vfdb']:.3f} | {r['sens_all']:.3f} | "
            f"{r['far_per_24h']:.3f} | {r['class']} |"
        )
    lines.append("")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def collect_event_packs(*, include_sddb44: bool = True) -> List[Dict[str, Any]]:
    packs: List[Dict[str, Any]] = []
    sddb_ids = list(ANALYTIC_SDDB)
    if include_sddb44:
        sddb_ids += EXTRA_SDDB
    for rec in sddb_ids:
        path = DATA / f"rr_{rec}_clean.npz"
        if not path.exists():
            print(f"  WARN missing SDDB {rec}", flush=True)
            continue
        pack = load_event_pack("sddb", rec, path)
        if pack is None:
            print(f"  WARN empty symbols SDDB {rec}", flush=True)
            continue
        packs.append(pack)
    for path in sorted(RR_EXT.glob("vfdb_*_clean.npz")):
        rec = path.stem.replace("vfdb_", "").replace("_clean", "")
        pack = load_event_pack("vfdb", rec, path)
        if pack is None:
            print(f"  WARN empty symbols VFDB {rec}", flush=True)
            continue
        packs.append(pack)
    return packs


def collect_control_packs(
    *,
    basal_hours: float = BASAL_HOURS,
    max_hours: float = CONTROL_MAX_HOURS,
) -> List[Dict[str, Any]]:
    packs: List[Dict[str, Any]] = []
    for path in sorted(RR_EXT.glob("nsrdb_*_clean.npz")):
        pack = load_control_pack(path, basal_hours=basal_hours, max_hours=max_hours)
        if pack is None:
            print(f"  WARN skip control {path.name}", flush=True)
            continue
        packs.append(pack)
    return packs


def run_exploration(
    grid: Sequence[Tuple[int, float, int]],
    event_packs: Sequence[Dict[str, Any]],
    control_packs: Sequence[Dict[str, Any]],
    *,
    refractory_h: float = REFRACTORY_H,
    verbose: bool = True,
) -> Tuple[List[Dict[str, Any]], Dict[str, Any], List[str]]:
    raw_rows: List[Dict[str, Any]] = []
    n_cells = len(grid)
    for i, (L, theta_D, theta_R) in enumerate(grid, 1):
        if verbose:
            print(
                f"  [{i}/{n_cells}] OPC L={L} θ_D={theta_D} θ_R={theta_R} ...",
                flush=True,
            )
        cell = evaluate_grid_cell(
            event_packs,
            control_packs,
            L=L,
            theta_D=theta_D,
            theta_R=theta_R,
            refractory_h=refractory_h,
        )
        if verbose:
            print(
                f"      sens_all={cell['sens_all']:.3f} "
                f"({cell['n_detected_all']}/{cell['n_all']})  "
                f"FAR={cell['far_per_24h']:.3f}  baseline={cell['is_baseline']}",
                flush=True,
            )
        raw_rows.append(cell)

    cell_rows = finalize_cell_rows(raw_rows)
    recommendation = build_recommendation(cell_rows)
    notes = qualitative_parameter_effects(cell_rows)
    return cell_rows, recommendation, notes


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(
        description="Modest OPC L/θ_D/θ_R exploration vs baseline L=50"
    )
    ap.add_argument("--out-dir", default=str(RES), help="Output directory")
    ap.add_argument("--write-report", action="store_true", help="Write markdown report")
    ap.add_argument(
        "--smoke",
        action="store_true",
        help="Tiny grid + first 2 events + first 2 controls (for tests)",
    )
    ap.add_argument("--refractory-h", type=float, default=REFRACTORY_H)
    ap.add_argument("--control-max-hours", type=float, default=CONTROL_MAX_HOURS)
    ap.add_argument("--basal-hours", type=float, default=BASAL_HOURS)
    args = ap.parse_args(argv)

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print("=== OPC parameter exploration (exploratory) ===")
    print("Grid justification:")
    for k, v in GRID_JUSTIFICATION.items():
        print(f"  {k}: {v}")
    print(
        f"Frozen abs-z unchanged: z={FROZEN_Z_THRESHOLD}, "
        f"min_run={FROZEN_MIN_CONSECUTIVE}"
    )

    if args.smoke:
        grid = build_param_grid(
            L_values=(50, 60),
            theta_D_values=(0.35,),
            theta_R_values=(5,),
        )
        # ensure baseline present
        if (BASELINE_L, BASELINE_THETA_D, BASELINE_THETA_R) not in grid:
            grid = [(BASELINE_L, BASELINE_THETA_D, BASELINE_THETA_R)] + grid
        print("SMOKE mode: reduced grid and records")
    else:
        grid = build_param_grid()

    print(f"Grid cells: {len(grid)}")
    print("Loading event packs (symbolize once)...")
    event_packs = collect_event_packs(include_sddb44=True)
    if args.smoke:
        event_packs = event_packs[:2]
    print(f"  events: {len(event_packs)}")
    print("Loading control packs (symbolize once)...")
    control_packs = collect_control_packs(
        basal_hours=args.basal_hours, max_hours=args.control_max_hours
    )
    if args.smoke:
        control_packs = control_packs[:2]
    print(f"  controls: {len(control_packs)}")

    if not event_packs:
        print("ERROR: no event packs", file=sys.stderr)
        return 1
    if not control_packs:
        print("ERROR: no control packs", file=sys.stderr)
        return 1

    cell_rows, recommendation, notes = run_exploration(
        grid,
        event_packs,
        control_packs,
        refractory_h=args.refractory_h,
        verbose=True,
    )

    grid_path = out_dir / "ordinal_opc_param_explore_grid.csv"
    write_csv(grid_path, cell_rows)
    print(f"Wrote {grid_path}")

    baseline_row = next((r for r in cell_rows if r.get("is_baseline")), None)
    summary = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "phase": "ordinal_opc_param_explore",
        "exploratory_only": True,
        "clinical_claim": False,
        "superiority_claim": False,
        "s5_claim": False,
        "fusion": False,
        "absz_frozen": {
            "z_threshold": FROZEN_Z_THRESHOLD,
            "min_consecutive": FROZEN_MIN_CONSECUTIVE,
            "retuned": False,
        },
        "grid": {
            "L_values": list(DEFAULT_L_VALUES) if not args.smoke else sorted({g[0] for g in grid}),
            "theta_D_values": list(DEFAULT_THETA_D_VALUES)
            if not args.smoke
            else sorted({g[1] for g in grid}),
            "theta_R_values": list(DEFAULT_THETA_R_VALUES)
            if not args.smoke
            else sorted({g[2] for g in grid}),
            "n_cells": len(cell_rows),
            "justification": GRID_JUSTIFICATION,
            "baseline": {
                "L": BASELINE_L,
                "theta_D": BASELINE_THETA_D,
                "theta_R": BASELINE_THETA_R,
            },
        },
        "methodology": {
            "encoding": "joint_bivariate_bandt_pompe_m3",
            "K": K_JOINT,
            "event_sensitivity": (
                "first OPC alarm in post-basal pre-event symbol window "
                "(opc_first_alarm_index)"
            ),
            "far_formula": "FAR = total_episodes / total_search_hours * 24",
            "episode_definition": (
                "Binary OPC alarm; each True sample outside refractory starts episode; "
                f"refractory {args.refractory_h} h (count_binary_alarm_episodes)"
            ),
            "basal_window_controls": "Phase-2 style: (0.25, min(basal_hours, 0.25*total_h))",
            "search_window_controls": "After basal end through end of capped recording",
            "control_cap_h": args.control_max_hours,
            "refractory_h": args.refractory_h,
            "n_sddb_expected": 11 if not args.smoke else None,
            "n_vfdb_expected": 22 if not args.smoke else None,
            "n_nsrdb_expected": 18 if not args.smoke else None,
            "device_mismatch": True,
            "device_mismatch_note": (
                "NSRDB is rhythm-healthy Holter ECG — NOT device-matched to VFDB/CU. "
                "FAR is interim public upper-bound estimate only."
            ),
        },
        "n_events_processed": {
            "sddb": sum(1 for p in event_packs if p["source"] == "sddb"),
            "vfdb": sum(1 for p in event_packs if p["source"] == "vfdb"),
            "all": len(event_packs),
        },
        "n_controls_processed": len(control_packs),
        "baseline_row": baseline_row,
        "cells": cell_rows,
        "recommendation": recommendation,
        "qualitative_parameter_notes": notes,
        "smoke": bool(args.smoke),
        "expected_baseline_reference": {
            "sens_all": 14 / 33,
            "far_per_24h": 3.7333769938233505,
            "source": "results/ordinal_exploratory_summary.json opc_L50_companion + "
            "ordinal_nsrdb_far_summary.json pooled_far.opc_L50",
        },
    }

    summary_path = out_dir / "ordinal_opc_param_explore_summary.json"
    with summary_path.open("w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2, allow_nan=True)
        f.write("\n")
    print(f"Wrote {summary_path}")

    if args.write_report or not args.smoke:
        report_path = out_dir / "ordinal_opc_param_explore_report.md"
        # Always write report for full runs; for smoke only if requested
        if not args.smoke or args.write_report:
            write_report_md(report_path, cell_rows, recommendation, notes)
            print(f"Wrote {report_path}")

    print("\n=== Recommendation ===")
    print(f"  decision: {recommendation.get('decision')}")
    print(f"  reason:   {recommendation.get('reason')}")
    if baseline_row:
        print(
            f"  baseline: sens_all={baseline_row['sens_all']:.4f} "
            f"FAR={baseline_row['far_per_24h']:.4f}"
        )
    print("=== Qualitative notes ===")
    for n in notes:
        print(f"  - {n}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
