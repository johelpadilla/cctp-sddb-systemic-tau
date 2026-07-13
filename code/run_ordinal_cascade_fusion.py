#!/usr/bin/env python3
"""
Exploratory light cascade fusion evaluation: SDD → OPC L=50 confirm (±5 min).

Head-to-head on the same cohorts/methods as existing ordinal work:
  - Sensitivity on SDDB (11) and VFDB (22): detect before event; lead-time when detected
  - FAR on NSRDB (18): episode counting, refractory 0.5 h,
    FAR = total_episodes / total_search_hours × 24
  - Arms: cascade | OPC L=50 alone | SDD alone | abs-z τ_s frozen

Cascade rule (fixed, not retuned):
  SDD candidate (L_c=50, θ_TV=0.35, θ_S=1) becomes cascade alarm only if
  OPC L=50 (θ_D=0.35, θ_R=5) also alarms within ±5 minutes (closed interval)
  on the shared joint-bivariate symbol timebase.

Does NOT modify opc_detect / sdd_detect / abs-z production thresholds.
Does NOT claim clinical validity, S5, or superiority.

Outputs under results/:
  ordinal_cascade_per_record.csv
  ordinal_cascade_nsrdb_far_per_record.csv
  ordinal_cascade_comparison.csv
  ordinal_cascade_gain_loss.csv
  ordinal_cascade_summary.json

Write-up: docs/ORDINAL_CASCADE_FUSION.md (with --write-doc).
"""
from __future__ import annotations

import argparse
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
    detect_lead_time,
    false_alarm_rate,
    get_event_and_windows,
    short_db_windows,
)
from ordinal_detectors.cascade_fusion import (
    CONFIRM_WINDOW_H,
    CONFIRM_WINDOW_MIN,
    cascade_first_causal_detection,
    cascade_sdd_confirm_opc,
)
from ordinal_detectors.opc_detector import opc_detect, opc_first_alarm_index
from ordinal_detectors.sdd_detector import sdd_detect, sdd_first_alarm_index
from recd_ordinal_levels import generate_multivariate_symbols

BASE = Path(__file__).resolve().parent.parent
DATA = BASE / "data"
RR_EXT = DATA / "rr_external"
RES = BASE / "results"
RES.mkdir(parents=True, exist_ok=True)

ANALYTIC_SDDB = ["30", "31", "32", "35", "36", "38", "45", "47", "50", "51"]
EXTRA_SDDB = ["44"]

M_EMB = 3
DELAY = 1
BP_ALPHABET = 6
K_JOINT = BP_ALPHABET * BP_ALPHABET  # 36

# Fixed params — same as bake-off L=50 companion / FAR comparison (not retuned)
OPC_L = 50
OPC_THETA_D = 0.35
OPC_THETA_R = 5
SDD_L_C = 50
SDD_THETA_TV = 0.35
SDD_THETA_S = 1

BASAL_HOURS = 2.0
CONTROL_MAX_HOURS = 12.0
REFRACTORY_H = 0.5

# This arm is exploratory fusion; singleton FAR scripts remain fusion=False.
EXPERIMENT_FUSION = True
EXPERIMENT_LABEL = "cascade_sdd_confirm_opc_L50"


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


def lead_from_symbol_alarm(
    t_hr_sym: np.ndarray,
    alarm_idx: Optional[int],
    event_hr: float,
) -> Tuple[float, float, int]:
    if alarm_idx is None:
        return float("nan"), float("nan"), 0
    det_hr = float(t_hr_sym[alarm_idx])
    return float(event_hr - det_hr), det_hr, 1


def control_basal(total_h: float, basal_hours: float = BASAL_HOURS) -> Tuple[float, float]:
    basal = (0.25, min(basal_hours, total_h * 0.25))
    if basal[1] <= basal[0]:
        basal = (0.0, max(total_h * 0.2, 0.1))
    return basal


def write_csv(path: Path, rows: List[dict]) -> None:
    if not rows:
        path.write_text("")
        return
    keys: List[str] = []
    seen = set()
    for r in rows:
        for k in r.keys():
            if k not in seen:
                seen.add(k)
                keys.append(k)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in keys})


def process_event_record(
    source: str,
    record: str,
    path: Path,
    *,
    confirm_window_h: float = CONFIRM_WINDOW_H,
) -> dict:
    """Sensitivity row: OPC L=50, SDD, cascade, abs-z on one event record."""
    d = load_npz(path)
    rr = np.asarray(d["rr_ms"], dtype=float)
    t_sec = np.asarray(d["t_sec"], dtype=float)
    t_hr = t_sec / 3600.0
    vfon_hr = float(d["vfon_sec"]) / 3600.0
    total_h = float(d.get("total_hours", float(np.nanmax(t_hr))))
    event_label = str(d.get("event_label", "VF/VT" if source == "vfdb" else "SDDB"))

    event_hr, basal, approach, stratum = resolve_windows(
        source, record, t_hr, vfon_hr, total_h
    )
    b0, b1 = basal

    sigma, K, offset = joint_bivariate_symbols(rr)
    base_meta = {
        "source": source,
        "record": record,
        "event_label": event_label,
        "duration_stratum": stratum,
        "event_hr": event_hr,
        "basal_start": b0,
        "basal_end": b1,
        "approach_start": approach[0],
        "approach_end": approach[1],
        "confirm_window_h": confirm_window_h,
        "confirm_window_min": confirm_window_h * 60.0,
        "fusion": EXPERIMENT_FUSION,
        "experiment": EXPERIMENT_LABEL,
    }

    if len(sigma) == 0:
        return {
            **base_meta,
            "error": "empty_symbol_stream",
            "opc_alarmed": 0,
            "sdd_alarmed": 0,
            "cascade_alarmed": 0,
            "absz_alarmed": 0,
        }

    t_hr_sym = t_hr[offset : offset + len(sigma)]
    if len(t_hr_sym) < len(sigma):
        sigma = sigma[: len(t_hr_sym)]
    elif len(t_hr_sym) > len(sigma):
        t_hr_sym = t_hr_sym[: len(sigma)]

    basal_i0 = hours_to_symbol_index(t_hr_sym, b0)
    basal_i1 = hours_to_symbol_index_right(t_hr_sym, b1)
    if basal_i1 <= basal_i0:
        event_i = hours_to_symbol_index_right(t_hr_sym, event_hr)
        basal_i0 = 0
        basal_i1 = max(1, event_i // 3)

    search_start = hours_to_symbol_index(t_hr_sym, b1)
    search_end = hours_to_symbol_index_right(t_hr_sym, event_hr)

    # --- OPC L=50 alone ---
    opc_idx, opc_out = opc_first_alarm_index(
        sigma,
        L=OPC_L,
        theta_D=OPC_THETA_D,
        theta_R=OPC_THETA_R,
        K=K,
        search_start=search_start,
        search_end=search_end,
    )
    opc_lead, opc_det_hr, opc_alarmed = lead_from_symbol_alarm(t_hr_sym, opc_idx, event_hr)

    # --- SDD alone ---
    basal_range = (basal_i0, basal_i1)
    if basal_range[1] - basal_range[0] < 2:
        sdd_idx, sdd_out = None, {
            "alarm": np.zeros(len(sigma), dtype=np.int8),
            "tv": np.full(len(sigma), np.nan),
        }
        sdd_lead, sdd_det_hr, sdd_alarmed = float("nan"), float("nan"), 0
        sdd_error = "basal_too_short_for_sdd"
    else:
        sdd_idx, sdd_out = sdd_first_alarm_index(
            sigma,
            basal_range,
            L_c=SDD_L_C,
            theta_TV=SDD_THETA_TV,
            theta_S=SDD_THETA_S,
            K=K,
            search_start=search_start,
            search_end=search_end,
            mask_basal=True,
        )
        sdd_lead, sdd_det_hr, sdd_alarmed = lead_from_symbol_alarm(
            t_hr_sym, sdd_idx, event_hr
        )
        sdd_error = ""

    # Full streams for cascade (same detectors). Causal confirm: OPC only if t < event_hr;
    # decision_time = max(t_SDD, t_OPC) must be < event_hr (no post-event look-ahead).
    opc_full = opc_detect(
        sigma, L=OPC_L, theta_D=OPC_THETA_D, theta_R=OPC_THETA_R, K=K
    )
    if basal_range[1] - basal_range[0] >= 2:
        sdd_full = sdd_detect(
            sigma,
            basal_range,
            L_c=SDD_L_C,
            theta_TV=SDD_THETA_TV,
            theta_S=SDD_THETA_S,
            K=K,
            mask_basal=True,
        )
    else:
        sdd_full = sdd_out

    cas_det = cascade_first_causal_detection(
        sdd_full["alarm"],
        opc_full["alarm"],
        t_hr_sym,
        confirm_window_h=confirm_window_h,
        search_start_hr=float(b1),
        event_hr=float(event_hr),
    )
    cas_out = cas_det["cascade_out"]
    cas_idx = cas_det["sdd_idx"]
    cas_lead = cas_det["lead_time_h"]
    cas_det_hr = cas_det["detection_hr"]
    cas_alarmed = int(cas_det["alarmed"])

    # --- abs-z frozen ---
    tau = compute_tau_series(rr)
    absz = detect_lead_time(
        tau,
        t_hr,
        event_hr,
        basal,
        z_threshold=FROZEN_Z_THRESHOLD,
        min_consecutive=FROZEN_MIN_CONSECUTIVE,
        use_abs=True,
    )

    return {
        **base_meta,
        "n_symbols": int(len(sigma)),
        "K": int(K),
        "basal_idx_start": int(basal_i0),
        "basal_idx_end": int(basal_i1),
        "search_start_idx": int(search_start),
        "search_end_idx": int(search_end),
        "sdd_error": sdd_error,
        # OPC L=50
        "opc_L": OPC_L,
        "opc_theta_D": OPC_THETA_D,
        "opc_theta_R": OPC_THETA_R,
        "opc_alarmed": int(opc_alarmed),
        "opc_lead_time_h": opc_lead,
        "opc_detection_hr": opc_det_hr,
        "opc_first_alarm_idx": opc_idx if opc_idx is not None else "",
        # SDD
        "sdd_L_c": SDD_L_C,
        "sdd_theta_TV": SDD_THETA_TV,
        "sdd_theta_S": SDD_THETA_S,
        "sdd_alarmed": int(sdd_alarmed),
        "sdd_lead_time_h": sdd_lead,
        "sdd_detection_hr": sdd_det_hr,
        "sdd_first_alarm_idx": sdd_idx if sdd_idx is not None else "",
        # Cascade
        "cascade_alarmed": int(cas_alarmed),
        "cascade_lead_time_h": cas_lead,
        "cascade_detection_hr": cas_det_hr,
        "cascade_first_alarm_idx": cas_idx if cas_idx is not None else "",
        "cascade_sdd_hr": cas_det.get("sdd_hr", float("nan")),
        "cascade_opc_confirm_hr": cas_det.get("opc_confirm_hr", float("nan")),
        "cascade_causal": True,
        "cascade_opc_max_hr": float(event_hr),
        "cascade_sdd_candidates": int(cas_out.get("sdd_candidate_count", 0)),
        "cascade_opc_alarms": int(cas_out.get("opc_alarm_count", 0)),
        "cascade_confirmed_samples": int(cas_out.get("cascade_alarm_count", 0)),
        # abs-z
        "absz_z_threshold": FROZEN_Z_THRESHOLD,
        "absz_min_consecutive": FROZEN_MIN_CONSECUTIVE,
        "absz_alarmed": int(absz["alarmed"]),
        "absz_lead_time_h": absz["lead_time_h"],
        "absz_detection_hr": absz["detection_hr"],
    }


def process_control(
    path: Path,
    *,
    basal_hours: float = BASAL_HOURS,
    max_hours: float = CONTROL_MAX_HOURS,
    refractory_h: float = REFRACTORY_H,
    confirm_window_h: float = CONFIRM_WINDOW_H,
) -> Dict[str, Any]:
    """NSRDB control: episode FAR for cascade + three singleton arms."""
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
        "fusion": EXPERIMENT_FUSION,
    }
    if len(rr) < W_TAU + 50 or total_h < basal_hours + 0.5:
        return empty

    basal = control_basal(total_h, basal_hours)
    b0, b1 = basal
    search_start = b1
    search_end = total_h

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
        }

    t_hr_sym = t_hr[offset : offset + len(sigma)]
    if len(t_hr_sym) < len(sigma):
        sigma = sigma[: len(t_hr_sym)]
    n_sym = min(len(sigma), len(t_hr_sym))
    sigma = sigma[:n_sym]
    t_hr_sym = t_hr_sym[:n_sym]

    basal_i0 = hours_to_symbol_index(t_hr_sym, b0)
    basal_i1 = hours_to_symbol_index_right(t_hr_sym, b1)
    if basal_i1 - basal_i0 < max(2, SDD_L_C // 5):
        basal_i0 = 0
        basal_i1 = max(SDD_L_C, min(len(sigma) // 5, len(sigma) // 4))
        if basal_i1 <= basal_i0:
            basal_i1 = min(len(sigma), basal_i0 + SDD_L_C)

    opc_out = opc_detect(
        sigma, L=OPC_L, theta_D=OPC_THETA_D, theta_R=OPC_THETA_R, K=K
    )
    opc_ep = count_binary_alarm_episodes(
        opc_out["alarm"],
        t_hr_sym,
        search_start_hr=search_start,
        search_end_hr=search_end,
        refractory_h=refractory_h,
    )

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

    cas_out = cascade_sdd_confirm_opc(
        sdd_out["alarm"],
        opc_out["alarm"],
        t_hr_sym,
        confirm_window_h=confirm_window_h,
    )
    cas_ep = count_binary_alarm_episodes(
        cas_out["alarm"],
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
        "confirm_window_h": confirm_window_h,
        "confirm_window_min": confirm_window_h * 60.0,
        "n_symbols": int(len(sigma)),
        "K": int(K),
        "basal_idx_start": int(basal_i0),
        "basal_idx_end": int(basal_i1),
        "fusion": EXPERIMENT_FUSION,
        "experiment": EXPERIMENT_LABEL,
        "opc_n_episodes": opc_ep["n_episodes"],
        "opc_search_hours": opc_ep["search_hours"],
        "opc_first_alarm_hr": opc_ep["first_alarm_hr"],
        "opc_alarmed": opc_ep["alarmed"],
        "opc_far_per_24h": _far(opc_ep),
        "sdd_n_episodes": sdd_ep["n_episodes"],
        "sdd_search_hours": sdd_ep["search_hours"],
        "sdd_first_alarm_hr": sdd_ep["first_alarm_hr"],
        "sdd_alarmed": sdd_ep["alarmed"],
        "sdd_far_per_24h": _far(sdd_ep),
        "cascade_n_episodes": cas_ep["n_episodes"],
        "cascade_search_hours": cas_ep["search_hours"],
        "cascade_first_alarm_hr": cas_ep["first_alarm_hr"],
        "cascade_alarmed": cas_ep["alarmed"],
        "cascade_far_per_24h": _far(cas_ep),
        "cascade_sdd_candidates": cas_out["sdd_candidate_count"],
        "cascade_confirmed_samples": cas_out["cascade_alarm_count"],
        "absz_n_episodes": absz_ep["n_episodes"],
        "absz_search_hours": absz_ep["search_hours"],
        "absz_first_alarm_hr": absz_ep["first_alarm_hr"],
        "absz_alarmed": absz_ep["alarmed"],
        "absz_far_per_24h": _far(absz_ep),
        "absz_z_threshold": FROZEN_Z_THRESHOLD,
        "absz_min_consecutive": FROZEN_MIN_CONSECUTIVE,
    }


def _sens_block(rows: List[dict], arm: str, source: Optional[str] = None) -> dict:
    key_a = f"{arm}_alarmed"
    key_l = f"{arm}_lead_time_h"
    sub = [r for r in rows if source is None or r.get("source") == source]
    n = len(sub)
    n_det = sum(1 for r in sub if int(r.get(key_a, 0) or 0) == 1)
    leads = [
        float(r[key_l])
        for r in sub
        if int(r.get(key_a, 0) or 0) == 1
        and np.isfinite(float(r.get(key_l, float("nan"))))
    ]
    return {
        "n": n,
        "n_detected": n_det,
        "sensitivity": float(n_det / n) if n else float("nan"),
        "median_lead_h": float(np.median(leads)) if leads else float("nan"),
        "mean_lead_h": float(np.mean(leads)) if leads else float("nan"),
    }


def gain_loss_rows(event_rows: List[dict]) -> List[dict]:
    """Per-record and cohort gain/loss of cascade vs each singleton."""
    out: List[dict] = []
    for r in event_rows:
        cas = int(r.get("cascade_alarmed", 0) or 0)
        opc = int(r.get("opc_alarmed", 0) or 0)
        sdd = int(r.get("sdd_alarmed", 0) or 0)
        absz = int(r.get("absz_alarmed", 0) or 0)
        out.append(
            {
                "source": r.get("source"),
                "record": r.get("record"),
                "cascade": cas,
                "opc_L50": opc,
                "sdd": sdd,
                "absz": absz,
                "vs_opc_gained": int(cas == 1 and opc == 0),
                "vs_opc_lost": int(cas == 0 and opc == 1),
                "vs_sdd_gained": int(cas == 1 and sdd == 0),
                "vs_sdd_lost": int(cas == 0 and sdd == 1),
                "vs_absz_gained": int(cas == 1 and absz == 0),
                "vs_absz_lost": int(cas == 0 and absz == 1),
            }
        )
    return out


def summarize_gain_loss(gl: List[dict], source: Optional[str] = None) -> dict:
    sub = [r for r in gl if source is None or r.get("source") == source]
    return {
        "n": len(sub),
        "vs_opc": {
            "gained": sum(r["vs_opc_gained"] for r in sub),
            "lost": sum(r["vs_opc_lost"] for r in sub),
            "both_on": sum(1 for r in sub if r["cascade"] and r["opc_L50"]),
            "both_off": sum(1 for r in sub if not r["cascade"] and not r["opc_L50"]),
        },
        "vs_sdd": {
            "gained": sum(r["vs_sdd_gained"] for r in sub),
            "lost": sum(r["vs_sdd_lost"] for r in sub),
            "both_on": sum(1 for r in sub if r["cascade"] and r["sdd"]),
            "both_off": sum(1 for r in sub if not r["cascade"] and not r["sdd"]),
        },
        "vs_absz": {
            "gained": sum(r["vs_absz_gained"] for r in sub),
            "lost": sum(r["vs_absz_lost"] for r in sub),
            "both_on": sum(1 for r in sub if r["cascade"] and r["absz"]),
            "both_off": sum(1 for r in sub if not r["cascade"] and not r["absz"]),
        },
    }


def build_observations(
    sens: dict,
    far: dict,
    gl_all: dict,
) -> List[str]:
    notes: List[str] = []
    c_sens = sens["cascade"]["all_events"]["sensitivity"]
    o_sens = sens["opc_L50"]["all_events"]["sensitivity"]
    s_sens = sens["sdd"]["all_events"]["sensitivity"]
    a_sens = sens["absz"]["all_events"]["sensitivity"]
    c_far = far["cascade"]["far_per_24h"]
    o_far = far["opc_L50"]["far_per_24h"]
    s_far = far["sdd"]["far_per_24h"]
    a_far = far["absz"]["far_per_24h"]

    notes.append(
        f"Cascade all-event sensitivity = {c_sens:.3f} "
        f"(OPC L=50 {o_sens:.3f}, SDD {s_sens:.3f}, abs-z {a_sens:.3f})."
    )
    notes.append(
        f"Cascade pooled NSRDB FAR = {c_far:.3f}/24h "
        f"(OPC {o_far:.3f}, SDD {s_far:.3f}, abs-z {a_far:.3f})."
    )
    if np.isfinite(c_far) and np.isfinite(s_far):
        if c_far <= s_far + 1e-9:
            notes.append(
                "Cascade FAR ≤ SDD FAR on pooled NSRDB, as expected under OPC confirmation "
                "(confirmation can only drop or keep SDD episodes, never invent SDD-less alarms)."
            )
        else:
            notes.append(
                "UNEXPECTED: cascade FAR > SDD FAR on pooled NSRDB — investigate episode "
                "counting / timebase (should not happen under pure confirm-on-SDD rule)."
            )
    notes.append(
        f"Events vs SDD: gained {gl_all['vs_sdd']['gained']}, lost {gl_all['vs_sdd']['lost']} "
        f"(lost expected when SDD fires without local OPC confirmation)."
    )
    notes.append(
        f"Events vs OPC L=50: gained {gl_all['vs_opc']['gained']}, lost {gl_all['vs_opc']['lost']} "
        f"(gains require SDD candidate co-located with OPC; losses when OPC alone was enough)."
    )
    notes.append(
        f"Events vs abs-z: gained {gl_all['vs_absz']['gained']}, lost {gl_all['vs_absz']['lost']}."
    )
    notes.append(
        "Cascade is SDD-first: it cannot detect events that SDD misses; sensitivity ceiling "
        "is SDD's sensitivity (minus unconfirmed SDD hits)."
    )
    notes.append(
        "NSRDB is healthy Holter, not device-matched to VFDB/SDDB telemetry — FAR is an "
        "interim public upper-bound estimate only."
    )
    notes.append(
        "Exploratory only: no clinical claim, no S5 claim, no superiority claim, "
        "no production threshold change."
    )
    return notes


def recommendation_text(sens: dict, far: dict, gl_all: dict) -> dict:
    c_sens = sens["cascade"]["all_events"]["sensitivity"]
    s_sens = sens["sdd"]["all_events"]["sensitivity"]
    o_sens = sens["opc_L50"]["all_events"]["sensitivity"]
    a_sens = sens["absz"]["all_events"]["sensitivity"]
    c_far = far["cascade"]["far_per_24h"]
    s_far = far["sdd"]["far_per_24h"]
    o_far = far["opc_L50"]["far_per_24h"]
    a_far = far["absz"]["far_per_24h"]

    # Simple exploratory heuristics (not clinical decision rules)
    sens_drop_vs_sdd = (s_sens - c_sens) if np.isfinite(s_sens) and np.isfinite(c_sens) else float("nan")
    far_drop_vs_sdd = (s_far - c_far) if np.isfinite(s_far) and np.isfinite(c_far) else float("nan")
    better_balance_than_absz = (
        np.isfinite(c_sens)
        and np.isfinite(a_sens)
        and np.isfinite(c_far)
        and np.isfinite(a_far)
        and c_sens >= a_sens - 0.05
        and c_far < a_far
    )
    near_opc_far = np.isfinite(c_far) and np.isfinite(o_far) and c_far <= o_far * 1.5 + 0.5

    if (
        np.isfinite(sens_drop_vs_sdd)
        and sens_drop_vs_sdd <= 0.15
        and np.isfinite(far_drop_vs_sdd)
        and far_drop_vs_sdd >= 10.0
        and (better_balance_than_absz or near_opc_far)
    ):
        warrant = "yes_further_exploration"
        text = (
            "Cascade substantially lowers FAR vs SDD while retaining most sensitivity; "
            "worth limited further exploration (e.g. secondary ±10 min note, "
            "device-matched controls) — still not production."
        )
    elif np.isfinite(c_sens) and c_sens < o_sens and np.isfinite(c_far) and c_far >= o_far:
        warrant = "low_priority"
        text = (
            "Cascade does not improve on OPC L=50 in both sens and FAR under these fixed "
            "params; further cascade work is low priority unless a new confirm rule is "
            "justified a priori (avoid post-hoc window fishing)."
        )
    elif np.isfinite(sens_drop_vs_sdd) and sens_drop_vs_sdd > 0.25:
        warrant = "limited_value_as_is"
        text = (
            "Cascade loses substantial sensitivity vs SDD (many SDD hits lack local OPC "
            "confirmation). The ±5 min SDD→OPC confirm rule as-is has limited value; "
            "do not promote. Optional: document as negative/exploratory result."
        )
    else:
        warrant = "mixed_document_only"
        text = (
            "Results are mixed: cascade sits between singleton arms without a clear "
            "Pareto improvement under fixed params. Document honestly; do not retune "
            "windows to chase a win; prefer partner-gated Tier A data over more fusion knobs."
        )

    return {
        "further_exploration_warranted": warrant,
        "text": text,
        "sens_drop_vs_sdd": sens_drop_vs_sdd,
        "far_drop_vs_sdd": far_drop_vs_sdd,
        "better_balance_than_absz_heuristic": bool(better_balance_than_absz),
    }


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(
        description="Exploratory SDD→OPC cascade fusion evaluation (±5 min confirm)"
    )
    ap.add_argument("--out-dir", default=str(RES))
    ap.add_argument(
        "--confirm-window-min",
        type=float,
        default=CONFIRM_WINDOW_MIN,
        help="Confirmation half-window in minutes (default 5)",
    )
    ap.add_argument("--refractory-h", type=float, default=REFRACTORY_H)
    ap.add_argument("--control-max-hours", type=float, default=CONTROL_MAX_HOURS)
    ap.add_argument("--skip-vfdb", action="store_true")
    ap.add_argument("--skip-nsrdb", action="store_true")
    ap.add_argument("--sddb-only-analytic", action="store_true")
    ap.add_argument("--write-doc", action="store_true")
    args = ap.parse_args(argv)

    confirm_window_h = float(args.confirm_window_min) / 60.0
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print("=== Exploratory cascade: SDD → OPC L=50 confirm ===")
    print(f"  confirm_window = ±{args.confirm_window_min} min ({confirm_window_h:.6f} h)")
    print(
        f"  OPC L={OPC_L} θ_D={OPC_THETA_D} θ_R={OPC_THETA_R} | "
        f"SDD L_c={SDD_L_C} θ_TV={SDD_THETA_TV} θ_S={SDD_THETA_S}"
    )
    print(f"  abs-z frozen z≥{FROZEN_Z_THRESHOLD} min_run={FROZEN_MIN_CONSECUTIVE}")
    print(f"  NSRDB refractory={args.refractory_h} h, cap={args.control_max_hours} h")
    print("  EXPLORATORY ONLY — no clinical / S5 / superiority claims\n")

    event_rows: List[dict] = []
    sddb_ids = list(ANALYTIC_SDDB)
    if not args.sddb_only_analytic:
        sddb_ids += EXTRA_SDDB

    print("=== SDDB ===")
    for rec in sddb_ids:
        path = DATA / f"rr_{rec}_clean.npz"
        if not path.exists():
            print(f"  skip missing {rec}")
            continue
        print(f"  SDDB {rec} ...", flush=True)
        row = process_event_record(
            "sddb", rec, path, confirm_window_h=confirm_window_h
        )
        event_rows.append(row)
        print(
            f"    cas={row.get('cascade_alarmed')} lead={row.get('cascade_lead_time_h')} | "
            f"OPC={row.get('opc_alarmed')} SDD={row.get('sdd_alarmed')} "
            f"absz={row.get('absz_alarmed')}"
        )

    if not args.skip_vfdb:
        print("\n=== VFDB ===")
        for path in sorted(RR_EXT.glob("vfdb_*_clean.npz")):
            rec = path.stem.replace("vfdb_", "").replace("_clean", "")
            print(f"  VFDB {rec} ...", flush=True)
            row = process_event_record(
                "vfdb", rec, path, confirm_window_h=confirm_window_h
            )
            event_rows.append(row)
            print(
                f"    cas={row.get('cascade_alarmed')} | OPC={row.get('opc_alarmed')} | "
                f"SDD={row.get('sdd_alarmed')} | absz={row.get('absz_alarmed')}"
            )

    control_rows: List[dict] = []
    if not args.skip_nsrdb:
        print("\n=== NSRDB FAR ===")
        paths = sorted(RR_EXT.glob("nsrdb_*_clean.npz"))
        print(f"  controls found: {len(paths)}")
        for path in paths:
            print(f"  {path.name} ...", flush=True)
            row = process_control(
                path,
                max_hours=args.control_max_hours,
                refractory_h=args.refractory_h,
                confirm_window_h=confirm_window_h,
            )
            if row.get("skipped"):
                print(f"    skipped: {row.get('skip_reason')}")
                continue
            control_rows.append(row)
            print(
                f"    FAR cas={row['cascade_far_per_24h']:.2f} "
                f"OPC={row['opc_far_per_24h']:.2f} "
                f"SDD={row['sdd_far_per_24h']:.2f} "
                f"absz={row['absz_far_per_24h']:.2f} "
                f"(ep cas/opc/sdd/absz="
                f"{int(row['cascade_n_episodes'])}/"
                f"{int(row['opc_n_episodes'])}/"
                f"{int(row['sdd_n_episodes'])}/"
                f"{int(row['absz_n_episodes'])})"
            )

    # Sensitivity blocks
    arms = ("cascade", "opc", "sdd", "absz")
    arm_key = {"cascade": "cascade", "opc": "opc", "sdd": "sdd", "absz": "absz"}
    sens: Dict[str, Any] = {}
    for label, prefix in (
        ("cascade", "cascade"),
        ("opc_L50", "opc"),
        ("sdd", "sdd"),
        ("absz", "absz"),
    ):
        sens[label] = {
            "sddb": _sens_block(event_rows, prefix, "sddb"),
            "vfdb": _sens_block(event_rows, prefix, "vfdb"),
            "all_events": _sens_block(event_rows, prefix, None),
        }

    # FAR pools
    far: Dict[str, Any] = {}
    for label, ep_key, sh_key, al_key in (
        ("cascade", "cascade_n_episodes", "cascade_search_hours", "cascade_alarmed"),
        ("opc_L50", "opc_n_episodes", "opc_search_hours", "opc_alarmed"),
        ("sdd", "sdd_n_episodes", "sdd_search_hours", "sdd_alarmed"),
        ("absz", "absz_n_episodes", "absz_search_hours", "absz_alarmed"),
    ):
        far[label] = false_alarm_rate(
            [
                {
                    "n_episodes": r[ep_key],
                    "search_hours": r[sh_key],
                    "alarmed": r[al_key],
                }
                for r in control_rows
            ]
        )

    gl = gain_loss_rows(event_rows)
    gl_all = summarize_gain_loss(gl)
    gl_sddb = summarize_gain_loss(gl, "sddb")
    gl_vfdb = summarize_gain_loss(gl, "vfdb")
    observations = build_observations(sens, far, gl_all)
    rec_block = recommendation_text(sens, far, gl_all)

    # Comparison table rows
    comp_rows = []
    for label in ("cascade", "opc_L50", "sdd", "absz"):
        s = sens[label]
        f = far[label]
        comp_rows.append(
            {
                "detector": label,
                "sddb_n": s["sddb"]["n"],
                "sddb_n_detected": s["sddb"]["n_detected"],
                "sddb_sensitivity": s["sddb"]["sensitivity"],
                "sddb_median_lead_h": s["sddb"]["median_lead_h"],
                "vfdb_n": s["vfdb"]["n"],
                "vfdb_n_detected": s["vfdb"]["n_detected"],
                "vfdb_sensitivity": s["vfdb"]["sensitivity"],
                "vfdb_median_lead_h": s["vfdb"]["median_lead_h"],
                "all_n": s["all_events"]["n"],
                "all_n_detected": s["all_events"]["n_detected"],
                "all_sensitivity": s["all_events"]["sensitivity"],
                "all_median_lead_h": s["all_events"]["median_lead_h"],
                "nsrdb_n": f.get("n_controls", 0),
                "nsrdb_total_episodes": f.get("total_episodes", float("nan")),
                "nsrdb_search_hours": f.get("total_search_hours", float("nan")),
                "nsrdb_far_per_24h": f.get("far_per_24h", float("nan")),
                "nsrdb_fraction_alarmed": f.get("fraction_alarmed", float("nan")),
            }
        )

    write_csv(out_dir / "ordinal_cascade_per_record.csv", event_rows)
    write_csv(out_dir / "ordinal_cascade_nsrdb_far_per_record.csv", control_rows)
    write_csv(out_dir / "ordinal_cascade_comparison.csv", comp_rows)
    write_csv(out_dir / "ordinal_cascade_gain_loss.csv", gl)

    summary = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "phase": "ordinal_cascade_sdd_confirm_opc",
        "exploratory_only": True,
        "clinical_claim": False,
        "superiority_claim": False,
        "s5_claim": False,
        "fusion": True,
        "fusion_kind": "cascade_sdd_confirm_opc",
        "rule": {
            "description": (
                "SDD is the candidate filter; SDD alarm at t becomes cascade alarm "
                "only if OPC L=50 also alarms within closed window [t-W, t+W]. "
                "Event evaluation is causal: only OPC with t_opc < event_hr may confirm; "
                "decision_time = max(t_SDD, t_OPC_confirm) must be < event_hr "
                "(no post-event look-ahead)."
            ),
            "confirm_window_min": float(args.confirm_window_min),
            "confirm_window_h": confirm_window_h,
            "causal": True,
            "opc_max_hr": "event_hr (strictly before event)",
            "decision_time": "max(t_SDD, earliest pre-event OPC in ±W)",
            "window_justification": (
                "±5 min is a local co-occurrence confirm: short VFDB horizons and "
                "multi-hour SDDB leads both only need same short-term instability epoch; "
                "not a second independent long-horizon detector. "
                "Causality forbids post-event OPC from confirming pre-event SDD."
            ),
            "candidate": "SDD",
            "confirm": "OPC_L50",
        },
        "params": {
            "opc": {
                "L": OPC_L,
                "theta_D": OPC_THETA_D,
                "theta_R": OPC_THETA_R,
                "K": K_JOINT,
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
            "encoding": "joint_bivariate_BP_m3",
            "singleton_detectors_retuned": False,
        },
        "methodology": {
            "sensitivity": (
                "first causal cascade decision in post-basal pre-event search; "
                "decision_time = max(t_SDD, t_OPC_pre_event); lead = event_hr - decision_time; "
                "OPC with t >= event_hr cannot confirm"
            ),
            "far_formula": "FAR = total_episodes / total_search_hours * 24",
            "far_cascade_note": (
                "NSRDB controls have no event; cascade FAR uses bidirectional ±W on full "
                "control streams (no post-event look-ahead issue)."
            ),
            "refractory_h": args.refractory_h,
            "control_max_hours": args.control_max_hours,
            "basal_hours": BASAL_HOURS,
            "device_mismatch": True,
            "device_mismatch_note": (
                "NSRDB is rhythm-healthy Holter ECG — NOT device-matched to VFDB/SDDB. "
                "FAR is an interim public upper-bound estimate only."
            ),
        },
        "n_event_records": len(event_rows),
        "n_sddb": sum(1 for r in event_rows if r.get("source") == "sddb"),
        "n_vfdb": sum(1 for r in event_rows if r.get("source") == "vfdb"),
        "n_nsrdb_controls": len(control_rows),
        "sensitivity": sens,
        "pooled_far": far,
        "gain_loss": {
            "all_events": gl_all,
            "sddb": gl_sddb,
            "vfdb": gl_vfdb,
        },
        "observations": observations,
        "recommendation": rec_block,
        "limitations": [
            "Cascade cannot recover events that SDD misses (SDD is hard filter).",
            "Confirmation may drop true SDD positives when OPC does not co-fire locally pre-event.",
            "Causal rule: post-event OPC cannot confirm (no look-ahead); sensitivity can drop vs non-causal co-occurrence.",
            "VFDB short records: ±5 min window can miss staggered OPC/SDD peaks.",
            "NSRDB device mismatch; FAR not ICU/telemetry-matched.",
            "Timebases: abs-z on strided τ_s vs symbol endpoints for ordinal/cascade.",
            "Fixed params only; no threshold or window grid search (avoids fishing).",
            "Exploratory fusion arm — singleton experiments remain fusion=False.",
        ],
        "executive_summary": (
            f"Light cascade SDD→OPC (±{args.confirm_window_min:g} min): "
            f"all-event sens={sens['cascade']['all_events']['sensitivity']:.3f} "
            f"(SDDB {sens['cascade']['sddb']['sensitivity']:.3f}, "
            f"VFDB {sens['cascade']['vfdb']['sensitivity']:.3f}); "
            f"NSRDB FAR={far['cascade']['far_per_24h']:.3f}/24h "
            f"vs OPC {far['opc_L50']['far_per_24h']:.3f}, "
            f"SDD {far['sdd']['far_per_24h']:.3f}, "
            f"abs-z {far['absz']['far_per_24h']:.3f}. "
            f"Recommendation: {rec_block['further_exploration_warranted']}. "
            "Exploratory only — no clinical claims."
        ),
    }

    with open(out_dir / "ordinal_cascade_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    # Console report
    print("\n=== COMPARISON (sensitivity + FAR) ===")
    print(
        f"{'Detector':<12} {'SDDB sens':>10} {'VFDB sens':>10} {'All sens':>10} "
        f"{'FAR/24h':>10} {'Episodes':>10}"
    )
    for row in comp_rows:
        print(
            f"{row['detector']:<12} "
            f"{row['sddb_sensitivity']:>10.3f} "
            f"{row['vfdb_sensitivity']:>10.3f} "
            f"{row['all_sensitivity']:>10.3f} "
            f"{row['nsrdb_far_per_24h']:>10.3f} "
            f"{row['nsrdb_total_episodes']:>10.0f}"
        )

    print("\n=== GAIN / LOSS (cascade vs singletons, all events) ===")
    for vs in ("vs_opc", "vs_sdd", "vs_absz"):
        block = gl_all[vs]
        print(
            f"  {vs}: gained={block['gained']} lost={block['lost']} "
            f"both_on={block['both_on']} both_off={block['both_off']}"
        )

    print("\n=== OBSERVATIONS ===")
    for n in observations:
        print(f"  - {n}")

    print("\n=== EXECUTIVE SUMMARY ===")
    print(summary["executive_summary"])
    print(f"\nRecommendation: {rec_block['further_exploration_warranted']}")
    print(f"  {rec_block['text']}")

    print(f"\nWrote artifacts under {out_dir}")
    for name in (
        "ordinal_cascade_per_record.csv",
        "ordinal_cascade_nsrdb_far_per_record.csv",
        "ordinal_cascade_comparison.csv",
        "ordinal_cascade_gain_loss.csv",
        "ordinal_cascade_summary.json",
    ):
        print(f"  {out_dir / name}")

    if args.write_doc:
        write_doc(summary, comp_rows, gl_all, out_dir)

    return 0


def write_doc(
    summary: dict,
    comp_rows: List[dict],
    gl_all: dict,
    out_dir: Path,
) -> None:
    doc = BASE / "docs" / "ORDINAL_CASCADE_FUSION.md"
    c = next(r for r in comp_rows if r["detector"] == "cascade")
    o = next(r for r in comp_rows if r["detector"] == "opc_L50")
    s = next(r for r in comp_rows if r["detector"] == "sdd")
    a = next(r for r in comp_rows if r["detector"] == "absz")
    rule = summary["rule"]
    rec = summary["recommendation"]

    lines = [
        "# Exploratory light cascade: SDD → OPC L=50 (±5 min confirm)",
        "",
        f"**Generated:** {summary['generated_at']}",
        "",
        "**Status:** Exploratory only. No clinical claim, no S5 claim, "
        "no superiority claim, no production threshold change.",
        "",
        "## Cascade rule",
        "",
        "1. **Candidate filter:** SDD "
        f"(L_c={summary['params']['sdd']['L_c']}, "
        f"θ_TV={summary['params']['sdd']['theta_TV']}, "
        f"θ_S={summary['params']['sdd']['theta_S']}).",
        "2. **Confirmation:** An SDD alarm at time *t* becomes a cascade alarm only if "
        f"**OPC L={summary['params']['opc']['L']}** "
        f"(θ_D={summary['params']['opc']['theta_D']}, "
        f"θ_R={summary['params']['opc']['theta_R']}) "
        f"also alarms within a **closed window [t−W, t+W]** with "
        f"**W = ±{rule['confirm_window_min']:g} minutes** "
        f"({rule['confirm_window_h']:.6f} h).",
        "3. Cascade marks only confirmed **SDD-on** samples (OPC-only never starts a cascade alarm).",
        "4. **Causality (events):** only OPC with `t_opc < event_hr` may confirm; "
        "`decision_time = max(t_SDD, t_OPC_confirm) < event_hr` "
        "(no post-event look-ahead).",
        "",
        "### Window justification",
        "",
        rule["window_justification"],
        "",
        "## Head-to-head comparison",
        "",
        "| Detector | SDDB sens (n=11) | VFDB sens (n=22) | All sens (n=33) | NSRDB FAR (/24h) | Episodes |",
        "|----------|------------------|------------------|-----------------|------------------|----------|",
        f"| **Cascade SDD→OPC** | {c['sddb_sensitivity']:.3f} ({int(c['sddb_n_detected'])}/11) | "
        f"{c['vfdb_sensitivity']:.3f} ({int(c['vfdb_n_detected'])}/22) | "
        f"**{c['all_sensitivity']:.3f}** ({int(c['all_n_detected'])}/33) | "
        f"**{c['nsrdb_far_per_24h']:.3f}** | {int(c['nsrdb_total_episodes'])} |",
        f"| OPC L=50 alone | {o['sddb_sensitivity']:.3f} ({int(o['sddb_n_detected'])}/11) | "
        f"{o['vfdb_sensitivity']:.3f} ({int(o['vfdb_n_detected'])}/22) | "
        f"{o['all_sensitivity']:.3f} ({int(o['all_n_detected'])}/33) | "
        f"{o['nsrdb_far_per_24h']:.3f} | {int(o['nsrdb_total_episodes'])} |",
        f"| SDD alone | {s['sddb_sensitivity']:.3f} ({int(s['sddb_n_detected'])}/11) | "
        f"{s['vfdb_sensitivity']:.3f} ({int(s['vfdb_n_detected'])}/22) | "
        f"{s['all_sensitivity']:.3f} ({int(s['all_n_detected'])}/33) | "
        f"{s['nsrdb_far_per_24h']:.3f} | {int(s['nsrdb_total_episodes'])} |",
        f"| abs-z τ_s frozen | {a['sddb_sensitivity']:.3f} ({int(a['sddb_n_detected'])}/11) | "
        f"{a['vfdb_sensitivity']:.3f} ({int(a['vfdb_n_detected'])}/22) | "
        f"{a['all_sensitivity']:.3f} ({int(a['all_n_detected'])}/33) | "
        f"{a['nsrdb_far_per_24h']:.3f} | {int(a['nsrdb_total_episodes'])} |",
        "",
        "## Events gained / lost (cascade vs singletons)",
        "",
        f"| Comparison | Gained | Lost | Both on | Both off |",
        f"|------------|--------|------|---------|----------|",
        f"| vs OPC L=50 | {gl_all['vs_opc']['gained']} | {gl_all['vs_opc']['lost']} | "
        f"{gl_all['vs_opc']['both_on']} | {gl_all['vs_opc']['both_off']} |",
        f"| vs SDD | {gl_all['vs_sdd']['gained']} | {gl_all['vs_sdd']['lost']} | "
        f"{gl_all['vs_sdd']['both_on']} | {gl_all['vs_sdd']['both_off']} |",
        f"| vs abs-z | {gl_all['vs_absz']['gained']} | {gl_all['vs_absz']['lost']} | "
        f"{gl_all['vs_absz']['both_on']} | {gl_all['vs_absz']['both_off']} |",
        "",
        "## Observations",
        "",
    ]
    for n in summary["observations"]:
        lines.append(f"- {n}")
    lines += [
        "",
        "## Limitations",
        "",
    ]
    for lim in summary["limitations"]:
        lines.append(f"- {lim}")
    lines += [
        "",
        "## Executive summary",
        "",
        summary["executive_summary"],
        "",
        "## Recommendation",
        "",
        f"**{rec['further_exploration_warranted']}:** {rec['text']}",
        "",
        "## Reproduce",
        "",
        "```bash",
        "cd Investigaciones/Cardiac_CCTP_Pilot",
        "python3 code/run_ordinal_cascade_fusion.py --write-doc",
        "python3 -m pytest tests/test_ordinal_cascade_fusion.py -q",
        "```",
        "",
        "## Artifacts",
        "",
        "- `results/ordinal_cascade_comparison.csv`",
        "- `results/ordinal_cascade_per_record.csv`",
        "- `results/ordinal_cascade_nsrdb_far_per_record.csv`",
        "- `results/ordinal_cascade_gain_loss.csv`",
        "- `results/ordinal_cascade_summary.json`",
        "- Pure merger: `code/ordinal_detectors/cascade_fusion.py`",
        "- Entry: `code/run_ordinal_cascade_fusion.py`",
        "",
    ]
    doc.write_text("\n".join(lines) + "\n")
    print(f"Wrote {doc}")


if __name__ == "__main__":
    raise SystemExit(main())
