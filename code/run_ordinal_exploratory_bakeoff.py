#!/usr/bin/env python3
"""
Exploratory bake-off: OPC (Option 1) and SDD (Option 2) separately vs frozen abs-z.

- Builds 1-D ordinal symbol streams from cleaned RR (no pre-baked symbol NPZ).
- Runs shipped opc_detect and sdd_detect independently (no fusion).
- Compares to abs-z on continuous τ_s via detect_lead_time (frozen z=2, min_run=3).
- Does NOT retune production abs-z / CCTP thresholds.

Outputs under results/:
  ordinal_data_inventory.txt
  ordinal_opc_per_record.csv
  ordinal_sdd_per_record.csv
  ordinal_absz_per_record.csv
  ordinal_opc_sdd_absz_comparison.csv
  ordinal_exploratory_summary.json

Write-up: docs/ORDINAL_EXPLORATORY_BAKEOFF.md (generated separately or by --write-doc).

Exploratory only — no clinical / superiority claims.
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
    FROZEN_MIN_CONSECUTIVE,
    FROZEN_Z_THRESHOLD,
    STRIDE,
    W_TAU,
    build_bivariate_proxy,
    detect_lead_time,
    get_event_and_windows,
    short_db_windows,
)
from ordinal_detectors.opc_detector import opc_detect, opc_first_alarm_index
from ordinal_detectors.sdd_detector import sdd_detect, sdd_first_alarm_index
from recd_ordinal_levels import generate_multivariate_symbols

BASE = Path(__file__).resolve().parent.parent
DATA = BASE / "data"
RR_EXT = DATA / "rr_external"
RES = BASE / "results"
RES.mkdir(parents=True, exist_ok=True)

# SDDB analytic set (+ optional 44)
ANALYTIC_SDDB = ["30", "31", "32", "35", "36", "38", "45", "47", "50", "51"]
EXTRA_SDDB = ["44"]

# Symbol encoding (frozen for this bake-off)
M_EMB = 3
DELAY = 1
BP_ALPHABET = 6  # m!
K_JOINT = BP_ALPHABET * BP_ALPHABET  # 36 joint bivariate codes

# Suggested starting params (docs/ORDINAL_ALARM_DETECTORS.md)
OPC_L = 8
OPC_THETA_D = 0.35
OPC_THETA_R = 5
# K-aware companion L so support diversity can exceed θ_D when K=36:
# max D for window L is min(L,K)/K; need L > θ_D * K ≈ 12.6 → L>=13.
# Use L=50 (same order as SDD L_c) as a non-optimized scale companion, not a retune of abs-z.
OPC_L_SCALE = 50

SDD_L_C = 50
SDD_THETA_TV = 0.35
SDD_THETA_S = 1


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
    """
    1-D joint Bandt–Pompe symbols from project bivariate proxy.

    S[t,0], S[t,1] ∈ {0,...,5}; sigma = S0 * 6 + S1 ∈ {0,...,35}; K=36.
    Returns (sigma, K, offset) where offset maps symbol index → beat index.
    """
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
    """First symbol index with t_hr_sym >= hr (clip to range)."""
    if len(t_hr_sym) == 0:
        return 0
    idx = int(np.searchsorted(t_hr_sym, hr, side="left"))
    return max(0, min(idx, len(t_hr_sym)))


def hours_to_symbol_index_right(t_hr_sym: np.ndarray, hr: float) -> int:
    """Last exclusive index with t_hr_sym < hr."""
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


def basal_alarm_stats(
    alarm: np.ndarray,
    t_hr_sym: np.ndarray,
    basal: Tuple[float, float],
) -> Dict[str, float]:
    b0, b1 = basal
    mask = (t_hr_sym >= b0) & (t_hr_sym < b1)
    n = int(mask.sum())
    if n == 0:
        return {"basal_n": 0.0, "basal_alarm_count": 0.0, "basal_alarm_frac": float("nan")}
    c = int(np.sum(alarm[mask] > 0))
    return {
        "basal_n": float(n),
        "basal_alarm_count": float(c),
        "basal_alarm_frac": float(c / n),
    }


def inventory_files() -> List[str]:
    lines: List[str] = []
    lines.append("ORDINAL EXPLORATORY BAKE-OFF — DATA INVENTORY")
    lines.append(f"generated_at: {datetime.now(timezone.utc).isoformat()}")
    lines.append("")
    lines.append("NOTE: Repo stores cleaned RR NPZ, not pre-baked ordinal-symbol NPZ.")
    lines.append("Symbols are derived at runtime: bivariate proxy → Bandt–Pompe m=3,τ=1")
    lines.append(f"joint code σ = s0*{BP_ALPHABET}+s1, K={K_JOINT}.")
    lines.append("")
    lines.append("=== SDDB cleaned RR (analytic + optional) ===")
    for rec in ANALYTIC_SDDB + EXTRA_SDDB:
        p = DATA / f"rr_{rec}_clean.npz"
        if p.exists():
            d = load_npz(p)
            lines.append(
                f"  OK  {p.name}  n_beats={int(d['n_beats'])}  "
                f"total_h={float(d['total_hours']):.3f}  vfon_h={float(d['vfon_sec'])/3600:.3f}"
            )
        else:
            lines.append(f"  MISSING  {p.name}")
    lines.append("")
    lines.append("=== VFDB cleaned RR (data/rr_external/vfdb_*_clean.npz) ===")
    vf = sorted(RR_EXT.glob("vfdb_*_clean.npz"))
    lines.append(f"  count={len(vf)}")
    for p in vf:
        d = load_npz(p)
        lines.append(
            f"  OK  {p.name}  n_beats={int(d['n_beats'])}  total_h={float(d['total_hours']):.3f}  "
            f"pre_event_h={float(d.get('pre_event_hours', float(d['vfon_sec'])/3600)):.3f}  "
            f"label={d.get('event_label', '')}"
        )
    lines.append("")
    lines.append("=== Raw PhysioNet dirs (NOT primary symbol store) ===")
    for name in ("sddb", "vfdb"):
        d = DATA / name
        if d.is_dir():
            n = len(list(d.iterdir()))
            lines.append(f"  {name}/  entries={n}  (hea/atr/dat annotations)")
        else:
            lines.append(f"  {name}/  MISSING")
    lines.append("")
    lines.append("=== Optional NSRDB controls ===")
    ns = sorted(RR_EXT.glob("nsrdb_*_clean.npz"))
    lines.append(f"  count={len(ns)}")
    for p in ns[:5]:
        lines.append(f"  OK  {p.name}")
    if len(ns) > 5:
        lines.append(f"  ... +{len(ns)-5} more")
    lines.append("")
    lines.append("=== Pre-baked ordinal symbol NPZ search ===")
    hits = list(DATA.rglob("*symbol*")) + list(DATA.rglob("*ordinal*"))
    if hits:
        for h in hits[:20]:
            lines.append(f"  found {h.relative_to(DATA)}")
    else:
        lines.append("  none (expected — symbols derived from RR)")
    return lines


def process_event_record(
    source: str,
    record: str,
    path: Path,
    opc_L: int,
    opc_theta_D: float,
    opc_theta_R: int,
    sdd_L_c: int,
    sdd_theta_TV: float,
    sdd_theta_S: int,
) -> Tuple[dict, dict, dict]:
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

    sigma, K, offset = joint_bivariate_symbols(rr)
    if len(sigma) == 0:
        empty = {
            "source": source,
            "record": record,
            "path": str(path.name),
            "error": "empty_symbol_stream",
            "alarmed": 0,
            "lead_time_h": float("nan"),
        }
        return empty, empty, empty

    # Align symbol times to beat times
    t_hr_full = t_hr
    t_hr_sym = t_hr_full[offset : offset + len(sigma)]
    if len(t_hr_sym) < len(sigma):
        sigma = sigma[: len(t_hr_sym)]
    elif len(t_hr_sym) > len(sigma):
        t_hr_sym = t_hr_sym[: len(sigma)]

    b0, b1 = basal
    a0, a1 = approach
    # Basal indices on symbol stream [start, end)
    basal_i0 = hours_to_symbol_index(t_hr_sym, b0)
    basal_i1 = hours_to_symbol_index_right(t_hr_sym, b1)
    if basal_i1 <= basal_i0:
        # fallback: first 10% of pre-event symbols
        event_i = hours_to_symbol_index_right(t_hr_sym, event_hr)
        basal_i0 = 0
        basal_i1 = max(1, event_i // 3)

    search_start = hours_to_symbol_index(t_hr_sym, b1)  # after basal
    search_end = hours_to_symbol_index_right(t_hr_sym, event_hr)

    # --- OPC (independent) ---
    opc_idx, opc_out = opc_first_alarm_index(
        sigma,
        L=opc_L,
        theta_D=opc_theta_D,
        theta_R=opc_theta_R,
        K=K,
        search_start=search_start,
        search_end=search_end,
    )
    opc_lead, opc_det_hr, opc_alarmed = lead_from_symbol_alarm(t_hr_sym, opc_idx, event_hr)
    opc_basal = basal_alarm_stats(opc_out["alarm"], t_hr_sym, basal)
    # diversity diagnostics in approach
    app_mask = (t_hr_sym >= a0) & (t_hr_sym < a1)
    div = opc_out["diversity"]
    div_app = div[app_mask & np.isfinite(div)]
    opc_row = {
        "source": source,
        "record": record,
        "event_label": event_label,
        "duration_stratum": stratum,
        "event_hr": event_hr,
        "basal_start": b0,
        "basal_end": b1,
        "approach_start": a0,
        "approach_end": a1,
        "n_symbols": int(len(sigma)),
        "K": K,
        "encoding": "joint_bivariate_BP_m3",
        "L": opc_L,
        "theta_D": opc_theta_D,
        "theta_R": opc_theta_R,
        "alarmed": opc_alarmed,
        "lead_time_h": opc_lead,
        "detection_hr": opc_det_hr,
        "first_alarm_idx": opc_idx if opc_idx is not None else "",
        "search_start_idx": search_start,
        "search_end_idx": search_end,
        "mean_D_approach": float(np.mean(div_app)) if len(div_app) else float("nan"),
        "min_D_approach": float(np.min(div_app)) if len(div_app) else float("nan"),
        "frac_low_div_approach": float(np.mean(opc_out["low_div"][app_mask])) if app_mask.any() else float("nan"),
        **opc_basal,
    }

    # --- SDD (independent, TV) ---
    basal_range = (basal_i0, basal_i1)
    if basal_range[1] - basal_range[0] < 2:
        sdd_row = {
            "source": source,
            "record": record,
            "error": "basal_too_short_for_sdd",
            "alarmed": 0,
            "lead_time_h": float("nan"),
            "basal_start": b0,
            "basal_end": b1,
        }
    else:
        sdd_idx, sdd_out = sdd_first_alarm_index(
            sigma,
            basal_range,
            L_c=sdd_L_c,
            theta_TV=sdd_theta_TV,
            theta_S=sdd_theta_S,
            K=K,
            search_start=search_start,
            search_end=search_end,
            mask_basal=True,
        )
        sdd_lead, sdd_det_hr, sdd_alarmed = lead_from_symbol_alarm(
            t_hr_sym, sdd_idx, event_hr
        )
        # Basal diagnostic: high_tv fraction with mask_basal=False (not production alarm)
        sdd_diag = sdd_detect(
            sigma,
            basal_range,
            L_c=sdd_L_c,
            theta_TV=sdd_theta_TV,
            theta_S=sdd_theta_S,
            K=K,
            mask_basal=False,
        )
        tv = sdd_out["tv"]
        tv_app = tv[app_mask & np.isfinite(tv)]
        # basal high_tv diagnostic on non-masked run
        bmask = (t_hr_sym >= b0) & (t_hr_sym < b1)
        ht = sdd_diag["high_tv"]
        n_b = int(bmask.sum())
        c_b = int(np.sum(ht[bmask] > 0)) if n_b else 0
        sdd_row = {
            "source": source,
            "record": record,
            "event_label": event_label,
            "duration_stratum": stratum,
            "event_hr": event_hr,
            "basal_start": b0,
            "basal_end": b1,
            "approach_start": a0,
            "approach_end": a1,
            "n_symbols": int(len(sigma)),
            "K": K,
            "encoding": "joint_bivariate_BP_m3",
            "L_c": sdd_L_c,
            "theta_TV": sdd_theta_TV,
            "theta_S": sdd_theta_S,
            "basal_idx_start": basal_i0,
            "basal_idx_end": basal_i1,
            "basal_n_symbols": basal_i1 - basal_i0,
            "alarmed": sdd_alarmed,
            "lead_time_h": sdd_lead,
            "detection_hr": sdd_det_hr,
            "first_alarm_idx": sdd_idx if sdd_idx is not None else "",
            "search_start_idx": search_start,
            "search_end_idx": search_end,
            "mean_TV_approach": float(np.mean(tv_app)) if len(tv_app) else float("nan"),
            "max_TV_approach": float(np.max(tv_app)) if len(tv_app) else float("nan"),
            "basal_high_tv_count_diag": float(c_b),
            "basal_n": float(n_b),
            "basal_high_tv_frac_diag": float(c_b / n_b) if n_b else float("nan"),
            "note": "primary_alarm_mask_basal=True; basal_high_tv is diagnostic only",
        }

    # --- abs-z on τ_s (frozen) ---
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
    absz_row = {
        "source": source,
        "record": record,
        "event_label": event_label,
        "duration_stratum": stratum,
        "event_hr": event_hr,
        "basal_start": b0,
        "basal_end": b1,
        "approach_start": a0,
        "approach_end": a1,
        "metric": "tau_s",
        "z_threshold": FROZEN_Z_THRESHOLD,
        "min_consecutive": FROZEN_MIN_CONSECUTIVE,
        "alarmed": int(absz["alarmed"]),
        "lead_time_h": absz["lead_time_h"],
        "detection_hr": absz["detection_hr"],
        "z_at_detection": absz["z_at_detection"],
        "n_search_points": absz["n_search_points"],
    }
    return opc_row, sdd_row, absz_row


def process_control_nsrdb(
    path: Path,
    opc_L: int,
    opc_theta_D: float,
    opc_theta_R: int,
    sdd_L_c: int,
    sdd_theta_TV: float,
    sdd_theta_S: int,
    max_hours: float = 6.0,
) -> Tuple[dict, dict]:
    """Control flavor: alarm fraction on early segment (no VF event)."""
    d = load_npz(path)
    rec = str(d.get("record_id", path.stem.replace("nsrdb_", "").replace("_clean", "")))
    rr = np.asarray(d["rr_ms"], dtype=float)
    t_sec = np.asarray(d["t_sec"], dtype=float)
    t_hr = t_sec / 3600.0
    # restrict for speed
    keep = t_hr <= max_hours
    if keep.sum() < 500:
        keep = np.ones(len(rr), dtype=bool)
    rr = rr[keep]
    t_hr = t_hr[keep]
    sigma, K, offset = joint_bivariate_symbols(rr)
    t_hr_sym = t_hr[offset : offset + len(sigma)]
    if len(t_hr_sym) < len(sigma):
        sigma = sigma[: len(t_hr_sym)]
    else:
        t_hr_sym = t_hr_sym[: len(sigma)]
    if len(sigma) < 200:
        return {"record": rec, "error": "too_short"}, {"record": rec, "error": "too_short"}

    # basal-like: first hour; search: hours 1..max_hours
    b0, b1 = 0.1, min(1.1, float(t_hr_sym[-1]) * 0.3)
    basal_i0 = hours_to_symbol_index(t_hr_sym, b0)
    basal_i1 = hours_to_symbol_index_right(t_hr_sym, b1)
    if basal_i1 - basal_i0 < 50:
        basal_i0, basal_i1 = 0, min(200, len(sigma) // 5)

    opc_out = opc_detect(sigma, L=opc_L, theta_D=opc_theta_D, theta_R=opc_theta_R, K=K)
    # search after basal
    s0 = basal_i1
    opc_alarm_frac = float(np.mean(opc_out["alarm"][s0:] > 0)) if s0 < len(sigma) else float("nan")
    opc_row = {
        "source": "nsrdb",
        "record": rec,
        "role": "control",
        "n_symbols": len(sigma),
        "L": opc_L,
        "theta_D": opc_theta_D,
        "theta_R": opc_theta_R,
        "post_basal_alarm_frac": opc_alarm_frac,
        "post_basal_alarm_count": int(np.sum(opc_out["alarm"][s0:] > 0)),
        "post_basal_n": int(len(sigma) - s0),
    }

    sdd_out = sdd_detect(
        sigma,
        (basal_i0, basal_i1),
        L_c=sdd_L_c,
        theta_TV=sdd_theta_TV,
        theta_S=sdd_theta_S,
        K=K,
        mask_basal=True,
    )
    sdd_alarm_frac = float(np.mean(sdd_out["alarm"][s0:] > 0)) if s0 < len(sigma) else float("nan")
    sdd_row = {
        "source": "nsrdb",
        "record": rec,
        "role": "control",
        "n_symbols": len(sigma),
        "L_c": sdd_L_c,
        "theta_TV": sdd_theta_TV,
        "theta_S": sdd_theta_S,
        "post_basal_alarm_frac": sdd_alarm_frac,
        "post_basal_alarm_count": int(np.sum(sdd_out["alarm"][s0:] > 0)),
        "post_basal_n": int(len(sigma) - s0),
    }
    return opc_row, sdd_row


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
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def summarize(opc_rows, sdd_rows, absz_rows, opc_scale_rows=None) -> dict:
    def sens(rows, source=None):
        sub = [r for r in rows if "error" not in r or not r.get("error")]
        if source:
            sub = [r for r in sub if r.get("source") == source]
        n = len(sub)
        det = sum(1 for r in sub if int(r.get("alarmed", 0)) == 1)
        leads = [float(r["lead_time_h"]) for r in sub if int(r.get("alarmed", 0)) == 1 and np.isfinite(float(r.get("lead_time_h", float("nan"))))]
        return {
            "n": n,
            "n_detected": det,
            "sensitivity": float(det / n) if n else float("nan"),
            "median_lead_h": float(np.median(leads)) if leads else float("nan"),
            "mean_lead_h": float(np.mean(leads)) if leads else float("nan"),
        }

    out = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "params": {
            "encoding": "joint_bivariate_BP_m3",
            "K": K_JOINT,
            "m": M_EMB,
            "delay": DELAY,
            "opc": {"L": OPC_L, "theta_D": OPC_THETA_D, "theta_R": OPC_THETA_R},
            "opc_scale_companion": {
                "L": OPC_L_SCALE,
                "theta_D": OPC_THETA_D,
                "theta_R": OPC_THETA_R,
                "note": "L expanded so min(L,K)/K can exceed theta_D for K=36; not abs-z retune",
            },
            "sdd": {"L_c": SDD_L_C, "theta_TV": SDD_THETA_TV, "theta_S": SDD_THETA_S},
            "absz": {
                "metric": "tau_s",
                "z_threshold": FROZEN_Z_THRESHOLD,
                "min_consecutive": FROZEN_MIN_CONSECUTIVE,
            },
            "fusion": False,
        },
        "opc_default": {
            "sddb": sens(opc_rows, "sddb"),
            "vfdb": sens(opc_rows, "vfdb"),
            "all_events": sens(opc_rows),
        },
        "sdd": {
            "sddb": sens(sdd_rows, "sddb"),
            "vfdb": sens(sdd_rows, "vfdb"),
            "all_events": sens(sdd_rows),
        },
        "absz_tau_s": {
            "sddb": sens(absz_rows, "sddb"),
            "vfdb": sens(absz_rows, "vfdb"),
            "all_events": sens(absz_rows),
        },
        "clinical_claim": False,
        "superiority_claim": False,
        "exploratory_only": True,
    }
    if opc_scale_rows is not None:
        out["opc_L50_companion"] = {
            "sddb": sens(opc_scale_rows, "sddb"),
            "vfdb": sens(opc_scale_rows, "vfdb"),
            "all_events": sens(opc_scale_rows),
        }
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description="OPC/SDD exploratory bake-off vs abs-z")
    ap.add_argument("--skip-nsrdb", action="store_true")
    ap.add_argument("--skip-vfdb", action="store_true")
    ap.add_argument("--sddb-only-analytic", action="store_true", help="Skip record 44")
    ap.add_argument("--out-dir", default=str(RES))
    args = ap.parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    inv_lines = inventory_files()
    inv_path = out_dir / "ordinal_data_inventory.txt"
    inv_path.write_text("\n".join(inv_lines) + "\n")
    print("\n".join(inv_lines[:40]))
    print(f"... wrote {inv_path}")

    sddb_ids = list(ANALYTIC_SDDB)
    if not args.sddb_only_analytic:
        sddb_ids += EXTRA_SDDB

    opc_rows: List[dict] = []
    opc_scale_rows: List[dict] = []
    sdd_rows: List[dict] = []
    absz_rows: List[dict] = []

    print("\n=== SDDB ===")
    for rec in sddb_ids:
        path = DATA / f"rr_{rec}_clean.npz"
        if not path.exists():
            print(f"  skip missing {rec}")
            continue
        print(f"  SDDB {rec} ...", flush=True)
        opc, sdd, absz = process_event_record(
            "sddb",
            rec,
            path,
            OPC_L,
            OPC_THETA_D,
            OPC_THETA_R,
            SDD_L_C,
            SDD_THETA_TV,
            SDD_THETA_S,
        )
        opc_rows.append(opc)
        sdd_rows.append(sdd)
        absz_rows.append(absz)
        opc_s, _, _ = process_event_record(
            "sddb",
            rec,
            path,
            OPC_L_SCALE,
            OPC_THETA_D,
            OPC_THETA_R,
            SDD_L_C,
            SDD_THETA_TV,
            SDD_THETA_S,
        )
        opc_s["L_note"] = "scale_companion_L50"
        opc_scale_rows.append(opc_s)
        print(
            f"    OPC L={OPC_L} alarmed={opc.get('alarmed')} lead={opc.get('lead_time_h')} | "
            f"OPC L50 alarmed={opc_s.get('alarmed')} | "
            f"SDD alarmed={sdd.get('alarmed')} lead={sdd.get('lead_time_h')} | "
            f"absz alarmed={absz.get('alarmed')} lead={absz.get('lead_time_h')}"
        )

    if not args.skip_vfdb:
        print("\n=== VFDB ===")
        for path in sorted(RR_EXT.glob("vfdb_*_clean.npz")):
            rec = path.stem.replace("vfdb_", "").replace("_clean", "")
            d = load_npz(path)
            pre = float(d.get("pre_event_hours", float(d["vfon_sec"]) / 3600.0))
            # Process all with short windows; mark stratum
            print(f"  VFDB {rec} pre_h={pre:.3f} ...", flush=True)
            opc, sdd, absz = process_event_record(
                "vfdb",
                rec,
                path,
                OPC_L,
                OPC_THETA_D,
                OPC_THETA_R,
                SDD_L_C,
                SDD_THETA_TV,
                SDD_THETA_S,
            )
            opc_rows.append(opc)
            sdd_rows.append(sdd)
            absz_rows.append(absz)
            opc_s, _, _ = process_event_record(
                "vfdb",
                rec,
                path,
                OPC_L_SCALE,
                OPC_THETA_D,
                OPC_THETA_R,
                SDD_L_C,
                SDD_THETA_TV,
                SDD_THETA_S,
            )
            opc_s["L_note"] = "scale_companion_L50"
            opc_scale_rows.append(opc_s)
            print(
                f"    OPC alarmed={opc.get('alarmed')} | SDD alarmed={sdd.get('alarmed')} | "
                f"absz alarmed={absz.get('alarmed')} stratum={opc.get('duration_stratum')}"
            )

    nsrdb_opc: List[dict] = []
    nsrdb_sdd: List[dict] = []
    if not args.skip_nsrdb:
        print("\n=== NSRDB control flavor (optional) ===")
        for path in sorted(RR_EXT.glob("nsrdb_*_clean.npz"))[:8]:
            print(f"  {path.name} ...", flush=True)
            o, s = process_control_nsrdb(
                path, OPC_L, OPC_THETA_D, OPC_THETA_R, SDD_L_C, SDD_THETA_TV, SDD_THETA_S
            )
            nsrdb_opc.append(o)
            nsrdb_sdd.append(s)

    # Comparison table
    comp: List[dict] = []
    key = lambda r: (r.get("source"), str(r.get("record")))
    opc_map = {key(r): r for r in opc_rows}
    opc_s_map = {key(r): r for r in opc_scale_rows}
    sdd_map = {key(r): r for r in sdd_rows}
    absz_map = {key(r): r for r in absz_rows}
    for k in sorted(set(opc_map) | set(sdd_map) | set(absz_map)):
        o = opc_map.get(k, {})
        os_ = opc_s_map.get(k, {})
        s = sdd_map.get(k, {})
        a = absz_map.get(k, {})
        comp.append(
            {
                "source": k[0],
                "record": k[1],
                "stratum": o.get("duration_stratum", s.get("duration_stratum", a.get("duration_stratum", ""))),
                "event_hr": o.get("event_hr", a.get("event_hr", "")),
                "opc_L8_alarmed": o.get("alarmed", ""),
                "opc_L8_lead_h": o.get("lead_time_h", ""),
                "opc_L8_basal_alarm_frac": o.get("basal_alarm_frac", ""),
                "opc_L50_alarmed": os_.get("alarmed", ""),
                "opc_L50_lead_h": os_.get("lead_time_h", ""),
                "opc_L50_basal_alarm_frac": os_.get("basal_alarm_frac", ""),
                "sdd_alarmed": s.get("alarmed", ""),
                "sdd_lead_h": s.get("lead_time_h", ""),
                "sdd_mean_TV_approach": s.get("mean_TV_approach", ""),
                "sdd_max_TV_approach": s.get("max_TV_approach", ""),
                "sdd_basal_high_tv_frac_diag": s.get("basal_high_tv_frac_diag", ""),
                "absz_tau_alarmed": a.get("alarmed", ""),
                "absz_tau_lead_h": a.get("lead_time_h", ""),
            }
        )

    write_csv(out_dir / "ordinal_opc_per_record.csv", opc_rows)
    write_csv(out_dir / "ordinal_opc_L50_per_record.csv", opc_scale_rows)
    write_csv(out_dir / "ordinal_sdd_per_record.csv", sdd_rows)
    write_csv(out_dir / "ordinal_absz_per_record.csv", absz_rows)
    write_csv(out_dir / "ordinal_opc_sdd_absz_comparison.csv", comp)
    if nsrdb_opc:
        write_csv(out_dir / "ordinal_opc_nsrdb_controls.csv", nsrdb_opc)
        write_csv(out_dir / "ordinal_sdd_nsrdb_controls.csv", nsrdb_sdd)

    summary = summarize(opc_rows, sdd_rows, absz_rows, opc_scale_rows)
    if nsrdb_opc:
        fracs_o = [float(r["post_basal_alarm_frac"]) for r in nsrdb_opc if "post_basal_alarm_frac" in r and np.isfinite(float(r.get("post_basal_alarm_frac", "nan") or "nan"))]
        fracs_s = [float(r["post_basal_alarm_frac"]) for r in nsrdb_sdd if "post_basal_alarm_frac" in r and np.isfinite(float(r.get("post_basal_alarm_frac", "nan") or "nan"))]
        summary["nsrdb_control_flavor"] = {
            "n": len(nsrdb_opc),
            "opc_mean_post_basal_alarm_frac": float(np.mean(fracs_o)) if fracs_o else float("nan"),
            "sdd_mean_post_basal_alarm_frac": float(np.mean(fracs_s)) if fracs_s else float("nan"),
            "note": "fraction of post-basal symbol endpoints alarmed; not FAR/24h",
        }
    with open(out_dir / "ordinal_exploratory_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    print("\n=== SUMMARY (exploratory) ===")
    print(json.dumps(summary, indent=2))
    print(f"\nWrote results under {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
