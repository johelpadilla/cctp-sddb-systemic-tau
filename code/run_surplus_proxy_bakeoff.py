#!/usr/bin/env python3
"""
Phase α — Proxy bake-off diagnostic (independent of I0 / detectors).

Compares frozen bivariate proxy variants under *identical* Bandt–Pompe m=3,
delay=1, L=50, TV surplus definition. No threshold grids, no alarms, no
detector promotion.

Arms (see code/proxy_variants.py):
  R0, P1 (AR1 residual), P2 (lag-10 |ΔRR|), P3 (signed ΔRR), P4 (|ΔRR| residual on RR)

Metrics:
  - NSRDB floor: mean/median/p90 of S, frac S≥0.40, high-S episode rate, corr(S,MI)
  - Event contrast: NSRDB vs event basals (SDDB+VFDB); approach−basal ΔS; Cohen's d

Outputs under results/:
  surplus_bakeoff_nsrdb_per_record.csv
  surplus_bakeoff_event_per_record.csv
  surplus_bakeoff_proxy_summary.csv
  surplus_bakeoff_summary.json

Reproduce:
  PYTHONPATH=code python3 code/run_surplus_proxy_bakeoff.py
  PYTHONPATH=code python3 code/run_surplus_proxy_bakeoff.py --smoke

Observational only — no clinical claims.
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
from cctp_metrics_core import get_event_and_windows, short_db_windows
from ordinal_detectors.opc_refinements import (
    joint_symbols_from_factors,
    ordinal_synergy_surplus,
)
from proxy_variants import (
    ALL_PROXY_DESCRIPTIONS,
    P2_LAG,
    P2_LAG_IDS,
    PROXY_DESCRIPTIONS,
    PROXY_IDS,
    build_proxy,
)

# Phase-α default set; lag-screen IDs (P2-5/P2-10/P2-20) also accepted via --proxies
KNOWN_PROXY_IDS = tuple(PROXY_IDS) + tuple(P2_LAG_IDS)
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
K_JOINT = BP_ALPHABET * BP_ALPHABET
L_SURPLUS = 50

BASAL_HOURS = 2.0
CONTROL_MAX_HOURS = 12.0

# Same diagnostic episode definition as SURPLUS_NSRDB_DIAGNOSTIC
HIGH_S_ABS = 0.40
MIN_EPISODE_LEN = 5
MERGE_GAP = 10
MI_STRIDE = 25  # sample every 25 surplus samples for corr(S, MI)

# Pre-registered exploratory “promising” bars (not final decision)
FLOOR_MEAN_DROP = 0.25  # ≥25% drop in NSRDB mean S vs R0
FLOOR_FRAC40_DROP = 0.30  # ≥30% drop in frac S≥0.40 vs R0


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


def control_basal(total_h: float, basal_hours: float = BASAL_HOURS) -> Tuple[float, float]:
    basal = (0.25, min(basal_hours, total_h * 0.25))
    if basal[1] <= basal[0]:
        basal = (0.0, max(total_h * 0.2, 0.1))
    return basal


def _safe_stats(x: np.ndarray) -> Dict[str, float]:
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if x.size == 0:
        return {
            "n": 0,
            "mean": float("nan"),
            "median": float("nan"),
            "std": float("nan"),
            "p10": float("nan"),
            "p90": float("nan"),
            "p95": float("nan"),
            "min": float("nan"),
            "max": float("nan"),
            "frac_ge_0p40": float("nan"),
            "frac_ge_0p50": float("nan"),
        }
    qs = np.percentile(x, [10, 50, 90, 95])
    return {
        "n": int(x.size),
        "mean": float(np.mean(x)),
        "median": float(qs[1]),
        "std": float(np.std(x)),
        "p10": float(qs[0]),
        "p90": float(qs[2]),
        "p95": float(qs[3]),
        "min": float(np.min(x)),
        "max": float(np.max(x)),
        "frac_ge_0p40": float(np.mean(x >= 0.40)),
        "frac_ge_0p50": float(np.mean(x >= 0.50)),
    }


def _shannon_entropy_from_probs(p: np.ndarray) -> float:
    p = p[p > 0]
    if p.size == 0:
        return float("nan")
    return float(-np.sum(p * np.log2(p)))


def window_mi_bits(pi1_win: np.ndarray, pi2_win: np.ndarray, k1: int = 6, k2: int = 6) -> float:
    joint = np.zeros((k1, k2), dtype=float)
    for a, b in zip(pi1_win, pi2_win):
        joint[int(a), int(b)] += 1.0
    joint /= float(len(pi1_win))
    p1 = joint.sum(axis=1)
    p2 = joint.sum(axis=0)
    h1 = _shannon_entropy_from_probs(p1)
    h2 = _shannon_entropy_from_probs(p2)
    hj = _shannon_entropy_from_probs(joint.ravel())
    return h1 + h2 - hj


def find_episodes(
    high: np.ndarray,
    t_hr: np.ndarray,
    S: np.ndarray,
    *,
    min_len: int = MIN_EPISODE_LEN,
    merge_gap: int = MERGE_GAP,
) -> List[Dict[str, Any]]:
    high = np.asarray(high, dtype=np.int8)
    T = high.size
    raw: List[Tuple[int, int]] = []
    i = 0
    while i < T:
        if high[i] == 0:
            i += 1
            continue
        j = i
        while j < T and high[j] == 1:
            j += 1
        if j - i >= min_len:
            raw.append((i, j - 1))
        i = j
    if not raw:
        return []
    merged = [raw[0]]
    for a, b in raw[1:]:
        pa, pb = merged[-1]
        if a - pb <= merge_gap:
            merged[-1] = (pa, b)
        else:
            merged.append((a, b))
    eps = []
    for a, b in merged:
        seg = S[a : b + 1]
        peak_i = a + int(np.nanargmax(seg))
        eps.append(
            {
                "start_idx": int(a),
                "end_idx": int(b),
                "peak_idx": int(peak_i),
                "duration_samples": int(b - a + 1),
                "start_hr": float(t_hr[a]) if a < len(t_hr) else float("nan"),
                "end_hr": float(t_hr[b]) if b < len(t_hr) else float("nan"),
                "peak_S": float(S[peak_i]),
                "mean_S": float(np.nanmean(seg)),
            }
        )
    return eps


def cohens_d(a: Sequence[float], b: Sequence[float]) -> float:
    """Cohen's d for a vs b (a_mean - b_mean) / pooled_sd."""
    aa = np.asarray(a, dtype=float)
    bb = np.asarray(b, dtype=float)
    aa = aa[np.isfinite(aa)]
    bb = bb[np.isfinite(bb)]
    if aa.size < 2 or bb.size < 2:
        return float("nan")
    ma, mb = float(np.mean(aa)), float(np.mean(bb))
    sa, sb = float(np.std(aa, ddof=1)), float(np.std(bb, ddof=1))
    # pooled
    n1, n2 = aa.size, bb.size
    sp2 = ((n1 - 1) * sa ** 2 + (n2 - 1) * sb ** 2) / (n1 + n2 - 2)
    if sp2 <= 0:
        return float("nan")
    return (ma - mb) / float(np.sqrt(sp2))


def factors_and_surplus(
    rr: np.ndarray,
    proxy_id: str,
    *,
    L: int = L_SURPLUS,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, int]:
    """Return (pi1, pi2, S, offset)."""
    X = build_proxy(rr, proxy_id)
    S_sym = generate_multivariate_symbols(X, m=M_EMB, delay=DELAY)
    offset = (M_EMB - 1) * DELAY
    if S_sym.size == 0 or S_sym.shape[1] < 2:
        empty = np.array([], dtype=np.int64)
        return empty, empty, np.array([], dtype=float), offset
    pi1 = S_sym[:, 0].astype(np.int64)
    pi2 = S_sym[:, 1].astype(np.int64)
    S = ordinal_synergy_surplus(pi1, pi2, L=L, k1=BP_ALPHABET, k2=BP_ALPHABET)
    return pi1, pi2, S, offset


def corr_s_mi(pi1: np.ndarray, pi2: np.ndarray, S: np.ndarray, L: int) -> float:
    """Pearson corr of S vs ordinal MI on strided windows."""
    T = len(S)
    ss: List[float] = []
    mis: List[float] = []
    for t in range(L - 1, T, MI_STRIDE):
        if not np.isfinite(S[t]):
            continue
        sl = slice(t - L + 1, t + 1)
        mi = window_mi_bits(pi1[sl], pi2[sl])
        if not np.isfinite(mi):
            continue
        ss.append(float(S[t]))
        mis.append(float(mi))
    if len(ss) < 10:
        return float("nan")
    a = np.asarray(ss)
    b = np.asarray(mis)
    if np.std(a) < 1e-12 or np.std(b) < 1e-12:
        return float("nan")
    return float(np.corrcoef(a, b)[0, 1])


def process_nsrdb_record(
    path: Path,
    proxy_id: str,
    *,
    L: int = L_SURPLUS,
    max_hours: float = CONTROL_MAX_HOURS,
    high_abs: float = HIGH_S_ABS,
) -> Dict[str, Any]:
    d = load_npz(path)
    rr = np.asarray(d["rr_ms"], dtype=float)
    t_sec = np.asarray(d["t_sec"], dtype=float)
    t_hr = t_sec / 3600.0
    total_h_full = float(d.get("total_hours", float(np.nanmax(t_hr))))
    # cap like Phase-2 FAR / surplus diagnostic
    if max_hours is not None and total_h_full > max_hours:
        mask = t_hr <= max_hours
        rr = rr[mask]
        t_hr = t_hr[mask]
        total_h = float(max_hours)
    else:
        total_h = float(np.nanmax(t_hr)) if t_hr.size else 0.0

    rec_id = path.stem.replace("nsrdb_", "").replace("_clean", "")
    meta = {
        "proxy_id": proxy_id,
        "cohort": "nsrdb",
        "record_id": rec_id,
        "source_path": path.name,
        "total_hours_used": total_h,
        "L": L,
        "m": M_EMB,
        "high_abs": high_abs,
        "proxy_desc": ALL_PROXY_DESCRIPTIONS.get(
            proxy_id, PROXY_DESCRIPTIONS.get(proxy_id, "")
        ),
    }
    if len(rr) < L + 100 or total_h < BASAL_HOURS + 0.5:
        return {**meta, "skipped": True, "skip_reason": "too_short"}

    pi1, pi2, S, offset = factors_and_surplus(rr, proxy_id, L=L)
    if len(pi1) < L + 10:
        return {**meta, "skipped": True, "skip_reason": "too_few_symbols"}

    t_sym = t_hr[offset : offset + len(pi1)]
    n = min(len(pi1), len(t_sym), len(S))
    pi1, pi2, S, t_sym = pi1[:n], pi2[:n], S[:n], t_sym[:n]
    finite = np.isfinite(S)
    stats_all = _safe_stats(S[finite])

    b0, b1 = control_basal(total_h)
    basal_mask = finite & (t_sym >= b0) & (t_sym < b1)
    search_mask = finite & (t_sym >= b1)
    stats_basal = _safe_stats(S[basal_mask])
    stats_search = _safe_stats(S[search_mask])

    high = (S >= high_abs).astype(np.int8)
    high[~finite] = 0
    high_search = high.copy()
    high_search[t_sym < b1] = 0
    eps = find_episodes(high_search, t_sym, S)
    search_h = float(max(0.0, total_h - b1))
    ep_rate_24h = (len(eps) / search_h * 24.0) if search_h > 0 else float("nan")

    r_s_mi = corr_s_mi(pi1, pi2, S, L)
    # mean MI on strided windows
    mis = []
    for t in range(L - 1, n, MI_STRIDE):
        if not np.isfinite(S[t]):
            continue
        mis.append(window_mi_bits(pi1[t - L + 1 : t + 1], pi2[t - L + 1 : t + 1]))
    mean_mi = float(np.nanmean(mis)) if mis else float("nan")

    return {
        **meta,
        "skipped": False,
        "skip_reason": "",
        "n_symbols": int(n),
        "basal_start": b0,
        "basal_end": b1,
        "search_hours": search_h,
        "n_episodes_S_ge_high_abs": len(eps),
        "episode_rate_per_24h": ep_rate_24h,
        "mean_S": stats_all["mean"],
        "median_S": stats_all["median"],
        "std_S": stats_all["std"],
        "p10_S": stats_all["p10"],
        "p90_S": stats_all["p90"],
        "p95_S": stats_all["p95"],
        "min_S": stats_all["min"],
        "max_S": stats_all["max"],
        "frac_S_ge_0p40": stats_all["frac_ge_0p40"],
        "frac_S_ge_0p50": stats_all["frac_ge_0p50"],
        "n_S": stats_all["n"],
        "mean_S_basal": stats_basal["mean"],
        "mean_S_search": stats_search["mean"],
        "p90_S_basal": stats_basal["p90"],
        "p90_S_search": stats_search["p90"],
        "corr_S_MI": r_s_mi,
        "mean_MI_bits": mean_mi,
    }


def resolve_event_windows(
    source: str, record: str, t_hr: np.ndarray, vfon_hr: float, total_h: float
) -> Tuple[float, Tuple[float, float], Tuple[float, float], str]:
    if source == "sddb":
        event_hr, basal, approach = get_event_and_windows(record, t_hr, vfon_hr)
        return float(event_hr), basal, approach, "holter_analytic"
    event_hr, basal, approach, stratum = short_db_windows(vfon_hr, total_h)
    return float(event_hr), basal, approach, stratum


def process_event_record(
    source: str,
    record: str,
    path: Path,
    proxy_id: str,
    *,
    L: int = L_SURPLUS,
) -> Dict[str, Any]:
    d = load_npz(path)
    rr = np.asarray(d["rr_ms"], dtype=float)
    t_sec = np.asarray(d["t_sec"], dtype=float)
    t_hr = t_sec / 3600.0
    vfon_hr = float(d["vfon_sec"]) / 3600.0
    total_h = float(d.get("total_hours", float(np.nanmax(t_hr))))
    event_label = str(d.get("event_label", "VF/VT" if source == "vfdb" else "SDDB"))

    event_hr, basal, approach, stratum = resolve_event_windows(
        source, record, t_hr, vfon_hr, total_h
    )
    b0, b1 = basal
    a0, a1 = approach

    meta = {
        "proxy_id": proxy_id,
        "cohort": source,
        "record_id": record,
        "source_path": path.name,
        "event_label": event_label,
        "duration_stratum": stratum,
        "event_hr": event_hr,
        "basal_start": b0,
        "basal_end": b1,
        "approach_start": a0,
        "approach_end": a1,
        "total_hours": total_h,
        "L": L,
        "m": M_EMB,
        "proxy_desc": ALL_PROXY_DESCRIPTIONS.get(
            proxy_id, PROXY_DESCRIPTIONS.get(proxy_id, "")
        ),
    }

    pi1, pi2, S, offset = factors_and_surplus(rr, proxy_id, L=L)
    if len(pi1) < L + 10:
        return {**meta, "skipped": True, "skip_reason": "too_few_symbols"}

    t_sym = t_hr[offset : offset + len(pi1)]
    n = min(len(pi1), len(t_sym), len(S))
    pi1, pi2, S, t_sym = pi1[:n], pi2[:n], S[:n], t_sym[:n]
    finite = np.isfinite(S)

    basal_mask = finite & (t_sym >= b0) & (t_sym < b1)
    approach_mask = finite & (t_sym >= a0) & (t_sym < a1)
    pre_mask = finite & (t_sym < event_hr)

    if basal_mask.sum() < max(20, L):
        pre_idx = np.where(pre_mask)[0]
        if pre_idx.size < L:
            return {**meta, "skipped": True, "skip_reason": "no_basal"}
        cut = max(L, pre_idx.size // 3)
        basal_mask = np.zeros_like(finite, dtype=bool)
        basal_mask[pre_idx[:cut]] = True

    stats_basal = _safe_stats(S[basal_mask])
    stats_approach = _safe_stats(S[approach_mask])

    return {
        **meta,
        "skipped": False,
        "skip_reason": "",
        "n_symbols": int(n),
        "n_S_basal": stats_basal["n"],
        "mean_S_basal": stats_basal["mean"],
        "median_S_basal": stats_basal["median"],
        "std_S_basal": stats_basal["std"],
        "p90_S_basal": stats_basal["p90"],
        "frac_S_ge_0p40_basal": stats_basal["frac_ge_0p40"],
        "n_S_approach": stats_approach["n"],
        "mean_S_approach": stats_approach["mean"],
        "median_S_approach": stats_approach["median"],
        "p90_S_approach": stats_approach["p90"],
        "delta_mean_approach_minus_basal": (
            stats_approach["mean"] - stats_basal["mean"]
            if np.isfinite(stats_approach["mean"]) and np.isfinite(stats_basal["mean"])
            else float("nan")
        ),
    }


def list_event_records() -> List[Tuple[str, str, Path]]:
    out: List[Tuple[str, str, Path]] = []
    for rec in ANALYTIC_SDDB + EXTRA_SDDB:
        p = DATA / f"rr_{rec}_clean.npz"
        if p.exists():
            out.append(("sddb", rec, p))
    for p in sorted(RR_EXT.glob("vfdb_*_clean.npz")):
        rec = p.stem.replace("vfdb_", "").replace("_clean", "")
        out.append(("vfdb", rec, p))
    return out


def list_nsrdb() -> List[Path]:
    return sorted(RR_EXT.glob("nsrdb_*_clean.npz"))


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


def aggregate_proxy(
    proxy_id: str,
    nsrdb_rows: List[dict],
    event_rows: List[dict],
) -> Dict[str, Any]:
    nsr = [r for r in nsrdb_rows if r.get("proxy_id") == proxy_id and not r.get("skipped")]
    ev = [r for r in event_rows if r.get("proxy_id") == proxy_id and not r.get("skipped")]
    sddb = [r for r in ev if r.get("cohort") == "sddb"]
    vfdb = [r for r in ev if r.get("cohort") == "vfdb"]

    def mean_of(rows: List[dict], key: str) -> float:
        vals = [float(r[key]) for r in rows if np.isfinite(float(r.get(key, np.nan)))]
        return float(np.mean(vals)) if vals else float("nan")

    def median_of(rows: List[dict], key: str) -> float:
        vals = [float(r[key]) for r in rows if np.isfinite(float(r.get(key, np.nan)))]
        return float(np.median(vals)) if vals else float("nan")

    nsr_means = [float(r["mean_S"]) for r in nsr if np.isfinite(float(r.get("mean_S", np.nan)))]
    nsr_medians = [float(r["median_S"]) for r in nsr if np.isfinite(float(r.get("median_S", np.nan)))]
    nsr_p90 = [float(r["p90_S"]) for r in nsr if np.isfinite(float(r.get("p90_S", np.nan)))]
    nsr_frac40 = [float(r["frac_S_ge_0p40"]) for r in nsr if np.isfinite(float(r.get("frac_S_ge_0p40", np.nan)))]
    nsr_ep = [float(r["episode_rate_per_24h"]) for r in nsr if np.isfinite(float(r.get("episode_rate_per_24h", np.nan)))]
    nsr_corr = [float(r["corr_S_MI"]) for r in nsr if np.isfinite(float(r.get("corr_S_MI", np.nan)))]
    nsr_mi = [float(r["mean_MI_bits"]) for r in nsr if np.isfinite(float(r.get("mean_MI_bits", np.nan)))]

    ev_basal = [float(r["mean_S_basal"]) for r in ev if np.isfinite(float(r.get("mean_S_basal", np.nan)))]
    ev_app = [float(r["mean_S_approach"]) for r in ev if np.isfinite(float(r.get("mean_S_approach", np.nan)))]
    ev_delta = [
        float(r["delta_mean_approach_minus_basal"])
        for r in ev
        if np.isfinite(float(r.get("delta_mean_approach_minus_basal", np.nan)))
    ]
    n_app_gt_basal = sum(1 for d in ev_delta if d > 0)

    d_nsrdb_vs_event = cohens_d(nsr_means, ev_basal)  # nsrdb - event; negative ⇒ event higher
    d_event_vs_nsrdb = -d_nsrdb_vs_event if np.isfinite(d_nsrdb_vs_event) else float("nan")

    return {
        "proxy_id": proxy_id,
        "proxy_desc": ALL_PROXY_DESCRIPTIONS.get(
            proxy_id, PROXY_DESCRIPTIONS.get(proxy_id, "")
        ),
        "n_nsrdb": len(nsr),
        "n_event": len(ev),
        "n_sddb": len(sddb),
        "n_vfdb": len(vfdb),
        # NSRDB floor
        "nsrdb_mean_of_mean_S": float(np.mean(nsr_means)) if nsr_means else float("nan"),
        "nsrdb_mean_of_median_S": float(np.mean(nsr_medians)) if nsr_medians else float("nan"),
        "nsrdb_mean_of_p90_S": float(np.mean(nsr_p90)) if nsr_p90 else float("nan"),
        "nsrdb_mean_frac_S_ge_0p40": float(np.mean(nsr_frac40)) if nsr_frac40 else float("nan"),
        "nsrdb_mean_episode_rate_per_24h": float(np.mean(nsr_ep)) if nsr_ep else float("nan"),
        "nsrdb_median_episode_rate_per_24h": float(np.median(nsr_ep)) if nsr_ep else float("nan"),
        "nsrdb_mean_corr_S_MI": float(np.mean(nsr_corr)) if nsr_corr else float("nan"),
        "nsrdb_mean_MI_bits": float(np.mean(nsr_mi)) if nsr_mi else float("nan"),
        # Event
        "event_mean_of_mean_S_basal": float(np.mean(ev_basal)) if ev_basal else float("nan"),
        "event_mean_of_mean_S_approach": float(np.mean(ev_app)) if ev_app else float("nan"),
        "event_mean_delta_approach_minus_basal": float(np.mean(ev_delta)) if ev_delta else float("nan"),
        "event_median_delta_approach_minus_basal": float(np.median(ev_delta)) if ev_delta else float("nan"),
        "event_frac_approach_gt_basal": (
            float(n_app_gt_basal / len(ev_delta)) if ev_delta else float("nan")
        ),
        "sddb_mean_of_mean_S_basal": mean_of(sddb, "mean_S_basal"),
        "vfdb_mean_of_mean_S_basal": mean_of(vfdb, "mean_S_basal"),
        # Contrasts
        "diff_event_basal_minus_nsrdb_mean": (
            float(np.mean(ev_basal) - np.mean(nsr_means))
            if ev_basal and nsr_means
            else float("nan")
        ),
        "cohens_d_event_basal_vs_nsrdb": d_event_vs_nsrdb,
        "event_basal_gt_nsrdb": (
            bool(np.mean(ev_basal) > np.mean(nsr_means))
            if ev_basal and nsr_means
            else False
        ),
    }


def evaluate_promising(summary_by_proxy: Dict[str, Dict[str, Any]]) -> Dict[str, Any]:
    """Apply pre-registered exploratory promising criteria vs R0."""
    r0 = summary_by_proxy.get("R0")
    if not r0:
        return {"error": "R0 missing"}
    r0_mean = float(r0["nsrdb_mean_of_mean_S"])
    r0_frac = float(r0["nsrdb_mean_frac_S_ge_0p40"])
    out: Dict[str, Any] = {"r0_mean_S": r0_mean, "r0_frac_ge_0p40": r0_frac, "arms": {}}
    for pid, row in summary_by_proxy.items():
        mean_s = float(row["nsrdb_mean_of_mean_S"])
        frac = float(row["nsrdb_mean_frac_S_ge_0p40"])
        mean_drop = (r0_mean - mean_s) / r0_mean if r0_mean > 1e-12 else float("nan")
        frac_drop = (r0_frac - frac) / r0_frac if r0_frac > 1e-12 else float("nan")
        floor_ok = (
            (np.isfinite(mean_drop) and mean_drop >= FLOOR_MEAN_DROP)
            or (np.isfinite(frac_drop) and frac_drop >= FLOOR_FRAC40_DROP)
        )
        contrast_ok = bool(row.get("event_basal_gt_nsrdb")) or (
            # “at least not inferior”: event basal mean ≥ nsrdb mean − epsilon
            np.isfinite(float(row.get("diff_event_basal_minus_nsrdb_mean", np.nan)))
            and float(row["diff_event_basal_minus_nsrdb_mean"]) >= -1e-6
        )
        # User criterion: event basals still higher OR at least not inferior
        contrast_ok = (
            np.isfinite(float(row.get("diff_event_basal_minus_nsrdb_mean", np.nan)))
            and float(row["diff_event_basal_minus_nsrdb_mean"]) >= -1e-6
        )
        promising = bool(floor_ok and contrast_ok) if pid != "R0" else False
        out["arms"][pid] = {
            "mean_S_drop_frac_vs_r0": mean_drop if pid != "R0" else 0.0,
            "frac40_drop_frac_vs_r0": frac_drop if pid != "R0" else 0.0,
            "floor_criterion_pass": bool(floor_ok) if pid != "R0" else None,
            "contrast_preserved": bool(contrast_ok),
            "promising": promising,
        }
    return out


def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser(description="Phase α proxy bake-off diagnostic")
    ap.add_argument("--smoke", action="store_true", help="2 NSRDB + 2 events, all proxies")
    ap.add_argument("--L", type=int, default=L_SURPLUS)
    ap.add_argument("--max-hours", type=float, default=CONTROL_MAX_HOURS)
    ap.add_argument("--high-abs", type=float, default=HIGH_S_ABS)
    ap.add_argument(
        "--proxies",
        type=str,
        default=",".join(PROXY_IDS),
        help="Comma-separated proxy ids (default: all Phase-α)",
    )
    args = ap.parse_args(argv)

    proxy_ids = [p.strip().upper() for p in args.proxies.split(",") if p.strip()]
    for p in proxy_ids:
        if p not in KNOWN_PROXY_IDS:
            print(f"Unknown proxy {p}; known: {KNOWN_PROXY_IDS}", file=sys.stderr)
            return 2

    nsrdb_paths = list_nsrdb()
    events = list_event_records()
    if args.smoke:
        nsrdb_paths = nsrdb_paths[:2]
        events = events[:2]

    print(f"[bakeoff] proxies={proxy_ids}")
    print(f"[bakeoff] NSRDB n={len(nsrdb_paths)} events n={len(events)}")
    print(f"[bakeoff] L={args.L} lag_P2={P2_LAG} high_abs={args.high_abs} max_h={args.max_hours}")

    nsrdb_rows: List[dict] = []
    event_rows: List[dict] = []

    for pid in proxy_ids:
        print(
            f"[bakeoff] --- proxy {pid}: "
            f"{ALL_PROXY_DESCRIPTIONS.get(pid, PROXY_DESCRIPTIONS.get(pid, ''))}"
        )
        for path in nsrdb_paths:
            row = process_nsrdb_record(
                path, pid, L=args.L, max_hours=args.max_hours, high_abs=args.high_abs
            )
            nsrdb_rows.append(row)
            if not row.get("skipped"):
                print(
                    f"  nsrdb/{row['record_id']}: meanS={row['mean_S']:.4f} "
                    f"p90={row['p90_S']:.4f} f40={row['frac_S_ge_0p40']:.3f} "
                    f"eps24={row['episode_rate_per_24h']:.1f} r(S,MI)={row['corr_S_MI']:.3f}"
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

    # Aggregate
    summaries = [aggregate_proxy(pid, nsrdb_rows, event_rows) for pid in proxy_ids]
    by_proxy = {s["proxy_id"]: s for s in summaries}
    promising = evaluate_promising(by_proxy)

    # Relative columns vs R0
    r0 = by_proxy.get("R0")
    for s in summaries:
        if r0 and s["proxy_id"] != "R0":
            r0m = float(r0["nsrdb_mean_of_mean_S"])
            r0f = float(r0["nsrdb_mean_frac_S_ge_0p40"])
            m = float(s["nsrdb_mean_of_mean_S"])
            f = float(s["nsrdb_mean_frac_S_ge_0p40"])
            s["rel_mean_S_vs_r0"] = m / r0m if r0m > 1e-12 else float("nan")
            s["mean_S_drop_frac_vs_r0"] = (r0m - m) / r0m if r0m > 1e-12 else float("nan")
            s["frac40_drop_frac_vs_r0"] = (r0f - f) / r0f if r0f > 1e-12 else float("nan")
        else:
            s["rel_mean_S_vs_r0"] = 1.0
            s["mean_S_drop_frac_vs_r0"] = 0.0
            s["frac40_drop_frac_vs_r0"] = 0.0
        arm = promising.get("arms", {}).get(s["proxy_id"], {})
        s["promising"] = arm.get("promising", False)
        s["floor_criterion_pass"] = arm.get("floor_criterion_pass")
        s["contrast_preserved"] = arm.get("contrast_preserved")

    p_nsr = RES / "surplus_bakeoff_nsrdb_per_record.csv"
    p_ev = RES / "surplus_bakeoff_event_per_record.csv"
    p_sum = RES / "surplus_bakeoff_proxy_summary.csv"
    p_json = RES / "surplus_bakeoff_summary.json"

    write_csv(p_nsr, nsrdb_rows)
    write_csv(p_ev, event_rows)
    write_csv(p_sum, summaries)

    payload = {
        "exercise": "surplus_proxy_bakeoff_phase_alpha",
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
            "p2_lag": P2_LAG,
            "proxies": {
                pid: ALL_PROXY_DESCRIPTIONS.get(pid, PROXY_DESCRIPTIONS.get(pid, ""))
                for pid in proxy_ids
            },
            "promising_criteria": {
                "mean_S_drop_ge": FLOOR_MEAN_DROP,
                "frac40_drop_ge": FLOOR_FRAC40_DROP,
                "contrast": "event_basal_mean_S >= nsrdb_mean_S",
            },
        },
        "proxy_summaries": summaries,
        "promising_evaluation": promising,
        "artifacts": {
            "nsrdb_per_record": str(p_nsr.name),
            "event_per_record": str(p_ev.name),
            "proxy_summary": str(p_sum.name),
        },
        "disclaimer": "Diagnostic only. No clinical claims. Not a detector.",
    }
    p_json.write_text(json.dumps(payload, indent=2, default=str))

    print("\n[bakeoff] === PROXY SUMMARY ===")
    hdr = (
        f"{'proxy':4s} {'meanS':>7s} {'p90':>7s} {'f≥0.40':>7s} {'eps/24h':>8s} "
        f"{'r(S,MI)':>7s} {'evBas':>7s} {'Δapp':>7s} {'d_ev':>6s} {'prom':>5s}"
    )
    print(hdr)
    for s in summaries:
        print(
            f"{s['proxy_id']:4s} "
            f"{s['nsrdb_mean_of_mean_S']:7.4f} "
            f"{s['nsrdb_mean_of_p90_S']:7.4f} "
            f"{s['nsrdb_mean_frac_S_ge_0p40']:7.3f} "
            f"{s['nsrdb_mean_episode_rate_per_24h']:8.1f} "
            f"{s['nsrdb_mean_corr_S_MI']:7.3f} "
            f"{s['event_mean_of_mean_S_basal']:7.4f} "
            f"{s['event_mean_delta_approach_minus_basal']:7.4f} "
            f"{s['cohens_d_event_basal_vs_nsrdb']:6.2f} "
            f"{'Y' if s.get('promising') else 'n':>5s}"
        )
    print(f"[bakeoff] wrote {p_nsr.name}, {p_ev.name}, {p_sum.name}, {p_json.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
