#!/usr/bin/env python3
"""
Diagnostic analysis of synergistic surplus S_t on NSRDB controls
(and comparison to basal windows of VFDB/SDDB event records).

Independent exercise — NOT a detector, NOT parameter tuning, NOT I0 continuation.

S_t = TV(P_joint, P1 ⊗ P2) on Bandt–Pompe bivariate factors (m=3, K=36, L=50),
same encoding as the closed surplus-detector track.

Outputs under results/:
  surplus_diag_nsrdb_per_record.csv
  surplus_diag_nsrdb_episodes.csv
  surplus_diag_nsrdb_high_vs_low.csv
  surplus_diag_event_basal_per_record.csv
  surplus_diag_cohort_comparison.csv
  surplus_diag_summary.json

Observational only — no clinical claims.
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
from cctp_metrics_core import (
    build_bivariate_proxy,
    get_event_and_windows,
    short_db_windows,
)
from ordinal_detectors.opc_refinements import (
    joint_symbols_from_factors,
    ordinal_synergy_surplus,
)
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
L_SURPLUS = 50

# Phase-2–aligned control windowing (same family as FAR protocol)
BASAL_HOURS = 2.0
CONTROL_MAX_HOURS = 12.0

# Episode definition (diagnostic, not alarm design)
# Absolute TV threshold for "high surplus" peaks; also report relative markers.
HIGH_S_ABS = 0.40  # near historical θ_S absolute band; TV ∈ [0,1]
MIN_EPISODE_LEN = 5  # consecutive samples above threshold
MERGE_GAP = 10  # merge peaks separated by ≤ this many samples


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


def bivariate_factors(rr: np.ndarray) -> Tuple[np.ndarray, np.ndarray, int]:
    X = build_bivariate_proxy(np.asarray(rr, dtype=float))
    S = generate_multivariate_symbols(X, m=M_EMB, delay=DELAY)
    offset = (M_EMB - 1) * DELAY
    if S.size == 0 or S.shape[1] < 2:
        return np.array([], dtype=np.int64), np.array([], dtype=np.int64), offset
    return S[:, 0].astype(np.int64), S[:, 1].astype(np.int64), offset


def control_basal(total_h: float, basal_hours: float = BASAL_HOURS) -> Tuple[float, float]:
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


def _safe_stats(x: np.ndarray) -> Dict[str, float]:
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if x.size == 0:
        return {
            "n": 0,
            "mean": float("nan"),
            "median": float("nan"),
            "std": float("nan"),
            "min": float("nan"),
            "max": float("nan"),
            "p10": float("nan"),
            "p25": float("nan"),
            "p75": float("nan"),
            "p90": float("nan"),
            "p95": float("nan"),
            "p99": float("nan"),
            "iqr": float("nan"),
            "cv": float("nan"),
            "frac_ge_0p15": float("nan"),
            "frac_ge_0p25": float("nan"),
            "frac_ge_0p40": float("nan"),
            "frac_ge_0p50": float("nan"),
        }
    qs = np.percentile(x, [10, 25, 50, 75, 90, 95, 99])
    mean = float(np.mean(x))
    std = float(np.std(x))
    return {
        "n": int(x.size),
        "mean": mean,
        "median": float(qs[2]),
        "std": std,
        "min": float(np.min(x)),
        "max": float(np.max(x)),
        "p10": float(qs[0]),
        "p25": float(qs[1]),
        "p75": float(qs[3]),
        "p90": float(qs[4]),
        "p95": float(qs[5]),
        "p99": float(qs[6]),
        "iqr": float(qs[3] - qs[1]),
        "cv": float(std / mean) if abs(mean) > 1e-12 else float("nan"),
        "frac_ge_0p15": float(np.mean(x >= 0.15)),
        "frac_ge_0p25": float(np.mean(x >= 0.25)),
        "frac_ge_0p40": float(np.mean(x >= 0.40)),
        "frac_ge_0p50": float(np.mean(x >= 0.50)),
    }


def _shannon_entropy(counts: np.ndarray) -> float:
    total = float(counts.sum())
    if total <= 0:
        return float("nan")
    p = counts[counts > 0] / total
    return float(-np.sum(p * np.log2(p)))


def window_covariates(
    rr_win: np.ndarray,
    pi1_win: np.ndarray,
    pi2_win: np.ndarray,
    sigma_win: np.ndarray,
    *,
    k1: int = BP_ALPHABET,
    k2: int = BP_ALPHABET,
    K: int = K_JOINT,
) -> Dict[str, float]:
    """Describe one L-window of RR + ordinal factors."""
    rr = np.asarray(rr_win, dtype=float)
    out: Dict[str, float] = {}
    out["rr_mean"] = float(np.mean(rr))
    out["rr_std"] = float(np.std(rr))
    out["rr_cv"] = float(out["rr_std"] / out["rr_mean"]) if out["rr_mean"] > 1e-9 else float("nan")
    out["rr_min"] = float(np.min(rr))
    out["rr_max"] = float(np.max(rr))
    out["rr_range"] = out["rr_max"] - out["rr_min"]
    if rr.size >= 3:
        drr = np.diff(rr)
        out["rmssd"] = float(np.sqrt(np.mean(drr ** 2)))
        if np.std(rr[:-1]) > 1e-9 and np.std(rr[1:]) > 1e-9:
            out["rr_ar1"] = float(np.corrcoef(rr[:-1], rr[1:])[0, 1])
        else:
            out["rr_ar1"] = float("nan")
    else:
        out["rmssd"] = float("nan")
        out["rr_ar1"] = float("nan")

    # Symbol support / diversity
    n_unique_joint = int(len(np.unique(sigma_win)))
    out["n_unique_joint"] = float(n_unique_joint)
    out["diversity"] = float(n_unique_joint / float(K))
    out["n_unique_pi1"] = float(len(np.unique(pi1_win)))
    out["n_unique_pi2"] = float(len(np.unique(pi2_win)))

    # Max relative frequency (repetition / locking proxy)
    _, counts_j = np.unique(sigma_win, return_counts=True)
    out["max_joint_frac"] = float(np.max(counts_j) / float(len(sigma_win)))
    out["top3_joint_frac"] = float(np.sum(np.sort(counts_j)[-3:]) / float(len(sigma_win)))

    # Histograms → entropy + TV components already in S; also mutual-info proxy
    joint = np.zeros((k1, k2), dtype=float)
    for a, b in zip(pi1_win, pi2_win):
        joint[int(a), int(b)] += 1.0
    joint /= float(len(pi1_win))
    p1 = joint.sum(axis=1)
    p2 = joint.sum(axis=0)
    out["H_joint"] = _shannon_entropy(joint.ravel() * len(pi1_win))
    out["H_pi1"] = _shannon_entropy(p1 * len(pi1_win))
    out["H_pi2"] = _shannon_entropy(p2 * len(pi1_win))
    # I(X;Y) = H(X)+H(Y)-H(X,Y) in bits (same as TV surplus is dependence measure)
    out["MI_bits"] = out["H_pi1"] + out["H_pi2"] - out["H_joint"]
    # Fraction of joint cells used
    out["joint_support_frac"] = float(np.sum(joint > 0) / float(k1 * k2))

    # Dominating factor: is surplus driven by concentrated margins or pure coupling?
    # Normalized MI-ish: MI / min(H1,H2)
    hmin = min(out["H_pi1"], out["H_pi2"])
    out["NMI"] = float(out["MI_bits"] / hmin) if hmin > 1e-9 else float("nan")

    return out


def find_episodes(
    high: np.ndarray,
    t_hr: np.ndarray,
    S: np.ndarray,
    *,
    min_len: int = MIN_EPISODE_LEN,
    merge_gap: int = MERGE_GAP,
) -> List[Dict[str, Any]]:
    """Merge runs of high==1 into episodes; return start/end indices and peak S."""
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
    # merge nearby
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
                "peak_hr": float(t_hr[peak_i]) if peak_i < len(t_hr) else float("nan"),
                "duration_h": float(t_hr[b] - t_hr[a]) if b < len(t_hr) else float("nan"),
                "peak_S": float(S[peak_i]),
                "mean_S": float(np.nanmean(seg)),
            }
        )
    return eps


def rolling_window_features(
    rr_sym: np.ndarray,
    pi1: np.ndarray,
    pi2: np.ndarray,
    sigma: np.ndarray,
    S: np.ndarray,
    *,
    L: int,
    stride: int,
) -> Tuple[np.ndarray, Dict[str, np.ndarray]]:
    """
    Sample windows at `stride` for association analysis.
    Returns (indices, feature dict of arrays).
    """
    T = len(S)
    idxs = list(range(L - 1, T, stride))
    keys = [
        "S",
        "rr_mean",
        "rr_std",
        "rr_cv",
        "rr_range",
        "rmssd",
        "rr_ar1",
        "n_unique_joint",
        "diversity",
        "n_unique_pi1",
        "n_unique_pi2",
        "max_joint_frac",
        "top3_joint_frac",
        "H_joint",
        "H_pi1",
        "H_pi2",
        "MI_bits",
        "NMI",
        "joint_support_frac",
    ]
    feats: Dict[str, List[float]] = {k: [] for k in keys}
    keep_idx: List[int] = []
    for t in idxs:
        if not np.isfinite(S[t]):
            continue
        sl = slice(t - L + 1, t + 1)
        cov = window_covariates(rr_sym[sl], pi1[sl], pi2[sl], sigma[sl])
        keep_idx.append(t)
        feats["S"].append(float(S[t]))
        for k in keys:
            if k == "S":
                continue
            feats[k].append(float(cov[k]))
    arr = {k: np.asarray(v, dtype=float) for k, v in feats.items()}
    return np.asarray(keep_idx, dtype=np.int64), arr


def high_vs_low_summary(
    feats: Dict[str, np.ndarray],
    *,
    high_frac: float = 0.10,
) -> List[Dict[str, Any]]:
    """Compare top high_frac vs bottom high_frac of S on each covariate."""
    S = feats["S"]
    if S.size < 50:
        return []
    thr_hi = float(np.percentile(S, 100.0 * (1.0 - high_frac)))
    thr_lo = float(np.percentile(S, 100.0 * high_frac))
    hi = S >= thr_hi
    lo = S <= thr_lo
    rows = []
    for k, v in feats.items():
        if k == "S":
            continue
        vh = v[hi]
        vl = v[lo]
        vh = vh[np.isfinite(vh)]
        vl = vl[np.isfinite(vl)]
        mh = float(np.mean(vh)) if vh.size else float("nan")
        ml = float(np.mean(vl)) if vl.size else float("nan")
        # Cohen's d (pooled)
        if vh.size > 2 and vl.size > 2:
            sp = np.sqrt(
                ((vh.size - 1) * np.var(vh, ddof=1) + (vl.size - 1) * np.var(vl, ddof=1))
                / max(vh.size + vl.size - 2, 1)
            )
            d = float((mh - ml) / sp) if sp > 1e-12 else float("nan")
        else:
            d = float("nan")
        # Spearman-ish rank corr of S with covariate (sample)
        if v.size > 20 and np.isfinite(v).sum() > 20:
            mask = np.isfinite(S) & np.isfinite(v)
            if mask.sum() > 20:
                rs = float(np.corrcoef(S[mask], v[mask])[0, 1])  # Pearson on values; ok for ranking
            else:
                rs = float("nan")
        else:
            rs = float("nan")
        rows.append(
            {
                "feature": k,
                "mean_high_S": mh,
                "mean_low_S": ml,
                "delta_high_minus_low": mh - ml if np.isfinite(mh) and np.isfinite(ml) else float("nan"),
                "cohens_d": d,
                "pearson_r_with_S": rs,
                "n_high": int(vh.size),
                "n_low": int(vl.size),
                "thr_high_S": thr_hi,
                "thr_low_S": thr_lo,
            }
        )
    rows.sort(key=lambda r: -abs(r.get("cohens_d") or 0.0) if np.isfinite(r.get("cohens_d", np.nan)) else 0)
    return rows


def process_nsrdb_record(
    path: Path,
    *,
    max_hours: float = CONTROL_MAX_HOURS,
    basal_hours: float = BASAL_HOURS,
    L: int = L_SURPLUS,
    high_abs: float = HIGH_S_ABS,
    cov_stride: int = 25,
) -> Tuple[Dict[str, Any], List[Dict[str, Any]], List[Dict[str, Any]], Dict[str, np.ndarray]]:
    d = load_npz(path)
    rr = np.asarray(d["rr_ms"], dtype=float)
    t_sec = np.asarray(d["t_sec"], dtype=float)
    t_hr = t_sec / 3600.0
    rec = str(d.get("record_id", path.stem))
    total_h = float(t_hr[-1]) if len(t_hr) else 0.0
    interp_frac = float(d.get("interp_frac", np.nan)) if "interp_frac" in d else float("nan")
    n_invalid = float(d.get("n_invalid", np.nan)) if "n_invalid" in d else float("nan")

    if max_hours is not None and total_h > max_hours:
        keep = t_hr <= max_hours
        rr = rr[keep]
        t_hr = t_hr[keep]
        total_h = float(t_hr[-1]) if len(t_hr) else 0.0

    meta_base = {
        "cohort": "nsrdb",
        "record_id": rec,
        "source_path": path.name,
        "total_hours_used": total_h,
        "interp_frac": interp_frac,
        "n_invalid": n_invalid,
        "L": L,
        "m": M_EMB,
        "K": K_JOINT,
    }

    if len(rr) < L + 100 or total_h < basal_hours + 0.5:
        return (
            {**meta_base, "skipped": True, "skip_reason": "too_short"},
            [],
            [],
            {},
        )

    pi1, pi2, offset = bivariate_factors(rr)
    if len(pi1) < L + 10:
        return (
            {**meta_base, "skipped": True, "skip_reason": "too_few_symbols"},
            [],
            [],
            {},
        )

    t_sym = t_hr[offset : offset + len(pi1)]
    n = min(len(pi1), len(t_sym))
    pi1, pi2, t_sym = pi1[:n], pi2[:n], t_sym[:n]
    rr_sym = rr[offset : offset + n]
    sigma = joint_symbols_from_factors(pi1, pi2, k2=BP_ALPHABET)

    S = ordinal_synergy_surplus(pi1, pi2, L=L, k1=BP_ALPHABET, k2=BP_ALPHABET)
    finite = np.isfinite(S)
    stats_all = _safe_stats(S[finite])

    basal = control_basal(total_h, basal_hours)
    b0, b1 = basal
    basal_mask = finite & (t_sym >= b0) & (t_sym < b1)
    search_mask = finite & (t_sym >= b1)
    stats_basal = _safe_stats(S[basal_mask])
    stats_search = _safe_stats(S[search_mask])

    # Episodes of absolute high S (diagnostic)
    high = (S >= high_abs).astype(np.int8)
    high[~finite] = 0
    # only count episodes in search (post-basal) for rate
    high_search = high.copy()
    high_search[t_sym < b1] = 0
    eps = find_episodes(high_search, t_sym, S)
    search_h = float(max(0.0, total_h - b1))
    ep_rate_24h = (len(eps) / search_h * 24.0) if search_h > 0 else float("nan")

    # Covariate sampling
    idxs, feats = rolling_window_features(
        rr_sym, pi1, pi2, sigma, S, L=L, stride=cov_stride
    )
    hvl = high_vs_low_summary(feats) if feats else []

    # Attach record id to episode / hvl rows
    ep_rows = []
    for e in eps:
        ep_rows.append({**meta_base, **e, "high_abs": high_abs, "region": "search"})
        # peak covariates
        pk = e["peak_idx"]
        if pk >= L - 1:
            sl = slice(pk - L + 1, pk + 1)
            cov = window_covariates(rr_sym[sl], pi1[sl], pi2[sl], sigma[sl])
            for k, v in cov.items():
                ep_rows[-1][f"peak_{k}"] = v

    hvl_rows = [{**meta_base, **r} for r in hvl]

    # Record-level summary row
    row = {
        **meta_base,
        "skipped": False,
        "skip_reason": "",
        "n_symbols": int(n),
        "basal_start": b0,
        "basal_end": b1,
        "search_hours": search_h,
        "n_episodes_S_ge_high_abs": len(eps),
        "episode_rate_per_24h": ep_rate_24h,
        "high_abs": high_abs,
        "mean_S_all": stats_all["mean"],
        "median_S_all": stats_all["median"],
        "std_S_all": stats_all["std"],
        "p10_S_all": stats_all["p10"],
        "p25_S_all": stats_all["p25"],
        "p75_S_all": stats_all["p75"],
        "p90_S_all": stats_all["p90"],
        "p95_S_all": stats_all["p95"],
        "p99_S_all": stats_all["p99"],
        "min_S_all": stats_all["min"],
        "max_S_all": stats_all["max"],
        "iqr_S_all": stats_all["iqr"],
        "frac_S_ge_0p15": stats_all["frac_ge_0p15"],
        "frac_S_ge_0p25": stats_all["frac_ge_0p25"],
        "frac_S_ge_0p40": stats_all["frac_ge_0p40"],
        "frac_S_ge_0p50": stats_all["frac_ge_0p50"],
        "mean_S_basal": stats_basal["mean"],
        "median_S_basal": stats_basal["median"],
        "std_S_basal": stats_basal["std"],
        "p90_S_basal": stats_basal["p90"],
        "p95_S_basal": stats_basal["p95"],
        "n_S_basal": stats_basal["n"],
        "mean_S_search": stats_search["mean"],
        "median_S_search": stats_search["median"],
        "std_S_search": stats_search["std"],
        "p90_S_search": stats_search["p90"],
        "p95_S_search": stats_search["p95"],
        "n_S_search": stats_search["n"],
        "delta_mean_search_minus_basal": (
            stats_search["mean"] - stats_basal["mean"]
            if np.isfinite(stats_search["mean"]) and np.isfinite(stats_basal["mean"])
            else float("nan")
        ),
        # mean of key covariates over all sampled windows
        "mean_diversity": float(np.nanmean(feats["diversity"])) if feats else float("nan"),
        "mean_MI_bits": float(np.nanmean(feats["MI_bits"])) if feats else float("nan"),
        "mean_rr_cv": float(np.nanmean(feats["rr_cv"])) if feats else float("nan"),
        "mean_rmssd": float(np.nanmean(feats["rmssd"])) if feats else float("nan"),
        "mean_max_joint_frac": float(np.nanmean(feats["max_joint_frac"])) if feats else float("nan"),
        "mean_NMI": float(np.nanmean(feats["NMI"])) if feats else float("nan"),
    }
    return row, ep_rows, hvl_rows, feats


def resolve_event_windows(
    source: str, record: str, t_hr: np.ndarray, vfon_hr: float, total_h: float
) -> Tuple[float, Tuple[float, float], Tuple[float, float], str]:
    if source == "sddb":
        event_hr, basal, approach = get_event_and_windows(record, t_hr, vfon_hr)
        return float(event_hr), basal, approach, "holter_analytic"
    event_hr, basal, approach, stratum = short_db_windows(vfon_hr, total_h)
    return float(event_hr), basal, approach, stratum


def process_event_basal(
    source: str,
    record: str,
    path: Path,
    *,
    L: int = L_SURPLUS,
    high_abs: float = HIGH_S_ABS,
    cov_stride: int = 10,
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
        "K": K_JOINT,
    }

    pi1, pi2, offset = bivariate_factors(rr)
    if len(pi1) < L + 10:
        return {**meta, "skipped": True, "skip_reason": "too_few_symbols"}

    t_sym = t_hr[offset : offset + len(pi1)]
    n = min(len(pi1), len(t_sym))
    pi1, pi2, t_sym = pi1[:n], pi2[:n], t_sym[:n]
    rr_sym = rr[offset : offset + n]
    sigma = joint_symbols_from_factors(pi1, pi2, k2=BP_ALPHABET)
    S = ordinal_synergy_surplus(pi1, pi2, L=L, k1=BP_ALPHABET, k2=BP_ALPHABET)
    finite = np.isfinite(S)

    basal_mask = finite & (t_sym >= b0) & (t_sym < b1)
    approach_mask = finite & (t_sym >= a0) & (t_sym < a1)
    # pre-event all finite before event
    pre_mask = finite & (t_sym < event_hr)

    if basal_mask.sum() < max(20, L):
        # fallback: first third of pre-event
        pre_idx = np.where(pre_mask)[0]
        if pre_idx.size < L:
            return {**meta, "skipped": True, "skip_reason": "no_basal"}
        cut = max(L, pre_idx.size // 3)
        basal_mask = np.zeros_like(finite)
        basal_mask[pre_idx[:cut]] = True

    stats_basal = _safe_stats(S[basal_mask])
    stats_approach = _safe_stats(S[approach_mask])
    stats_pre = _safe_stats(S[pre_mask])

    # light covariates on basal only
    idxs, feats = rolling_window_features(
        rr_sym, pi1, pi2, sigma, S, L=L, stride=cov_stride
    )
    if feats and idxs.size:
        # restrict to basal times
        basal_feat_mask = (t_sym[idxs] >= b0) & (t_sym[idxs] < b1)
        if basal_feat_mask.sum() < 5:
            basal_feat_mask = np.ones(len(idxs), dtype=bool)
        mean_div = float(np.nanmean(feats["diversity"][basal_feat_mask]))
        mean_mi = float(np.nanmean(feats["MI_bits"][basal_feat_mask]))
        mean_cv = float(np.nanmean(feats["rr_cv"][basal_feat_mask]))
        mean_rmssd = float(np.nanmean(feats["rmssd"][basal_feat_mask]))
        mean_maxf = float(np.nanmean(feats["max_joint_frac"][basal_feat_mask]))
    else:
        mean_div = mean_mi = mean_cv = mean_rmssd = mean_maxf = float("nan")

    high = (S >= high_abs) & basal_mask
    n_high_basal = int(high.sum())
    frac_high_basal = float(n_high_basal / basal_mask.sum()) if basal_mask.sum() else float("nan")

    return {
        **meta,
        "skipped": False,
        "skip_reason": "",
        "n_symbols": int(n),
        "n_S_basal": stats_basal["n"],
        "mean_S_basal": stats_basal["mean"],
        "median_S_basal": stats_basal["median"],
        "std_S_basal": stats_basal["std"],
        "p10_S_basal": stats_basal["p10"],
        "p90_S_basal": stats_basal["p90"],
        "p95_S_basal": stats_basal["p95"],
        "min_S_basal": stats_basal["min"],
        "max_S_basal": stats_basal["max"],
        "frac_S_ge_0p40_basal": stats_basal["frac_ge_0p40"],
        "frac_S_ge_0p50_basal": stats_basal["frac_ge_0p50"],
        "n_high_S_basal_samples": n_high_basal,
        "frac_high_S_basal": frac_high_basal,
        "mean_S_approach": stats_approach["mean"],
        "median_S_approach": stats_approach["median"],
        "std_S_approach": stats_approach["std"],
        "p90_S_approach": stats_approach["p90"],
        "delta_mean_approach_minus_basal": (
            stats_approach["mean"] - stats_basal["mean"]
            if np.isfinite(stats_approach["mean"]) and np.isfinite(stats_basal["mean"])
            else float("nan")
        ),
        "mean_S_pre_event": stats_pre["mean"],
        "mean_diversity_basal": mean_div,
        "mean_MI_bits_basal": mean_mi,
        "mean_rr_cv_basal": mean_cv,
        "mean_rmssd_basal": mean_rmssd,
        "mean_max_joint_frac_basal": mean_maxf,
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


def aggregate_hvl(hvl_rows: List[dict]) -> List[dict]:
    """Pool high-vs-low feature rows across NSRDB records (mean of per-record stats)."""
    by_feat: Dict[str, List[dict]] = {}
    for r in hvl_rows:
        by_feat.setdefault(str(r["feature"]), []).append(r)
    out = []
    for feat, rs in sorted(by_feat.items()):
        ds = [float(r["cohens_d"]) for r in rs if np.isfinite(float(r.get("cohens_d", np.nan)))]
        pr = [float(r["pearson_r_with_S"]) for r in rs if np.isfinite(float(r.get("pearson_r_with_S", np.nan)))]
        dh = [float(r["mean_high_S"]) for r in rs if np.isfinite(float(r.get("mean_high_S", np.nan)))]
        dl = [float(r["mean_low_S"]) for r in rs if np.isfinite(float(r.get("mean_low_S", np.nan)))]
        out.append(
            {
                "feature": feat,
                "n_records": len(rs),
                "mean_cohens_d": float(np.mean(ds)) if ds else float("nan"),
                "median_cohens_d": float(np.median(ds)) if ds else float("nan"),
                "mean_pearson_r": float(np.mean(pr)) if pr else float("nan"),
                "mean_of_mean_high": float(np.mean(dh)) if dh else float("nan"),
                "mean_of_mean_low": float(np.mean(dl)) if dl else float("nan"),
                "mean_delta": float(np.mean(dh) - np.mean(dl)) if dh and dl else float("nan"),
            }
        )
    out.sort(
        key=lambda r: -abs(r["mean_cohens_d"])
        if np.isfinite(r.get("mean_cohens_d", np.nan))
        else 0
    )
    return out


def cohort_comparison(
    nsrdb_rows: List[dict],
    event_rows: List[dict],
) -> List[dict]:
    """Compare NSRDB (all/search/basal) vs event basals (and approaches)."""
    def pack(name: str, vals: List[float], extra: Optional[dict] = None) -> dict:
        st = _safe_stats(np.asarray(vals, dtype=float))
        row = {"group": name, **{f"S_{k}": v for k, v in st.items()}}
        if extra:
            row.update(extra)
        return row

    rows = []
    # NSRDB per-record means of S (record-level distribution of means)
    nsr = [r for r in nsrdb_rows if not r.get("skipped")]
    rows.append(
        pack(
            "nsrdb_record_mean_S_all",
            [float(r["mean_S_all"]) for r in nsr],
            {"unit": "per_record_mean", "n_records": len(nsr)},
        )
    )
    rows.append(
        pack(
            "nsrdb_record_mean_S_basal",
            [float(r["mean_S_basal"]) for r in nsr],
            {"unit": "per_record_mean", "n_records": len(nsr)},
        )
    )
    rows.append(
        pack(
            "nsrdb_record_mean_S_search",
            [float(r["mean_S_search"]) for r in nsr],
            {"unit": "per_record_mean", "n_records": len(nsr)},
        )
    )
    rows.append(
        pack(
            "nsrdb_record_p90_S_all",
            [float(r["p90_S_all"]) for r in nsr],
            {"unit": "per_record_p90", "n_records": len(nsr)},
        )
    )
    rows.append(
        pack(
            "nsrdb_record_frac_S_ge_0p40",
            [float(r["frac_S_ge_0p40"]) for r in nsr],
            {"unit": "per_record_frac", "n_records": len(nsr)},
        )
    )
    rows.append(
        pack(
            "nsrdb_episode_rate_per_24h",
            [float(r["episode_rate_per_24h"]) for r in nsr],
            {"unit": "episodes_per_24h", "n_records": len(nsr)},
        )
    )

    ev = [r for r in event_rows if not r.get("skipped")]
    rows.append(
        pack(
            "event_basal_mean_S",
            [float(r["mean_S_basal"]) for r in ev],
            {"unit": "per_record_mean", "n_records": len(ev)},
        )
    )
    rows.append(
        pack(
            "event_basal_p90_S",
            [float(r["p90_S_basal"]) for r in ev],
            {"unit": "per_record_p90", "n_records": len(ev)},
        )
    )
    rows.append(
        pack(
            "event_approach_mean_S",
            [float(r["mean_S_approach"]) for r in ev if np.isfinite(float(r.get("mean_S_approach", np.nan)))],
            {
                "unit": "per_record_mean",
                "n_records": sum(
                    1
                    for r in ev
                    if np.isfinite(float(r.get("mean_S_approach", np.nan)))
                ),
            },
        )
    )
    for source in ("sddb", "vfdb"):
        sub = [r for r in ev if r.get("cohort") == source]
        if not sub:
            continue
        rows.append(
            pack(
                f"{source}_basal_mean_S",
                [float(r["mean_S_basal"]) for r in sub],
                {"unit": "per_record_mean", "n_records": len(sub)},
            )
        )
        rows.append(
            pack(
                f"{source}_approach_mean_S",
                [
                    float(r["mean_S_approach"])
                    for r in sub
                    if np.isfinite(float(r.get("mean_S_approach", np.nan)))
                ],
                {"unit": "per_record_mean", "n_records": len(sub)},
            )
        )

    # Covariate levels: NSRDB mean diversity/MI vs event basal
    rows.append(
        {
            "group": "nsrdb_mean_diversity",
            "S_mean": float(np.nanmean([float(r["mean_diversity"]) for r in nsr])),
            "S_median": float(np.nanmedian([float(r["mean_diversity"]) for r in nsr])),
            "n_records": len(nsr),
            "unit": "diversity",
        }
    )
    rows.append(
        {
            "group": "event_basal_mean_diversity",
            "S_mean": float(np.nanmean([float(r["mean_diversity_basal"]) for r in ev])),
            "S_median": float(np.nanmedian([float(r["mean_diversity_basal"]) for r in ev])),
            "n_records": len(ev),
            "unit": "diversity",
        }
    )
    rows.append(
        {
            "group": "nsrdb_mean_MI_bits",
            "S_mean": float(np.nanmean([float(r["mean_MI_bits"]) for r in nsr])),
            "S_median": float(np.nanmedian([float(r["mean_MI_bits"]) for r in nsr])),
            "n_records": len(nsr),
            "unit": "MI_bits",
        }
    )
    rows.append(
        {
            "group": "event_basal_mean_MI_bits",
            "S_mean": float(np.nanmean([float(r["mean_MI_bits_basal"]) for r in ev])),
            "S_median": float(np.nanmedian([float(r["mean_MI_bits_basal"]) for r in ev])),
            "n_records": len(ev),
            "unit": "MI_bits",
        }
    )
    rows.append(
        {
            "group": "nsrdb_mean_rr_cv",
            "S_mean": float(np.nanmean([float(r["mean_rr_cv"]) for r in nsr])),
            "S_median": float(np.nanmedian([float(r["mean_rr_cv"]) for r in nsr])),
            "n_records": len(nsr),
            "unit": "rr_cv",
        }
    )
    rows.append(
        {
            "group": "event_basal_mean_rr_cv",
            "S_mean": float(np.nanmean([float(r["mean_rr_cv_basal"]) for r in ev])),
            "S_median": float(np.nanmedian([float(r["mean_rr_cv_basal"]) for r in ev])),
            "n_records": len(ev),
            "unit": "rr_cv",
        }
    )
    return rows


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(description="Diagnostic surplus S_t analysis (NSRDB + event basals)")
    ap.add_argument("--control-max-hours", type=float, default=CONTROL_MAX_HOURS)
    ap.add_argument("--basal-hours", type=float, default=BASAL_HOURS)
    ap.add_argument("--L", type=int, default=L_SURPLUS)
    ap.add_argument("--high-abs", type=float, default=HIGH_S_ABS)
    ap.add_argument("--cov-stride", type=int, default=25, help="Stride for covariate sampling on NSRDB")
    ap.add_argument("--smoke", action="store_true", help="1 NSRDB + 2 event records only")
    args = ap.parse_args(argv)

    nsrdb_paths = list_nsrdb()
    events = list_event_records()
    if args.smoke:
        nsrdb_paths = nsrdb_paths[:1]
        events = events[:1] + [e for e in events if e[0] == "vfdb"][:1]

    print(f"[surplus-diag] NSRDB records: {len(nsrdb_paths)}")
    print(f"[surplus-diag] Event records: {len(events)}")
    print(f"[surplus-diag] L={args.L} high_abs={args.high_abs} max_h={args.control_max_hours}")

    nsrdb_rows: List[dict] = []
    ep_rows: List[dict] = []
    hvl_rows: List[dict] = []
    for i, p in enumerate(nsrdb_paths, 1):
        print(f"  NSRDB [{i}/{len(nsrdb_paths)}] {p.name} ...", flush=True)
        row, eps, hvl, _ = process_nsrdb_record(
            p,
            max_hours=args.control_max_hours,
            basal_hours=args.basal_hours,
            L=args.L,
            high_abs=args.high_abs,
            cov_stride=args.cov_stride,
        )
        nsrdb_rows.append(row)
        ep_rows.extend(eps)
        hvl_rows.extend(hvl)
        if not row.get("skipped"):
            print(
                f"    mean_S={row['mean_S_all']:.3f} p90={row['p90_S_all']:.3f} "
                f"eps={row['n_episodes_S_ge_high_abs']} rate24={row['episode_rate_per_24h']:.1f}",
                flush=True,
            )

    event_rows: List[dict] = []
    for i, (source, rec, p) in enumerate(events, 1):
        print(f"  EVENT [{i}/{len(events)}] {source}/{rec} ...", flush=True)
        er = process_event_basal(source, rec, p, L=args.L, high_abs=args.high_abs)
        event_rows.append(er)
        if not er.get("skipped"):
            print(
                f"    basal mean_S={er['mean_S_basal']:.3f} p90={er['p90_S_basal']:.3f} "
                f"approach mean_S={er.get('mean_S_approach', float('nan')):.3f}",
                flush=True,
            )

    hvl_agg = aggregate_hvl(hvl_rows)
    comp = cohort_comparison(nsrdb_rows, event_rows)

    # Paths
    p_nsrdb = RES / "surplus_diag_nsrdb_per_record.csv"
    p_eps = RES / "surplus_diag_nsrdb_episodes.csv"
    p_hvl = RES / "surplus_diag_nsrdb_high_vs_low.csv"
    p_hvl_agg = RES / "surplus_diag_nsrdb_high_vs_low_pooled.csv"
    p_ev = RES / "surplus_diag_event_basal_per_record.csv"
    p_comp = RES / "surplus_diag_cohort_comparison.csv"
    p_sum = RES / "surplus_diag_summary.json"

    write_csv(p_nsrdb, nsrdb_rows)
    write_csv(p_eps, ep_rows)
    write_csv(p_hvl, hvl_rows)
    write_csv(p_hvl_agg, hvl_agg)
    write_csv(p_ev, event_rows)
    write_csv(p_comp, comp)

    nsr = [r for r in nsrdb_rows if not r.get("skipped")]
    ev = [r for r in event_rows if not r.get("skipped")]
    means_nsrdb = [float(r["mean_S_all"]) for r in nsr]
    means_ev_basal = [float(r["mean_S_basal"]) for r in ev]
    means_ev_app = [
        float(r["mean_S_approach"])
        for r in ev
        if np.isfinite(float(r.get("mean_S_approach", np.nan)))
    ]

    # Simple two-sample effect (NSRDB mean vs event basal mean of means)
    def cohens_d(a: List[float], b: List[float]) -> float:
        a = np.asarray(a, dtype=float)
        b = np.asarray(b, dtype=float)
        a = a[np.isfinite(a)]
        b = b[np.isfinite(b)]
        if a.size < 2 or b.size < 2:
            return float("nan")
        sp = np.sqrt(
            ((a.size - 1) * np.var(a, ddof=1) + (b.size - 1) * np.var(b, ddof=1))
            / (a.size + b.size - 2)
        )
        return float((np.mean(a) - np.mean(b)) / sp) if sp > 1e-12 else float("nan")

    summary = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "exercise": "surplus_nsrdb_diagnostic",
        "not_a_detector": True,
        "params": {
            "L": args.L,
            "m": M_EMB,
            "delay": DELAY,
            "k1": BP_ALPHABET,
            "k2": BP_ALPHABET,
            "K_joint": K_JOINT,
            "high_abs": args.high_abs,
            "control_max_hours": args.control_max_hours,
            "basal_hours": args.basal_hours,
            "cov_stride": args.cov_stride,
            "min_episode_len": MIN_EPISODE_LEN,
            "merge_gap": MERGE_GAP,
            "definition_S": "TV(P_joint, P1⊗P2) Bandt-Pompe bivariate m=3",
        },
        "nsrdb": {
            "n_records": len(nsr),
            "n_skipped": sum(1 for r in nsrdb_rows if r.get("skipped")),
            "mean_of_mean_S": float(np.mean(means_nsrdb)) if means_nsrdb else float("nan"),
            "median_of_mean_S": float(np.median(means_nsrdb)) if means_nsrdb else float("nan"),
            "mean_of_p90_S": float(np.mean([float(r["p90_S_all"]) for r in nsr])) if nsr else float("nan"),
            "mean_frac_S_ge_0p40": float(np.mean([float(r["frac_S_ge_0p40"]) for r in nsr])) if nsr else float("nan"),
            "mean_frac_S_ge_0p50": float(np.mean([float(r["frac_S_ge_0p50"]) for r in nsr])) if nsr else float("nan"),
            "total_episodes": int(sum(int(r["n_episodes_S_ge_high_abs"]) for r in nsr)),
            "mean_episode_rate_per_24h": float(np.mean([float(r["episode_rate_per_24h"]) for r in nsr])) if nsr else float("nan"),
            "mean_of_mean_S_basal": float(np.mean([float(r["mean_S_basal"]) for r in nsr])) if nsr else float("nan"),
            "mean_of_mean_S_search": float(np.mean([float(r["mean_S_search"]) for r in nsr])) if nsr else float("nan"),
            "mean_diversity": float(np.nanmean([float(r["mean_diversity"]) for r in nsr])) if nsr else float("nan"),
            "mean_MI_bits": float(np.nanmean([float(r["mean_MI_bits"]) for r in nsr])) if nsr else float("nan"),
            "mean_rr_cv": float(np.nanmean([float(r["mean_rr_cv"]) for r in nsr])) if nsr else float("nan"),
        },
        "event_basal": {
            "n_records": len(ev),
            "n_sddb": sum(1 for r in ev if r.get("cohort") == "sddb"),
            "n_vfdb": sum(1 for r in ev if r.get("cohort") == "vfdb"),
            "mean_of_mean_S_basal": float(np.mean(means_ev_basal)) if means_ev_basal else float("nan"),
            "median_of_mean_S_basal": float(np.median(means_ev_basal)) if means_ev_basal else float("nan"),
            "mean_of_p90_S_basal": float(np.mean([float(r["p90_S_basal"]) for r in ev])) if ev else float("nan"),
            "mean_of_mean_S_approach": float(np.mean(means_ev_app)) if means_ev_app else float("nan"),
            "mean_diversity_basal": float(np.nanmean([float(r["mean_diversity_basal"]) for r in ev])) if ev else float("nan"),
            "mean_MI_bits_basal": float(np.nanmean([float(r["mean_MI_bits_basal"]) for r in ev])) if ev else float("nan"),
            "mean_rr_cv_basal": float(np.nanmean([float(r["mean_rr_cv_basal"]) for r in ev])) if ev else float("nan"),
        },
        "contrast": {
            "nsrdb_meanS_minus_event_basal_meanS": (
                float(np.mean(means_nsrdb) - np.mean(means_ev_basal))
                if means_nsrdb and means_ev_basal
                else float("nan")
            ),
            "cohens_d_nsrdb_vs_event_basal": cohens_d(means_nsrdb, means_ev_basal),
            "nsrdb_basal_meanS_minus_event_basal_meanS": (
                float(
                    np.mean([float(r["mean_S_basal"]) for r in nsr])
                    - np.mean(means_ev_basal)
                )
                if nsr and means_ev_basal
                else float("nan")
            ),
        },
        "top_covariates_pooled": hvl_agg[:12],
        "artifacts": {
            "nsrdb_per_record": str(p_nsrdb.name),
            "episodes": str(p_eps.name),
            "high_vs_low": str(p_hvl.name),
            "high_vs_low_pooled": str(p_hvl_agg.name),
            "event_basal": str(p_ev.name),
            "cohort_comparison": str(p_comp.name),
            "summary": str(p_sum.name),
        },
        "smoke": bool(args.smoke),
    }
    p_sum.write_text(json.dumps(summary, indent=2))
    print(f"[surplus-diag] wrote {p_sum}")
    print(
        f"[surplus-diag] NSRDB mean_of_mean_S={summary['nsrdb']['mean_of_mean_S']:.4f} "
        f"| event basal={summary['event_basal']['mean_of_mean_S_basal']:.4f} "
        f"| d={summary['contrast']['cohens_d_nsrdb_vs_event_basal']:.3f}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
