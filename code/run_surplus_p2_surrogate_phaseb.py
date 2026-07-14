#!/usr/bin/env python3
"""
Phase β — Surplus-over-surrogate on P2-10 (diagnostic only).

For each record, compute raw TV surplus S on the frozen P2-10 proxy and a
pointwise null by phase-randomizing (or IAAFT) channel 2 only:

  S^ex_t = S_t − median_k(S^{surr}_{k,t})
  S^z_t  = (S_t − mean_k) / std_k

Compare NSRDB floor and event basal / approach−basal contrast for S vs S^ex
(and report S^z). Not a detector; no I0; does not mutate build_bivariate_proxy.

Outputs under results/:
  surplus_p2_surrogate_nsrdb_per_record.csv
  surplus_p2_surrogate_event_per_record.csv
  surplus_p2_surrogate_summary.csv
  surplus_p2_surrogate_summary.json

Reproduce:
  PYTHONPATH=code python3 code/run_surplus_p2_surrogate_phaseb.py
  PYTHONPATH=code python3 code/run_surplus_p2_surrogate_phaseb.py --smoke

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

from proxy_variants import P2_LAG  # noqa: E402
from run_surplus_proxy_bakeoff import (  # noqa: E402
    BASAL_HOURS,
    CONTROL_MAX_HOURS,
    DELAY,
    HIGH_S_ABS,
    K_JOINT,
    L_SURPLUS,
    M_EMB,
    MERGE_GAP,
    MIN_EPISODE_LEN,
    RES,
    _safe_stats,
    cohens_d,
    control_basal,
    find_episodes,
    list_event_records,
    list_nsrdb,
    load_npz,
    resolve_event_windows,
    write_csv,
)
from surplus_surrogates import (  # noqa: E402
    N_SURR_DEFAULT,
    PROXY_DEFAULT,
    surplus_raw_and_excess,
)

# Pre-registered exploratory screening bars (Phase β; not validation)
FLOOR_EX_DROP_VS_RAW = 0.20  # ≥20% additional drop of NSRDB mean floor vs P2-10 raw
MIN_COHENS_D_EX = 0.60  # d(event basal vs NSRDB) on S^ex should not fall below this

# High-episode thresholds (diagnostic scales; not detector gates)
HIGH_SEX_ABS = 0.05  # exploratory absolute on excess
HIGH_SZ = 2.0  # z over surrogate ensemble at each t

PROXY_ID = PROXY_DEFAULT  # P2-10


def _stable_seed_offset(label: str) -> int:
    """Deterministic non-cryptographic offset from a record label (not Python hash)."""
    h = 5381
    for ch in str(label):
        h = ((h << 5) + h + ord(ch)) & 0xFFFFFFFF
    return int(h % 10_000)


def _mean_of(rows: List[dict], key: str) -> float:
    vals = [float(r[key]) for r in rows if np.isfinite(float(r.get(key, np.nan)))]
    return float(np.mean(vals)) if vals else float("nan")


def process_nsrdb_record(
    path: Path,
    *,
    L: int,
    max_hours: float,
    n_surr: int,
    method: str,
    base_seed: int,
    high_abs: float,
    high_sex: float,
    high_sz: float,
) -> Dict[str, Any]:
    d = load_npz(path)
    rr = np.asarray(d["rr_ms"], dtype=float)
    t_sec = np.asarray(d["t_sec"], dtype=float)
    t_hr = t_sec / 3600.0
    total_h_full = float(d.get("total_hours", float(np.nanmax(t_hr))))
    if max_hours is not None and total_h_full > max_hours:
        mask = t_hr <= max_hours
        rr = rr[mask]
        t_hr = t_hr[mask]
        total_h = float(max_hours)
    else:
        total_h = float(np.nanmax(t_hr)) if t_hr.size else 0.0

    rec_id = path.stem.replace("nsrdb_", "").replace("_clean", "")
    meta: Dict[str, Any] = {
        "proxy_id": PROXY_ID,
        "cohort": "nsrdb",
        "record_id": rec_id,
        "source_path": path.name,
        "total_hours_used": total_h,
        "L": L,
        "m": M_EMB,
        "n_surr": n_surr,
        "surr_method": method,
        "high_abs": high_abs,
        "high_sex": high_sex,
        "high_sz": high_sz,
    }
    if len(rr) < L + 100 or total_h < BASAL_HOURS + 0.5:
        return {**meta, "skipped": True, "skip_reason": "too_short"}

    out = surplus_raw_and_excess(
        rr,
        proxy_id=PROXY_ID,
        L=L,
        n_surr=n_surr,
        method=method,
        base_seed=base_seed + _stable_seed_offset(rec_id),
    )
    S = out["S"]
    S_ex = out["S_ex"]
    S_z = out["S_z"]
    offset = int(out["offset"])
    if len(S) < L + 10:
        return {**meta, "skipped": True, "skip_reason": "too_few_symbols"}

    t_sym = t_hr[offset : offset + len(S)]
    n = min(len(S), len(t_sym), len(S_ex), len(S_z))
    S, S_ex, S_z, t_sym = S[:n], S_ex[:n], S_z[:n], t_sym[:n]
    finite = np.isfinite(S) & np.isfinite(S_ex) & np.isfinite(S_z)

    b0, b1 = control_basal(total_h)
    search_h = float(max(0.0, total_h - b1))

    st_S = _safe_stats(S[finite])
    st_ex = _safe_stats(S_ex[finite])
    st_z = _safe_stats(S_z[finite])

    # Episodes on search segment (after basal)
    def ep_rate(series: np.ndarray, thr: float) -> tuple:
        high = (series >= thr).astype(np.int8)
        high[~finite] = 0
        high[t_sym < b1] = 0
        eps = find_episodes(high, t_sym, series)
        rate = (len(eps) / search_h * 24.0) if search_h > 0 else float("nan")
        return len(eps), rate

    n_ep_S, rate_S = ep_rate(S, high_abs)
    n_ep_ex, rate_ex = ep_rate(S_ex, high_sex)
    n_ep_z, rate_z = ep_rate(S_z, high_sz)

    surr_med_mean = float(np.nanmean(out["S_surr_med"][:n]))
    surr_mean_mean = float(np.nanmean(out["S_surr_mean"][:n]))

    return {
        **meta,
        "skipped": False,
        "skip_reason": "",
        "n_symbols": int(n),
        "basal_start": b0,
        "basal_end": b1,
        "search_hours": search_h,
        # raw S
        "mean_S": st_S["mean"],
        "median_S": st_S["median"],
        "std_S": st_S["std"],
        "p90_S": st_S["p90"],
        "frac_S_ge_0p40": st_S["frac_ge_0p40"],
        "n_episodes_S": n_ep_S,
        "episode_rate_S_per_24h": rate_S,
        # null level
        "mean_S_surr_med": surr_med_mean,
        "mean_S_surr_mean": surr_mean_mean,
        # excess
        "mean_Sex": st_ex["mean"],
        "median_Sex": st_ex["median"],
        "std_Sex": st_ex["std"],
        "p90_Sex": st_ex["p90"],
        "frac_Sex_ge_0": float(np.mean(S_ex[finite] >= 0.0)),
        "frac_Sex_ge_high": float(np.mean(S_ex[finite] >= high_sex)),
        "n_episodes_Sex": n_ep_ex,
        "episode_rate_Sex_per_24h": rate_ex,
        # z
        "mean_Sz": st_z["mean"],
        "median_Sz": st_z["median"],
        "std_Sz": st_z["std"],
        "p90_Sz": st_z["p90"],
        "frac_Sz_ge_2": float(np.mean(S_z[finite] >= high_sz)),
        "n_episodes_Sz": n_ep_z,
        "episode_rate_Sz_per_24h": rate_z,
    }


def process_event_record(
    source: str,
    record: str,
    path: Path,
    *,
    L: int,
    n_surr: int,
    method: str,
    base_seed: int,
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

    meta: Dict[str, Any] = {
        "proxy_id": PROXY_ID,
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
        "n_surr": n_surr,
        "surr_method": method,
    }

    out = surplus_raw_and_excess(
        rr,
        proxy_id=PROXY_ID,
        L=L,
        n_surr=n_surr,
        method=method,
        base_seed=base_seed + _stable_seed_offset(f"{source}_{record}"),
    )
    S = out["S"]
    S_ex = out["S_ex"]
    S_z = out["S_z"]
    offset = int(out["offset"])
    if len(S) < L + 10:
        return {**meta, "skipped": True, "skip_reason": "too_few_symbols"}

    t_sym = t_hr[offset : offset + len(S)]
    n = min(len(S), len(t_sym), len(S_ex), len(S_z))
    S, S_ex, S_z, t_sym = S[:n], S_ex[:n], S_z[:n], t_sym[:n]
    finite = np.isfinite(S) & np.isfinite(S_ex) & np.isfinite(S_z)

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

    def win_stats(series: np.ndarray, mask: np.ndarray) -> Dict[str, float]:
        return _safe_stats(series[mask])

    st_S_b = win_stats(S, basal_mask)
    st_S_a = win_stats(S, approach_mask)
    st_ex_b = win_stats(S_ex, basal_mask)
    st_ex_a = win_stats(S_ex, approach_mask)
    st_z_b = win_stats(S_z, basal_mask)
    st_z_a = win_stats(S_z, approach_mask)

    def delta(a: Dict[str, float], b: Dict[str, float]) -> float:
        if np.isfinite(a["mean"]) and np.isfinite(b["mean"]):
            return float(a["mean"] - b["mean"])
        return float("nan")

    return {
        **meta,
        "skipped": False,
        "skip_reason": "",
        "n_symbols": int(n),
        "n_S_basal": st_S_b["n"],
        "n_S_approach": st_S_a["n"],
        # raw
        "mean_S_basal": st_S_b["mean"],
        "mean_S_approach": st_S_a["mean"],
        "p90_S_basal": st_S_b["p90"],
        "delta_mean_S_approach_minus_basal": delta(st_S_a, st_S_b),
        # excess
        "mean_Sex_basal": st_ex_b["mean"],
        "mean_Sex_approach": st_ex_a["mean"],
        "p90_Sex_basal": st_ex_b["p90"],
        "delta_mean_Sex_approach_minus_basal": delta(st_ex_a, st_ex_b),
        # z
        "mean_Sz_basal": st_z_b["mean"],
        "mean_Sz_approach": st_z_a["mean"],
        "p90_Sz_basal": st_z_b["p90"],
        "delta_mean_Sz_approach_minus_basal": delta(st_z_a, st_z_b),
        "mean_S_surr_med": float(np.nanmean(out["S_surr_med"][:n])),
    }


def aggregate(
    nsrdb_rows: List[dict],
    event_rows: List[dict],
) -> Dict[str, Any]:
    nsr = [r for r in nsrdb_rows if not r.get("skipped")]
    ev = [r for r in event_rows if not r.get("skipped")]

    nsr_mean_S = [float(r["mean_S"]) for r in nsr if np.isfinite(float(r.get("mean_S", np.nan)))]
    nsr_mean_ex = [float(r["mean_Sex"]) for r in nsr if np.isfinite(float(r.get("mean_Sex", np.nan)))]
    nsr_mean_z = [float(r["mean_Sz"]) for r in nsr if np.isfinite(float(r.get("mean_Sz", np.nan)))]
    nsr_frac40 = [
        float(r["frac_S_ge_0p40"]) for r in nsr if np.isfinite(float(r.get("frac_S_ge_0p40", np.nan)))
    ]
    nsr_frac_ex_high = [
        float(r["frac_Sex_ge_high"])
        for r in nsr
        if np.isfinite(float(r.get("frac_Sex_ge_high", np.nan)))
    ]
    nsr_ep_S = [
        float(r["episode_rate_S_per_24h"])
        for r in nsr
        if np.isfinite(float(r.get("episode_rate_S_per_24h", np.nan)))
    ]
    nsr_ep_ex = [
        float(r["episode_rate_Sex_per_24h"])
        for r in nsr
        if np.isfinite(float(r.get("episode_rate_Sex_per_24h", np.nan)))
    ]
    nsr_ep_z = [
        float(r["episode_rate_Sz_per_24h"])
        for r in nsr
        if np.isfinite(float(r.get("episode_rate_Sz_per_24h", np.nan)))
    ]
    nsr_p90_S = [float(r["p90_S"]) for r in nsr if np.isfinite(float(r.get("p90_S", np.nan)))]
    nsr_p90_ex = [float(r["p90_Sex"]) for r in nsr if np.isfinite(float(r.get("p90_Sex", np.nan)))]

    ev_S_b = [
        float(r["mean_S_basal"]) for r in ev if np.isfinite(float(r.get("mean_S_basal", np.nan)))
    ]
    ev_ex_b = [
        float(r["mean_Sex_basal"])
        for r in ev
        if np.isfinite(float(r.get("mean_Sex_basal", np.nan)))
    ]
    ev_z_b = [
        float(r["mean_Sz_basal"])
        for r in ev
        if np.isfinite(float(r.get("mean_Sz_basal", np.nan)))
    ]
    ev_dS = [
        float(r["delta_mean_S_approach_minus_basal"])
        for r in ev
        if np.isfinite(float(r.get("delta_mean_S_approach_minus_basal", np.nan)))
    ]
    ev_dex = [
        float(r["delta_mean_Sex_approach_minus_basal"])
        for r in ev
        if np.isfinite(float(r.get("delta_mean_Sex_approach_minus_basal", np.nan)))
    ]
    ev_dz = [
        float(r["delta_mean_Sz_approach_minus_basal"])
        for r in ev
        if np.isfinite(float(r.get("delta_mean_Sz_approach_minus_basal", np.nan)))
    ]

    d_S = cohens_d(ev_S_b, nsr_mean_S)
    d_ex = cohens_d(ev_ex_b, nsr_mean_ex)
    d_z = cohens_d(ev_z_b, nsr_mean_z)

    mean_S = float(np.mean(nsr_mean_S)) if nsr_mean_S else float("nan")
    mean_ex = float(np.mean(nsr_mean_ex)) if nsr_mean_ex else float("nan")
    # Floor drop of excess vs raw: how much lower is mean |centered| floor?
    # Use signed mean: if S_ex mean is smaller than S mean (typical after subtracting
    # positive surr null), report relative reduction (raw - ex) / |raw|.
    if np.isfinite(mean_S) and abs(mean_S) > 1e-12 and np.isfinite(mean_ex):
        floor_drop_frac = (mean_S - mean_ex) / mean_S
    else:
        floor_drop_frac = float("nan")

    n_app_gt_S = sum(1 for x in ev_dS if x > 0)
    n_app_gt_ex = sum(1 for x in ev_dex if x > 0)

    return {
        "proxy_id": PROXY_ID,
        "n_nsrdb": len(nsr),
        "n_event": len(ev),
        # NSRDB raw
        "nsrdb_mean_of_mean_S": mean_S,
        "nsrdb_mean_of_p90_S": float(np.mean(nsr_p90_S)) if nsr_p90_S else float("nan"),
        "nsrdb_mean_frac_S_ge_0p40": float(np.mean(nsr_frac40)) if nsr_frac40 else float("nan"),
        "nsrdb_mean_episode_rate_S_per_24h": float(np.mean(nsr_ep_S)) if nsr_ep_S else float("nan"),
        "nsrdb_mean_of_mean_S_surr_med": _mean_of(nsr, "mean_S_surr_med"),
        # NSRDB excess
        "nsrdb_mean_of_mean_Sex": mean_ex,
        "nsrdb_mean_of_p90_Sex": float(np.mean(nsr_p90_ex)) if nsr_p90_ex else float("nan"),
        "nsrdb_mean_frac_Sex_ge_high": (
            float(np.mean(nsr_frac_ex_high)) if nsr_frac_ex_high else float("nan")
        ),
        "nsrdb_mean_episode_rate_Sex_per_24h": (
            float(np.mean(nsr_ep_ex)) if nsr_ep_ex else float("nan")
        ),
        "nsrdb_mean_of_mean_Sz": float(np.mean(nsr_mean_z)) if nsr_mean_z else float("nan"),
        "nsrdb_mean_episode_rate_Sz_per_24h": (
            float(np.mean(nsr_ep_z)) if nsr_ep_z else float("nan")
        ),
        "floor_drop_frac_Sex_vs_S": floor_drop_frac,
        # Event raw
        "event_mean_of_mean_S_basal": float(np.mean(ev_S_b)) if ev_S_b else float("nan"),
        "event_mean_delta_S_approach_minus_basal": (
            float(np.mean(ev_dS)) if ev_dS else float("nan")
        ),
        "event_frac_approach_gt_basal_S": (
            float(n_app_gt_S / len(ev_dS)) if ev_dS else float("nan")
        ),
        "cohens_d_event_basal_vs_nsrdb_S": d_S,
        # Event excess
        "event_mean_of_mean_Sex_basal": float(np.mean(ev_ex_b)) if ev_ex_b else float("nan"),
        "event_mean_delta_Sex_approach_minus_basal": (
            float(np.mean(ev_dex)) if ev_dex else float("nan")
        ),
        "event_frac_approach_gt_basal_Sex": (
            float(n_app_gt_ex / len(ev_dex)) if ev_dex else float("nan")
        ),
        "cohens_d_event_basal_vs_nsrdb_Sex": d_ex,
        "diff_event_basal_minus_nsrdb_Sex": (
            float(np.mean(ev_ex_b) - np.mean(nsr_mean_ex))
            if ev_ex_b and nsr_mean_ex
            else float("nan")
        ),
        # Event z
        "event_mean_of_mean_Sz_basal": float(np.mean(ev_z_b)) if ev_z_b else float("nan"),
        "event_mean_delta_Sz_approach_minus_basal": (
            float(np.mean(ev_dz)) if ev_dz else float("nan")
        ),
        "cohens_d_event_basal_vs_nsrdb_Sz": d_z,
    }


def evaluate_screening(summary: Dict[str, Any]) -> Dict[str, Any]:
    """Pre-registered exploratory Phase-β bars (not validation)."""
    drop = float(summary.get("floor_drop_frac_Sex_vs_S", np.nan))
    d_ex = float(summary.get("cohens_d_event_basal_vs_nsrdb_Sex", np.nan))
    delta_ex = float(summary.get("event_mean_delta_Sex_approach_minus_basal", np.nan))
    floor_ok = bool(np.isfinite(drop) and drop >= FLOOR_EX_DROP_VS_RAW)
    d_ok = bool(np.isfinite(d_ex) and d_ex >= MIN_COHENS_D_EX)
    rise_ok = bool(np.isfinite(delta_ex) and delta_ex > 0.0)
    # interesting if floor drops ≥20% AND (d retained OR positive approach rise)
    interesting = bool(floor_ok and (d_ok or rise_ok))
    return {
        "floor_drop_frac_Sex_vs_S": drop,
        "floor_criterion_pass": floor_ok,
        "cohens_d_Sex": d_ex,
        "d_criterion_pass": d_ok,
        "delta_Sex_approach_minus_basal": delta_ex,
        "positive_approach_rise": rise_ok,
        "interesting_screening": interesting,
        "criteria": {
            "floor_drop_ge": FLOOR_EX_DROP_VS_RAW,
            "min_cohens_d": MIN_COHENS_D_EX,
            "approach_delta": "mean(approach-basal) > 0 on S^ex",
            "logic": "floor_ok AND (d_ok OR positive_approach_rise)",
        },
    }


def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser(description="Phase β surplus-over-surrogate on P2-10")
    ap.add_argument("--smoke", action="store_true", help="2 NSRDB + 2 events")
    ap.add_argument("--L", type=int, default=L_SURPLUS)
    ap.add_argument("--max-hours", type=float, default=CONTROL_MAX_HOURS)
    ap.add_argument("--n-surr", type=int, default=N_SURR_DEFAULT)
    ap.add_argument(
        "--method",
        type=str,
        default="phase",
        choices=["phase", "iaaft"],
        help="Channel-2 surrogate: phase (default) or iaaft",
    )
    ap.add_argument("--seed", type=int, default=20260714)
    ap.add_argument("--high-abs", type=float, default=HIGH_S_ABS)
    ap.add_argument("--high-sex", type=float, default=HIGH_SEX_ABS)
    ap.add_argument("--high-sz", type=float, default=HIGH_SZ)
    args = ap.parse_args(argv)

    if args.n_surr < 1:
        print("--n-surr must be >= 1", file=sys.stderr)
        return 2

    nsrdb_paths = list_nsrdb()
    events = list_event_records()
    if args.smoke:
        nsrdb_paths = nsrdb_paths[:2]
        events = events[:2]
        # still enough surrogates for median; allow smaller for speed if user set default
        if args.n_surr == N_SURR_DEFAULT:
            # keep 19 for API fidelity even in smoke (short records only if smoke limits n)
            pass

    print(
        f"[phaseb] proxy={PROXY_ID} (ℓ={P2_LAG}) method={args.method} "
        f"n_surr={args.n_surr} L={args.L}"
    )
    print(f"[phaseb] NSRDB n={len(nsrdb_paths)} events n={len(events)}")
    print(
        f"[phaseb] high_S={args.high_abs} high_Sex={args.high_sex} "
        f"high_Sz={args.high_sz} max_h={args.max_hours}"
    )

    nsrdb_rows: List[dict] = []
    event_rows: List[dict] = []

    for path in nsrdb_paths:
        row = process_nsrdb_record(
            path,
            L=args.L,
            max_hours=args.max_hours,
            n_surr=args.n_surr,
            method=args.method,
            base_seed=args.seed,
            high_abs=args.high_abs,
            high_sex=args.high_sex,
            high_sz=args.high_sz,
        )
        nsrdb_rows.append(row)
        if not row.get("skipped"):
            print(
                f"  nsrdb/{row['record_id']}: S={row['mean_S']:.4f} "
                f"Sex={row['mean_Sex']:.4f} Sz={row['mean_Sz']:.2f} "
                f"surr_med={row['mean_S_surr_med']:.4f} "
                f"f40={row['frac_S_ge_0p40']:.3f}"
            )
        else:
            print(f"  nsrdb/{row['record_id']}: SKIP {row.get('skip_reason')}")

    for source, rec, path in events:
        row = process_event_record(
            source,
            rec,
            path,
            L=args.L,
            n_surr=args.n_surr,
            method=args.method,
            base_seed=args.seed,
        )
        event_rows.append(row)
        if not row.get("skipped"):
            print(
                f"  {source}/{rec}: S_b={row['mean_S_basal']:.4f} "
                f"Sex_b={row['mean_Sex_basal']:.4f} "
                f"ΔS={row['delta_mean_S_approach_minus_basal']:.4f} "
                f"ΔSex={row['delta_mean_Sex_approach_minus_basal']:.4f}"
            )
        else:
            print(f"  {source}/{rec}: SKIP {row.get('skip_reason')}")

    summary = aggregate(nsrdb_rows, event_rows)
    screening = evaluate_screening(summary)

    # Long-form comparison table rows (raw vs excess vs z)
    table_rows = [
        {
            "metric_family": "raw_S",
            "nsrdb_mean": summary["nsrdb_mean_of_mean_S"],
            "nsrdb_p90": summary["nsrdb_mean_of_p90_S"],
            "nsrdb_frac_high": summary["nsrdb_mean_frac_S_ge_0p40"],
            "nsrdb_eps_per_24h": summary["nsrdb_mean_episode_rate_S_per_24h"],
            "event_basal_mean": summary["event_mean_of_mean_S_basal"],
            "delta_approach_basal": summary["event_mean_delta_S_approach_minus_basal"],
            "cohens_d_event_vs_nsrdb": summary["cohens_d_event_basal_vs_nsrdb_S"],
            "frac_approach_gt_basal": summary["event_frac_approach_gt_basal_S"],
        },
        {
            "metric_family": "S_ex",
            "nsrdb_mean": summary["nsrdb_mean_of_mean_Sex"],
            "nsrdb_p90": summary["nsrdb_mean_of_p90_Sex"],
            "nsrdb_frac_high": summary["nsrdb_mean_frac_Sex_ge_high"],
            "nsrdb_eps_per_24h": summary["nsrdb_mean_episode_rate_Sex_per_24h"],
            "event_basal_mean": summary["event_mean_of_mean_Sex_basal"],
            "delta_approach_basal": summary["event_mean_delta_Sex_approach_minus_basal"],
            "cohens_d_event_vs_nsrdb": summary["cohens_d_event_basal_vs_nsrdb_Sex"],
            "frac_approach_gt_basal": summary["event_frac_approach_gt_basal_Sex"],
        },
        {
            "metric_family": "S_z",
            "nsrdb_mean": summary["nsrdb_mean_of_mean_Sz"],
            "nsrdb_p90": float("nan"),
            "nsrdb_frac_high": float("nan"),
            "nsrdb_eps_per_24h": summary["nsrdb_mean_episode_rate_Sz_per_24h"],
            "event_basal_mean": summary["event_mean_of_mean_Sz_basal"],
            "delta_approach_basal": summary["event_mean_delta_Sz_approach_minus_basal"],
            "cohens_d_event_vs_nsrdb": summary["cohens_d_event_basal_vs_nsrdb_Sz"],
            "frac_approach_gt_basal": float("nan"),
        },
    ]

    p_nsr = RES / "surplus_p2_surrogate_nsrdb_per_record.csv"
    p_ev = RES / "surplus_p2_surrogate_event_per_record.csv"
    p_sum = RES / "surplus_p2_surrogate_summary.csv"
    p_json = RES / "surplus_p2_surrogate_summary.json"

    write_csv(p_nsr, nsrdb_rows)
    write_csv(p_ev, event_rows)
    write_csv(p_sum, table_rows)

    payload: Dict[str, Any] = {
        "exercise": "surplus_p2_surrogate_phase_beta",
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "smoke": bool(args.smoke),
        "methods": {
            "proxy": PROXY_ID,
            "p2_lag": P2_LAG,
            "m": M_EMB,
            "delay": DELAY,
            "L": args.L,
            "K": K_JOINT,
            "statistic": "TV(P_joint, P1⊗P2)",
            "surr_method": args.method,
            "surr_target": "channel_2_only (lagged |ΔRR| after z-score in proxy)",
            "n_surr": args.n_surr,
            "S_ex": "S - median_k(S_surr_k) pointwise",
            "S_z": "(S - mean_k) / std_k pointwise over surrogates",
            "high_S": args.high_abs,
            "high_Sex": args.high_sex,
            "high_Sz": args.high_sz,
            "min_episode_len": MIN_EPISODE_LEN,
            "merge_gap": MERGE_GAP,
            "control_max_hours": args.max_hours,
            "seed": args.seed,
            "note": (
                "Diagnostic surplus-over-surrogate on P2-10 only. "
                "Not a detector; not I0; production build_bivariate_proxy untouched."
            ),
        },
        "summary": summary,
        "screening": screening,
        "comparison_table": table_rows,
        "artifacts": {
            "nsrdb_per_record": str(p_nsr.name),
            "event_per_record": str(p_ev.name),
            "summary": str(p_sum.name),
        },
        "disclaimer": "Diagnostic only. No clinical claims. Not a detector.",
    }
    p_json.write_text(json.dumps(payload, indent=2, default=str))

    print("\n[phaseb] === SUMMARY (means of per-record stats) ===")
    print(
        f"{'family':8s} {'nsrdb':>8s} {'p90':>8s} {'f_high':>8s} "
        f"{'eps/24h':>8s} {'evBas':>8s} {'Δapp':>8s} {'d':>6s}"
    )
    for r in table_rows:
        print(
            f"{r['metric_family']:8s} "
            f"{r['nsrdb_mean']:8.4f} "
            f"{r['nsrdb_p90']:8.4f} "
            f"{r['nsrdb_frac_high']:8.3f} "
            f"{r['nsrdb_eps_per_24h']:8.1f} "
            f"{r['event_basal_mean']:8.4f} "
            f"{r['delta_approach_basal']:8.4f} "
            f"{r['cohens_d_event_vs_nsrdb']:6.2f}"
        )
    print(
        f"[phaseb] floor_drop Sex vs S = {summary['floor_drop_frac_Sex_vs_S']:.3f} "
        f"| screening interesting={screening['interesting_screening']}"
    )
    print(f"[phaseb] wrote {p_nsr.name}, {p_ev.name}, {p_sum.name}, {p_json.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
