#!/usr/bin/env python3
"""
Holter evaluation of integrated OPSP detector (I0 primary; I-mod / I-confirm secondary).

Protocol aligned with ordinal exploratory bake-off (sensitivity) and Phase-2 NSRDB FAR:
  - First alarm in post-basal pre-event window → sensitivity (SDDB + VFDB)
  - Episode × 0.5 h refractory → FAR/24h (NSRDB controls, cap 12 h)

Configs (docs/INTEGRATED_N2_PERSIST_N3_SURPLUS_DETECTOR.md):
  I0:        collapse_role=none,     L=50, θ_ΔS=0.08, θ_R^S=5
  I-mod:     collapse_role=modulate, same surplus gates + θ_D=0.35, δ_R=2, δ_S=0.02
  I-confirm: collapse_role=confirm,  θ_ΔS=0.12, θ_R^S=8, θ_D=0.35, w=5  (secondary FAR-leaning)
  OPS-S1:    ops_detect basal-relative same gates as I0 (Mode-S / surplus-persist twin)

Does NOT retune frozen abs-z. Observational only — no clinical claims.

Outputs under results/:
  opsp_integrated_sens_per_record.csv
  opsp_integrated_nsrdb_far_per_record.csv
  opsp_integrated_summary.json
  opsp_integrated_comparison.csv
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
    count_binary_alarm_episodes,
    get_event_and_windows,
    short_db_windows,
    build_bivariate_proxy,
)
from ordinal_detectors.opc_refinements import ops_detect, opsp_integrated_detect
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

# Integrated configs
I0 = dict(
    name="I0",
    collapse_role="none",
    L=50,
    theta_delta_S=0.08,
    theta_R=5,
    theta_D=0.35,
    confirm_window=5,
    modulate_delta_R=2,
    modulate_delta_S=0.02,
)
I_MOD = dict(
    name="I-mod",
    collapse_role="modulate",
    L=50,
    theta_delta_S=0.08,
    theta_R=5,
    theta_D=0.35,
    confirm_window=5,
    modulate_delta_R=2,
    modulate_delta_S=0.02,
)
I_CONFIRM = dict(
    name="I-confirm",
    collapse_role="confirm",
    L=50,
    theta_delta_S=0.12,
    theta_R=8,
    theta_D=0.35,
    confirm_window=5,
    modulate_delta_R=2,
    modulate_delta_S=0.02,
)
OPS_S1 = dict(
    name="OPS-S1",
    L=50,
    theta_delta_S=0.08,
    theta_R=5,
)

# Phase-2 control window defaults
BASAL_HOURS = 2.0
CONTROL_MAX_HOURS = 12.0
REFRACTORY_H = 0.5

# Published baseline anchors (same protocol family; not re-invented)
BASELINE_ANCHORS = {
    "opc_L50": {
        "sens_sddb": 0.5454545455,
        "sens_vfdb": 0.3636363636,
        "sens_all": 0.4242424242,
        "n_sddb": 11,
        "n_detected_sddb": 6,
        "n_vfdb": 22,
        "n_detected_vfdb": 8,
        "n_all": 33,
        "n_detected_all": 14,
        "far_per_24h": 3.733376994,
        "total_episodes": 28,
        "total_search_hours": 179.997895,
        "n_controls": 18,
        "source": "results/ordinal_sensitivity_specificity_tradeoff.csv + ordinal_nsrdb_far_summary.json",
    },
    "absz_tau_s": {
        "sens_sddb": 1.0,
        "sens_vfdb": 0.8636363636,
        "sens_all": 0.9090909091,
        "n_sddb": 11,
        "n_detected_sddb": 11,
        "n_vfdb": 22,
        "n_detected_vfdb": 19,
        "n_all": 33,
        "n_detected_all": 30,
        "far_per_24h": 33.73372784,
        "total_episodes": 253,
        "total_search_hours": 179.997895,
        "n_controls": 18,
        "source": "results/ordinal_exploratory_summary.json + ordinal_nsrdb_far_summary.json",
    },
    "sdd": {
        "sens_sddb": 1.0,
        "sens_vfdb": 0.9545454545,
        "sens_all": 0.9696969697,
        "far_per_24h": 46.26720774,
        "source": "results/ordinal_sensitivity_specificity_tradeoff.csv",
    },
}


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
    """Return (pi1, pi2, offset) Bandt–Pompe factors; empty if too short."""
    X = build_bivariate_proxy(np.asarray(rr, dtype=float))
    S = generate_multivariate_symbols(X, m=M_EMB, delay=DELAY)
    offset = (M_EMB - 1) * DELAY
    if S.size == 0 or S.shape[1] < 2:
        return np.array([], dtype=np.int64), np.array([], dtype=np.int64), offset
    pi1 = S[:, 0].astype(np.int64)
    pi2 = S[:, 1].astype(np.int64)
    return pi1, pi2, offset


def hours_to_index_left(t_hr: np.ndarray, hr: float) -> int:
    if len(t_hr) == 0:
        return 0
    return max(0, min(int(np.searchsorted(t_hr, hr, side="left")), len(t_hr)))


def hours_to_index_right(t_hr: np.ndarray, hr: float) -> int:
    if len(t_hr) == 0:
        return 0
    return max(0, min(int(np.searchsorted(t_hr, hr, side="left")), len(t_hr)))


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


def control_basal(total_h: float, basal_hours: float = BASAL_HOURS) -> Tuple[float, float]:
    basal = (0.25, min(basal_hours, total_h * 0.25))
    if basal[1] <= basal[0]:
        basal = (0.0, max(total_h * 0.2, 0.1))
    return basal


def first_alarm_in_range(
    alarm: np.ndarray, search_start: int, search_end: int
) -> Optional[int]:
    """First index with alarm==1 in [search_start, search_end)."""
    lo = max(0, int(search_start))
    hi = min(len(alarm), int(search_end))
    if hi <= lo:
        return None
    hits = np.flatnonzero(alarm[lo:hi] > 0)
    if hits.size == 0:
        return None
    return int(lo + hits[0])


def run_opsp(
    pi1: np.ndarray,
    pi2: np.ndarray,
    cfg: dict,
    basal_end: int,
) -> Dict[str, np.ndarray]:
    """Run opsp_integrated_detect; optional structural keys: basal_mode, basal_q, mad_kappa, hop."""
    kwargs: Dict[str, Any] = dict(
        L=int(cfg["L"]),
        theta_R=int(cfg["theta_R"]),
        k1=BP_ALPHABET,
        k2=BP_ALPHABET,
        basal_end=int(basal_end),
        collapse_role=str(cfg.get("collapse_role", "none")),
        theta_D=float(cfg.get("theta_D", 0.35)),
        confirm_window=int(cfg.get("confirm_window", 5)),
        modulate_delta_R=int(cfg.get("modulate_delta_R", 2)),
        modulate_delta_S=float(cfg.get("modulate_delta_S", 0.02)),
        basal_mode=str(cfg.get("basal_mode", "mean_delta")),
        basal_q=float(cfg.get("basal_q", 90.0)),
        mad_kappa=float(cfg.get("mad_kappa", 2.5)),
        hop=int(cfg.get("hop", 1)),
    )
    if cfg.get("theta_delta_S") is not None:
        kwargs["theta_delta_S"] = float(cfg["theta_delta_S"])
    if cfg.get("theta_S") is not None:
        kwargs["theta_S"] = float(cfg["theta_S"])
    return opsp_integrated_detect(pi1, pi2, **kwargs)


def run_ops_s1(
    pi1: np.ndarray,
    pi2: np.ndarray,
    basal_end: int,
    cfg: dict = OPS_S1,
) -> Dict[str, np.ndarray]:
    kwargs: Dict[str, Any] = dict(
        L=int(cfg["L"]),
        theta_R=int(cfg["theta_R"]),
        k1=BP_ALPHABET,
        k2=BP_ALPHABET,
        basal_end=int(basal_end),
        basal_mode=str(cfg.get("basal_mode", "mean_delta")),
        basal_q=float(cfg.get("basal_q", 90.0)),
        mad_kappa=float(cfg.get("mad_kappa", 2.5)),
        hop=int(cfg.get("hop", 1)),
    )
    if cfg.get("theta_delta_S") is not None:
        kwargs["theta_delta_S"] = float(cfg["theta_delta_S"])
    if cfg.get("theta_S") is not None:
        kwargs["theta_S"] = float(cfg["theta_S"])
    return ops_detect(pi1, pi2, **kwargs)


def process_event_record(
    source: str,
    record: str,
    path: Path,
    configs: List[dict],
) -> List[dict]:
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
    a0, a1 = approach

    pi1, pi2, offset = bivariate_factors(rr)
    rows: List[dict] = []
    base_meta = {
        "source": source,
        "record": record,
        "event_label": event_label,
        "duration_stratum": stratum,
        "event_hr": event_hr,
        "basal_start": b0,
        "basal_end": b1,
        "approach_start": a0,
        "approach_end": a1,
        "path": path.name,
    }
    if len(pi1) == 0:
        for cfg in configs:
            rows.append(
                {
                    **base_meta,
                    "detector": cfg["name"],
                    "error": "empty_symbol_stream",
                    "alarmed": 0,
                    "lead_time_h": float("nan"),
                }
            )
        return rows

    t_hr_sym = t_hr[offset : offset + len(pi1)]
    n = min(len(pi1), len(t_hr_sym))
    pi1, pi2, t_hr_sym = pi1[:n], pi2[:n], t_hr_sym[:n]

    basal_i0 = hours_to_index_left(t_hr_sym, b0)
    basal_i1 = hours_to_index_right(t_hr_sym, b1)
    if basal_i1 <= basal_i0:
        event_i = hours_to_index_right(t_hr_sym, event_hr)
        basal_i0 = 0
        basal_i1 = max(1, event_i // 3)

    # OPS needs basal_end >= L
    L_need = max(int(c["L"]) for c in configs)
    if basal_i1 < L_need:
        basal_i1 = min(len(pi1), max(L_need, basal_i1))

    search_start = hours_to_index_left(t_hr_sym, b1)
    search_end = hours_to_index_right(t_hr_sym, event_hr)

    for cfg in configs:
        name = cfg["name"]
        try:
            if name == "OPS-S1":
                out = run_ops_s1(pi1, pi2, basal_i1, cfg)
                core = out["alarm"]
            else:
                out = run_opsp(pi1, pi2, cfg, basal_i1)
                core = out.get("core_alarm", out["alarm"])
            alarm = out["alarm"]
            idx = first_alarm_in_range(alarm, search_start, search_end)
            if idx is None:
                lead, det_hr, alarmed = float("nan"), float("nan"), 0
            else:
                det_hr = float(t_hr_sym[idx])
                lead = float(event_hr - det_hr)
                alarmed = 1

            # Diagnostics: surplus / collapse in approach window
            app_mask = (t_hr_sym >= a0) & (t_hr_sym < a1)
            syn = out.get("synergy")
            div = out.get("diversity")
            low = out.get("low_div")
            high = out.get("high_syn")
            mean_S_app = (
                float(np.nanmean(syn[app_mask]))
                if syn is not None and app_mask.any()
                else float("nan")
            )
            mean_D_app = (
                float(np.nanmean(div[app_mask]))
                if div is not None and app_mask.any()
                else float("nan")
            )
            frac_low_app = (
                float(np.mean(low[app_mask]))
                if low is not None and app_mask.any()
                else float("nan")
            )
            frac_high_syn_app = (
                float(np.mean(high[app_mask]))
                if high is not None and app_mask.any()
                else float("nan")
            )
            n_core_search = int(np.sum(core[search_start:search_end] > 0)) if core is not None else -1
            n_alarm_search = int(np.sum(alarm[search_start:search_end] > 0))

            rows.append(
                {
                    **base_meta,
                    "detector": name,
                    "collapse_role": cfg.get("collapse_role", "n/a"),
                    "L": cfg["L"],
                    "theta_delta_S": cfg.get("theta_delta_S", ""),
                    "theta_R": cfg["theta_R"],
                    "basal_mode": cfg.get("basal_mode", "mean_delta"),
                    "basal_q": cfg.get("basal_q", ""),
                    "mad_kappa": cfg.get("mad_kappa", ""),
                    "hop": cfg.get("hop", 1),
                    "n_symbols": int(len(pi1)),
                    "basal_idx_start": int(basal_i0),
                    "basal_idx_end": int(basal_i1),
                    "search_start_idx": int(search_start),
                    "search_end_idx": int(search_end),
                    "alarmed": alarmed,
                    "lead_time_h": lead,
                    "detection_hr": det_hr,
                    "first_alarm_idx": idx if idx is not None else "",
                    "n_alarm_samples_search": n_alarm_search,
                    "n_core_alarm_samples_search": n_core_search,
                    "mean_S_approach": mean_S_app,
                    "mean_D_approach": mean_D_app,
                    "frac_low_div_approach": frac_low_app,
                    "frac_high_syn_approach": frac_high_syn_app,
                    "error": "",
                }
            )
        except Exception as exc:  # noqa: BLE001 — record-level isolation
            rows.append(
                {
                    **base_meta,
                    "detector": name,
                    "collapse_role": cfg.get("collapse_role", "n/a"),
                    "error": str(exc),
                    "alarmed": 0,
                    "lead_time_h": float("nan"),
                }
            )
    return rows


def process_control(
    path: Path,
    configs: List[dict],
    *,
    basal_hours: float = BASAL_HOURS,
    max_hours: float = CONTROL_MAX_HOURS,
    refractory_h: float = REFRACTORY_H,
) -> List[dict]:
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

    rows: List[dict] = []
    empty_base = {
        "record_id": rec,
        "source_path": path.name,
        "total_hours_used": total_h,
        "skipped": True,
    }
    if len(rr) < 200 or total_h < basal_hours + 0.5:
        for cfg in configs:
            rows.append(
                {
                    **empty_base,
                    "detector": cfg["name"],
                    "skip_reason": "too_short",
                }
            )
        return rows

    basal = control_basal(total_h, basal_hours)
    b0, b1 = basal
    search_start = b1
    search_end = total_h

    pi1, pi2, offset = bivariate_factors(rr)
    if len(pi1) < 60:
        for cfg in configs:
            rows.append(
                {
                    **empty_base,
                    "detector": cfg["name"],
                    "skip_reason": "symbol_stream_too_short",
                }
            )
        return rows

    t_hr_sym = t_hr[offset : offset + len(pi1)]
    n = min(len(pi1), len(t_hr_sym))
    pi1, pi2, t_hr_sym = pi1[:n], pi2[:n], t_hr_sym[:n]

    basal_i0 = hours_to_index_left(t_hr_sym, b0)
    basal_i1 = hours_to_index_right(t_hr_sym, b1)
    L_need = max(int(c["L"]) for c in configs)
    if basal_i1 - basal_i0 < max(2, L_need // 5):
        basal_i0 = 0
        basal_i1 = max(L_need, min(len(pi1) // 5, len(pi1) // 4))
        if basal_i1 <= basal_i0:
            basal_i1 = min(len(pi1), basal_i0 + L_need)
    if basal_i1 < L_need:
        basal_i1 = min(len(pi1), L_need)

    for cfg in configs:
        name = cfg["name"]
        try:
            if name == "OPS-S1":
                out = run_ops_s1(pi1, pi2, basal_i1, cfg)
            else:
                out = run_opsp(pi1, pi2, cfg, basal_i1)
            ep = count_binary_alarm_episodes(
                out["alarm"],
                t_hr_sym,
                search_start_hr=search_start,
                search_end_hr=search_end,
                refractory_h=refractory_h,
            )
            h = float(ep["search_hours"])
            far = float(ep["n_episodes"]) / h * 24.0 if h > 0 else float("nan")
            rows.append(
                {
                    "record_id": rec,
                    "source_path": path.name,
                    "skipped": False,
                    "skip_reason": "",
                    "detector": name,
                    "collapse_role": cfg.get("collapse_role", "n/a"),
                    "total_hours_used": total_h,
                    "basal_start": b0,
                    "basal_end": b1,
                    "search_start": search_start,
                    "search_end": search_end,
                    "refractory_h": refractory_h,
                    "n_symbols": int(len(pi1)),
                    "basal_idx_start": int(basal_i0),
                    "basal_idx_end": int(basal_i1),
                    "L": cfg["L"],
                    "theta_delta_S": cfg.get("theta_delta_S", ""),
                    "theta_R": cfg["theta_R"],
                    "basal_mode": cfg.get("basal_mode", "mean_delta"),
                    "basal_q": cfg.get("basal_q", ""),
                    "mad_kappa": cfg.get("mad_kappa", ""),
                    "hop": cfg.get("hop", 1),
                    "n_episodes": ep["n_episodes"],
                    "search_hours": ep["search_hours"],
                    "first_alarm_hr": ep["first_alarm_hr"],
                    "alarmed": ep["alarmed"],
                    "far_per_24h": far,
                    "error": "",
                }
            )
        except Exception as exc:  # noqa: BLE001
            rows.append(
                {
                    **empty_base,
                    "detector": name,
                    "skip_reason": f"error:{exc}",
                    "skipped": True,
                }
            )
    return rows


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


def sens_block(rows: List[dict], detector: str, source: Optional[str] = None) -> dict:
    sub = [r for r in rows if r.get("detector") == detector and not r.get("error")]
    if source is not None:
        sub = [r for r in sub if r.get("source") == source]
    n = len(sub)
    n_det = sum(int(r.get("alarmed") or 0) for r in sub)
    leads = [
        float(r["lead_time_h"])
        for r in sub
        if int(r.get("alarmed") or 0) == 1 and np.isfinite(float(r.get("lead_time_h", np.nan)))
    ]
    return {
        "n": n,
        "n_detected": n_det,
        "sensitivity": float(n_det / n) if n else float("nan"),
        "median_lead_h": float(np.median(leads)) if leads else float("nan"),
        "mean_lead_h": float(np.mean(leads)) if leads else float("nan"),
    }


def far_block(rows: List[dict], detector: str) -> dict:
    sub = [r for r in rows if r.get("detector") == detector and not r.get("skipped")]
    n = len(sub)
    if n == 0:
        return {
            "n_controls": 0,
            "total_episodes": 0.0,
            "total_search_hours": 0.0,
            "far_per_24h": float("nan"),
            "fraction_alarmed": float("nan"),
            "median_far_per_record": float("nan"),
            "mean_far_per_record": float("nan"),
        }
    tot_ep = sum(float(r.get("n_episodes") or 0) for r in sub)
    tot_h = sum(float(r.get("search_hours") or 0) for r in sub)
    fars = [float(r["far_per_24h"]) for r in sub if np.isfinite(float(r.get("far_per_24h", np.nan)))]
    n_alarmed = sum(int(r.get("alarmed") or 0) for r in sub)
    return {
        "n_controls": n,
        "total_episodes": tot_ep,
        "total_search_hours": tot_h,
        "far_per_24h": float(tot_ep / tot_h * 24.0) if tot_h > 0 else float("nan"),
        "fraction_alarmed": float(n_alarmed / n),
        "median_far_per_record": float(np.median(fars)) if fars else float("nan"),
        "mean_far_per_record": float(np.mean(fars)) if fars else float("nan"),
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


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(description="OPSP integrated Holter evaluation (I0 / I-mod / I-confirm)")
    ap.add_argument("--control-max-hours", type=float, default=CONTROL_MAX_HOURS)
    ap.add_argument("--basal-hours", type=float, default=BASAL_HOURS)
    ap.add_argument("--refractory-h", type=float, default=REFRACTORY_H)
    ap.add_argument(
        "--smoke",
        action="store_true",
        help="Process only first 2 SDDB + 2 VFDB + 2 NSRDB (fast path check)",
    )
    ap.add_argument(
        "--scratch-dir",
        type=str,
        default="",
        help="Optional directory to copy logs/summaries",
    )
    args = ap.parse_args(argv)

    configs = [I0, I_MOD, I_CONFIRM, OPS_S1]
    events = list_event_records()
    controls = list_nsrdb()
    if args.smoke:
        sddb = [e for e in events if e[0] == "sddb"][:2]
        vfdb = [e for e in events if e[0] == "vfdb"][:2]
        events = sddb + vfdb
        controls = controls[:2]
        print(f"SMOKE mode: {len(events)} events, {len(controls)} controls")

    print(f"Event records: {len(events)}  NSRDB controls: {len(controls)}")
    print(f"Configs: {[c['name'] for c in configs]}")

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

    # Summaries
    det_names = [c["name"] for c in configs]
    summary: Dict[str, Any] = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "phase": "opsp_integrated_holter_eval",
        "smoke": bool(args.smoke),
        "clinical_claim": False,
        "params": {
            "encoding": "bivariate_BandtPompe_m3_factors",
            "K_joint": K_JOINT,
            "I0": {k: v for k, v in I0.items()},
            "I-mod": {k: v for k, v in I_MOD.items()},
            "I-confirm": {k: v for k, v in I_CONFIRM.items()},
            "OPS-S1": {k: v for k, v in OPS_S1.items()},
            "control_max_hours": args.control_max_hours,
            "basal_hours": args.basal_hours,
            "refractory_h": args.refractory_h,
        },
        "methodology": {
            "sensitivity": (
                "First alarm in post-basal pre-event search window "
                "(same family as ordinal exploratory bake-off)"
            ),
            "far": (
                "FAR = total_episodes / total_search_hours * 24; "
                "episode refractory 0.5 h (count_binary_alarm_episodes); "
                "Phase-2 NSRDB control basal/search"
            ),
            "primary_arm": "I0 (collapse_role=none; surplus-persist primary)",
            "secondary_arms": ["I-mod", "I-confirm"],
            "ops_note": (
                "OPS-S1 uses ops_detect with same basal-relative gates as I0; "
                "I0 core is ontologically the same surplus-persist path"
            ),
            "proxy_label": "S_t = TV(joint, product margins) is Level-3–consistent proxy, not excess3",
        },
        "sensitivity": {},
        "far": {},
        "baseline_anchors": BASELINE_ANCHORS,
    }

    for name in det_names:
        summary["sensitivity"][name] = {
            "sddb": sens_block(sens_rows, name, "sddb"),
            "vfdb": sens_block(sens_rows, name, "vfdb"),
            "all_events": sens_block(sens_rows, name, None),
        }
        summary["far"][name] = far_block(far_rows, name)

    # Comparison table rows
    comparison: List[dict] = []
    for name in det_names:
        s = summary["sensitivity"][name]
        f = summary["far"][name]
        role = "primary" if name == "I0" else (
            "mode_s_twin" if name == "OPS-S1" else "secondary"
        )
        comparison.append(
            {
                "detector": name,
                "role": role,
                "source": "this_run",
                "sens_sddb": s["sddb"]["sensitivity"],
                "n_sddb": s["sddb"]["n"],
                "n_detected_sddb": s["sddb"]["n_detected"],
                "sens_vfdb": s["vfdb"]["sensitivity"],
                "n_vfdb": s["vfdb"]["n"],
                "n_detected_vfdb": s["vfdb"]["n_detected"],
                "sens_all": s["all_events"]["sensitivity"],
                "n_all": s["all_events"]["n"],
                "n_detected_all": s["all_events"]["n_detected"],
                "median_lead_h_all": s["all_events"]["median_lead_h"],
                "far_per_24h": f["far_per_24h"],
                "total_episodes": f["total_episodes"],
                "total_search_hours": f["total_search_hours"],
                "n_controls": f["n_controls"],
                "fraction_controls_alarmed": f["fraction_alarmed"],
            }
        )
    for bname, b in BASELINE_ANCHORS.items():
        if bname == "sdd":
            continue  # high-FAR contrast only if needed; keep table focused
        comparison.append(
            {
                "detector": bname,
                "role": "baseline_anchor",
                "source": b.get("source", "existing_results"),
                "sens_sddb": b.get("sens_sddb"),
                "n_sddb": b.get("n_sddb"),
                "n_detected_sddb": b.get("n_detected_sddb"),
                "sens_vfdb": b.get("sens_vfdb"),
                "n_vfdb": b.get("n_vfdb"),
                "n_detected_vfdb": b.get("n_detected_vfdb"),
                "sens_all": b.get("sens_all"),
                "n_all": b.get("n_all"),
                "n_detected_all": b.get("n_detected_all"),
                "median_lead_h_all": "",
                "far_per_24h": b.get("far_per_24h"),
                "total_episodes": b.get("total_episodes"),
                "total_search_hours": b.get("total_search_hours"),
                "n_controls": b.get("n_controls"),
                "fraction_controls_alarmed": "",
            }
        )

    summary["comparison"] = comparison

    # Paths
    sens_path = RES / "opsp_integrated_sens_per_record.csv"
    far_path = RES / "opsp_integrated_nsrdb_far_per_record.csv"
    sum_path = RES / "opsp_integrated_summary.json"
    cmp_path = RES / "opsp_integrated_comparison.csv"

    write_csv(sens_path, sens_rows)
    write_csv(far_path, far_rows)
    write_csv(cmp_path, comparison)
    sum_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print("\n=== SUMMARY ===")
    for name in det_names:
        s = summary["sensitivity"][name]["all_events"]
        f = summary["far"][name]
        print(
            f"  {name}: sens_all={s['sensitivity']:.3f} "
            f"({s['n_detected']}/{s['n']})  "
            f"FAR={f['far_per_24h']:.3f}/24h  "
            f"ep={f['total_episodes']}"
        )
    print(f"\nWrote:\n  {sens_path}\n  {far_path}\n  {sum_path}\n  {cmp_path}")

    if args.scratch_dir:
        scratch = Path(args.scratch_dir)
        scratch.mkdir(parents=True, exist_ok=True)
        for p in (sens_path, far_path, sum_path, cmp_path):
            (scratch / p.name).write_bytes(p.read_bytes())
        # also write a short text log
        log = scratch / "opsp_eval_console_summary.txt"
        lines = [
            f"generated_at={summary['generated_at']}",
            f"smoke={args.smoke}",
            f"n_events={len(events)} n_controls={len(controls)}",
        ]
        for name in det_names:
            s = summary["sensitivity"][name]
            f = summary["far"][name]
            lines.append(
                f"{name}: sddb={s['sddb']['sensitivity']:.4f} "
                f"vfdb={s['vfdb']['sensitivity']:.4f} "
                f"all={s['all_events']['sensitivity']:.4f} "
                f"FAR={f['far_per_24h']}"
            )
        log.write_text("\n".join(lines) + "\n", encoding="utf-8")
        print(f"Scratch copies → {scratch}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
