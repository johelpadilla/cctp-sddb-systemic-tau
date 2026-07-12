#!/usr/bin/env python3
"""
Pure CCTP metric helpers: windows, rolling EWS, lead-time detection, detector scoring.

These functions are I/O-free and unit-tested. CLI wrappers load RR/JSON and call them.
Manuscript conventions: W_TAU=101, W_EWS=501, stride=5, bivariate proxy [z(RR), z(|dRR|)].
"""
from __future__ import annotations

from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np

# Paper defaults
W_EWS = 501
W_TAU = 101
STRIDE = 5


def get_event_and_windows(
    record: str,
    t_hr: np.ndarray,
    vfon_hr: float,
) -> Tuple[float, Tuple[float, float], Tuple[float, float]]:
    """Basal / approach windows matching analyze_cctp_pilot.py and RECD scripts."""
    total_h = float(np.nanmax(t_hr))
    if str(record) == "35" or (vfon_hr > total_h - 0.8):
        event_hr = total_h
    else:
        event_hr = float(vfon_hr)
    approach_start = max(0.0, event_hr - 3.0)
    approach_end = event_hr

    if str(record) == "35":
        basal_start, basal_end = 6.0, 16.0
    elif str(record) == "30":
        basal_start, basal_end = 0.5, 3.5
    else:
        b_end = max(3.5, approach_start - 4.0)
        basal_start = max(0.5, b_end - 3.0)
        basal_end = b_end
    return event_hr, (basal_start, basal_end), (approach_start, approach_end)


def build_bivariate_proxy(rr: np.ndarray) -> np.ndarray:
    """Minimal bivariate proxy used for τ_s and RECD: [z(RR), z(|ΔRR|)]."""
    rr = np.asarray(rr, dtype=float)
    drr = np.abs(np.diff(rr, prepend=rr[0]))
    rr_z = (rr - np.mean(rr)) / (np.std(rr) + 1e-12)
    drr_z = (drr - np.mean(drr)) / (np.std(drr) + 1e-12)
    return np.column_stack([rr_z, drr_z])


def rolling_var(x: np.ndarray, w: int = W_EWS, stride: int = STRIDE) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    n = len(x)
    out = np.full(n, np.nan)
    for i in range(w - 1, n, stride):
        out[i] = np.var(x[i - w + 1 : i + 1])
    return out


def rolling_ar1(x: np.ndarray, w: int = W_EWS, stride: int = STRIDE) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    n = len(x)
    out = np.full(n, np.nan)
    for i in range(w - 1, n, stride):
        seg = x[i - w + 1 : i + 1]
        if len(seg) > 2 and np.std(seg) > 1e-9:
            r = np.corrcoef(seg[:-1], seg[1:])[0, 1]
            out[i] = r if not np.isnan(r) else np.nan
    return out


def rolling_mean_series(y: np.ndarray, t: np.ndarray, w_hours: float) -> Tuple[np.ndarray, np.ndarray]:
    """
    Causal rolling mean of a sparse (strided) metric series using a fixed time window (hours).
    Returns (t_valid, rolling_mean) only where enough points exist.
    """
    y = np.asarray(y, dtype=float)
    t = np.asarray(t, dtype=float)
    valid = ~np.isnan(y)
    t_v, y_v = t[valid], y[valid]
    if len(y_v) == 0:
        return t_v, y_v
    out = np.full(len(y_v), np.nan)
    j0 = 0
    for i in range(len(y_v)):
        while t_v[i] - t_v[j0] > w_hours:
            j0 += 1
        if i - j0 + 1 >= 3:
            out[i] = np.mean(y_v[j0 : i + 1])
    return t_v, out


def regime_delta(
    metric: np.ndarray,
    t_hr: np.ndarray,
    basal: Tuple[float, float],
    approach: Tuple[float, float],
) -> Dict[str, float]:
    """Mean(approach) - mean(basal) with sample counts."""
    metric = np.asarray(metric, dtype=float)
    t_hr = np.asarray(t_hr, dtype=float)
    b0, b1 = basal
    a0, a1 = approach
    m_b = metric[(t_hr >= b0) & (t_hr <= b1)]
    m_a = metric[(t_hr >= a0) & (t_hr < a1)]
    m_b = m_b[~np.isnan(m_b)]
    m_a = m_a[~np.isnan(m_a)]
    if len(m_b) == 0 or len(m_a) == 0:
        return {
            "basal_mean": float("nan"),
            "approach_mean": float("nan"),
            "delta": float("nan"),
            "n_basal": float(len(m_b)),
            "n_approach": float(len(m_a)),
        }
    mb, ma = float(np.mean(m_b)), float(np.mean(m_a))
    return {
        "basal_mean": mb,
        "approach_mean": ma,
        "delta": ma - mb,
        "n_basal": float(len(m_b)),
        "n_approach": float(len(m_a)),
    }


def basal_stats(
    metric: np.ndarray,
    t_hr: np.ndarray,
    basal: Tuple[float, float],
) -> Tuple[float, float]:
    """Return (mean, std) of metric in basal window; std is sample std (ddof=1) if n>=2."""
    metric = np.asarray(metric, dtype=float)
    t_hr = np.asarray(t_hr, dtype=float)
    b0, b1 = basal
    m = metric[(t_hr >= b0) & (t_hr <= b1)]
    m = m[~np.isnan(m)]
    if len(m) == 0:
        return float("nan"), float("nan")
    mu = float(np.mean(m))
    sd = float(np.std(m, ddof=1)) if len(m) >= 2 else 0.0
    return mu, sd


def detect_lead_time(
    metric: np.ndarray,
    t_hr: np.ndarray,
    event_hr: float,
    basal: Tuple[float, float],
    *,
    z_threshold: float = 2.0,
    min_consecutive: int = 3,
    use_abs: bool = True,
    search_start_hr: Optional[float] = None,
) -> Dict[str, float]:
    """
    Detect first sustained departure of `metric` from basal before VF.

    A sample at time t is "alarmed" if:
      |m(t) - basal_mean| / (basal_std + eps)  >= z_threshold   (if use_abs)
      or signed z in the direction of the approach delta.

    Lead-time (hours) = event_hr - t_first_alarm when min_consecutive alarmed
    samples occur in the pre-event search window after basal ends.

    Returns dict with lead_time_h, detection_hr, alarmed (0/1), z_at_detection,
    basal_mean, basal_std, n_search_points. Missing detection → lead_time_h = nan, alarmed=0.
    """
    metric = np.asarray(metric, dtype=float)
    t_hr = np.asarray(t_hr, dtype=float)
    mu, sd = basal_stats(metric, t_hr, basal)
    eps = 1e-12
    if not np.isfinite(mu):
        return {
            "lead_time_h": float("nan"),
            "detection_hr": float("nan"),
            "alarmed": 0.0,
            "z_at_detection": float("nan"),
            "basal_mean": float("nan"),
            "basal_std": float("nan"),
            "n_search_points": 0.0,
            "z_threshold": float(z_threshold),
            "min_consecutive": float(min_consecutive),
        }

    b_end = basal[1]
    t0 = b_end if search_start_hr is None else float(search_start_hr)
    # search in (t0, event_hr)
    mask = (t_hr > t0) & (t_hr < event_hr) & ~np.isnan(metric)
    t_s = t_hr[mask]
    m_s = metric[mask]
    n_search = len(m_s)
    if n_search == 0:
        return {
            "lead_time_h": float("nan"),
            "detection_hr": float("nan"),
            "alarmed": 0.0,
            "z_at_detection": float("nan"),
            "basal_mean": mu,
            "basal_std": sd,
            "n_search_points": 0.0,
            "z_threshold": float(z_threshold),
            "min_consecutive": float(min_consecutive),
        }

    z = (m_s - mu) / (sd + eps)
    if use_abs:
        alarm = np.abs(z) >= z_threshold
    else:
        # prefer direction of late-window mean vs basal
        late = m_s[t_s >= (event_hr - 1.0)] if np.any(t_s >= event_hr - 1.0) else m_s
        direction = 1.0 if (np.mean(late) - mu) >= 0 else -1.0
        alarm = (direction * z) >= z_threshold

    run = 0
    det_i = None
    for i, a in enumerate(alarm):
        if a:
            run += 1
            if run >= min_consecutive:
                det_i = i - min_consecutive + 1
                break
        else:
            run = 0

    if det_i is None:
        return {
            "lead_time_h": float("nan"),
            "detection_hr": float("nan"),
            "alarmed": 0.0,
            "z_at_detection": float("nan"),
            "basal_mean": mu,
            "basal_std": sd,
            "n_search_points": float(n_search),
            "z_threshold": float(z_threshold),
            "min_consecutive": float(min_consecutive),
        }

    det_hr = float(t_s[det_i])
    return {
        "lead_time_h": float(event_hr - det_hr),
        "detection_hr": det_hr,
        "alarmed": 1.0,
        "z_at_detection": float(z[det_i]),
        "basal_mean": mu,
        "basal_std": sd,
        "n_search_points": float(n_search),
        "z_threshold": float(z_threshold),
        "min_consecutive": float(min_consecutive),
    }


def detector_performance(
    lead_rows: Sequence[Dict[str, float]],
    *,
    metric_key: str = "lead_time_h",
) -> Dict[str, float]:
    """
    Cohort-level detector summary from per-record lead-time dicts (with 'alarmed').

    On an all-positive pre-VF cohort we report sensitivity (detection rate) and
    median/mean lead-time among detections. False-alarm rate requires control
    Holters (reported as nan here when not provided).
    """
    if not lead_rows:
        return {
            "n_records": 0.0,
            "n_detected": 0.0,
            "sensitivity": float("nan"),
            "median_lead_time_h": float("nan"),
            "mean_lead_time_h": float("nan"),
            "min_lead_time_h": float("nan"),
            "max_lead_time_h": float("nan"),
            "false_alarm_rate": float("nan"),
        }
    alarmed = np.array([float(r.get("alarmed", 0)) for r in lead_rows], dtype=float)
    leads = np.array(
        [float(r.get(metric_key, np.nan)) for r in lead_rows if float(r.get("alarmed", 0)) > 0.5],
        dtype=float,
    )
    leads = leads[np.isfinite(leads)]
    n = float(len(lead_rows))
    n_det = float(np.sum(alarmed > 0.5))
    return {
        "n_records": n,
        "n_detected": n_det,
        "sensitivity": n_det / n if n > 0 else float("nan"),
        "median_lead_time_h": float(np.median(leads)) if len(leads) else float("nan"),
        "mean_lead_time_h": float(np.mean(leads)) if len(leads) else float("nan"),
        "min_lead_time_h": float(np.min(leads)) if len(leads) else float("nan"),
        "max_lead_time_h": float(np.max(leads)) if len(leads) else float("nan"),
        "false_alarm_rate": float("nan"),  # requires negative controls
    }


def cumulative_detection_curve(
    lead_rows: Sequence[Dict[str, float]],
    horizons_h: Sequence[float] = (0.5, 1.0, 2.0, 3.0, 6.0, 12.0),
) -> List[Dict[str, float]]:
    """
    Fraction of records detected with lead_time >= h for each horizon h.
    Records not alarmed contribute 0 at all horizons.
    """
    n = len(lead_rows)
    if n == 0:
        return [{"horizon_h": float(h), "detection_rate": float("nan"), "n_detected": 0.0} for h in horizons_h]
    out = []
    for h in horizons_h:
        det = 0
        for r in lead_rows:
            lt = float(r.get("lead_time_h", np.nan))
            if float(r.get("alarmed", 0)) > 0.5 and np.isfinite(lt) and lt >= float(h):
                det += 1
        out.append(
            {
                "horizon_h": float(h),
                "detection_rate": det / n,
                "n_detected": float(det),
                "n_records": float(n),
            }
        )
    return out


def sign_concordance(a: Sequence[float], b: Sequence[float]) -> Dict[str, float]:
    """Direction concordance: same sign of deltas (zeros count as non-concordant)."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    mask = np.isfinite(a) & np.isfinite(b)
    a, b = a[mask], b[mask]
    if len(a) == 0:
        return {"n": 0.0, "n_concordant": 0.0, "concordance": float("nan")}
    sa = np.sign(a)
    sb = np.sign(b)
    # zeros: treat as non-matching unless both zero
    conc = (sa == sb) & (sa != 0)
    both_zero = (sa == 0) & (sb == 0)
    n_conc = float(np.sum(conc | both_zero))
    return {
        "n": float(len(a)),
        "n_concordant": n_conc,
        "concordance": n_conc / float(len(a)),
    }


def effect_size_cohens_d(x: Sequence[float], y: Sequence[float]) -> float:
    """Cohen's d between two independent samples (basal vs approach style arrays)."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    x = x[np.isfinite(x)]
    y = y[np.isfinite(y)]
    if len(x) < 2 or len(y) < 2:
        return float("nan")
    nx, ny = len(x), len(y)
    vx, vy = np.var(x, ddof=1), np.var(y, ddof=1)
    pooled = np.sqrt(((nx - 1) * vx + (ny - 1) * vy) / (nx + ny - 2) + 1e-12)
    return float((np.mean(y) - np.mean(x)) / pooled)


def substrate_label(rhythm: str, pacing: str) -> str:
    """Coarse substrate stratum for stratified tables."""
    r = (rhythm or "").lower()
    p = (pacing or "none").lower()
    if "fibrillation" in r or r.strip() == "af" or "atrial fib" in r:
        return "AF"
    if p not in ("none", "", "nan") and p != "false":
        return "paced"
    if "pac" in r:
        return "paced"
    return "sinus"


def event_type_label(record: str, event_hr: Optional[float], duration_h: Optional[float]) -> str:
    """Terminal (event near end) vs intermediate (clear pre/post).

    Rules (hours from recording start):
      - intermediate if event is >3 h before end of series (room for post-event context)
      - terminal if event is within 1 h of end, or in the ambiguous 1–3 h band
      - record 30 is always intermediate (manuscript intermediate-event example)
    Never treat duration_h itself as the event time.
    """
    if str(record) == "30":
        return "intermediate"
    if event_hr is None or duration_h is None:
        return "unknown"
    try:
        eh = float(event_hr)
        dh = float(duration_h)
    except (TypeError, ValueError):
        return "unknown"
    if not np.isfinite(eh) or not np.isfinite(dh) or dh <= 0:
        return "unknown"
    gap = dh - eh
    if gap > 3.0:
        return "intermediate"
    return "terminal"


def resolve_event_timing_from_npz(
    record: str,
    npz_path: str,
) -> dict:
    """
    Load cleaned RR npz and return event_hr / duration_h / event_type via
    get_event_and_windows (same alignment as the manuscript pipeline).

    Required npz keys: t_sec, vfon_sec (and optionally total_hours).
    """
    data = np.load(npz_path, allow_pickle=True)
    t_sec = np.asarray(data["t_sec"], dtype=float)
    t_hr = t_sec / 3600.0
    vfon_hr = float(data["vfon_sec"]) / 3600.0
    if "total_hours" in data:
        duration_h = float(data["total_hours"])
    else:
        duration_h = float(np.nanmax(t_hr))
    event_hr, basal, approach = get_event_and_windows(str(record), t_hr, vfon_hr)
    # Prefer series span for event_type gap (matches get_event_and_windows total_h)
    span_h = float(np.nanmax(t_hr)) if len(t_hr) else duration_h
    et = event_type_label(str(record), event_hr, span_h)
    return {
        "event_hr": float(event_hr),
        "duration_h": float(duration_h if np.isfinite(duration_h) else span_h),
        "vfon_hr": float(vfon_hr),
        "span_h": float(span_h),
        "basal_start": float(basal[0]),
        "basal_end": float(basal[1]),
        "approach_start": float(approach[0]),
        "approach_end": float(approach[1]),
        "event_type": et,
    }


def stratified_summary(
    rows: Sequence[dict],
    group_key: str,
    metric_cols: Sequence[str] = ("delta_tau", "delta_excess3", "delta_var", "delta_ar1"),
) -> List[dict]:
    """Group rows by stratum and report mean/median of metric columns + n."""
    from collections import defaultdict

    groups: Dict[str, List[dict]] = defaultdict(list)
    for r in rows:
        groups[str(r.get(group_key, "unknown"))].append(r)
    out = []
    for g, items in sorted(groups.items()):
        rec = {"stratum": g, "group_key": group_key, "n": len(items)}
        for col in metric_cols:
            vals = np.array([float(it[col]) for it in items if col in it and it[col] not in (None, "")], dtype=float)
            vals = vals[np.isfinite(vals)]
            rec[f"{col}_mean"] = float(np.mean(vals)) if len(vals) else float("nan")
            rec[f"{col}_median"] = float(np.median(vals)) if len(vals) else float("nan")
            rec[f"{col}_n"] = int(len(vals))
        out.append(rec)
    return out


# ---------------------------------------------------------------------------
# External-validation helpers (short DBs + FAR on negative controls)
# ---------------------------------------------------------------------------

# Frozen discovery defaults (Jul-12 2026) — do not retune on validation data
FROZEN_THETA3 = 0.08
FROZEN_HIGH_THRESHOLD = 0.65
FROZEN_Z_THRESHOLD = 2.0
FROZEN_MIN_CONSECUTIVE = 3


def short_db_windows(
    event_hr: float,
    total_h: float,
    *,
    min_pre_event_h: float = 0.25,
) -> Tuple[float, Tuple[float, float], Tuple[float, float], str]:
    """
    Basal / approach windows for short external DBs (VFDB/CU, often << 6 h).

    Rules (pre-specified for Phase 1; independent of discovery Holter heuristics):
      - Prefer Holter-style windows when pre-event duration ≥ 6 h (delegate to
        get_event_and_windows with a dummy record id).
      - Else: split pre-event interval [0, event_hr] into thirds:
          basal = first third, approach = last third (middle is buffer).
      - If pre-event < min_pre_event_h, still return windows but label
        duration_stratum = 'too_short'.

    Returns (event_hr, basal, approach, duration_stratum) where
    duration_stratum ∈ {'holter_ge6h', 'short_15_60min', 'short_ge60min', 'too_short'}.
    """
    event_hr = float(event_hr)
    total_h = float(total_h)
    pre = max(0.0, min(event_hr, total_h))
    if pre >= 6.0:
        # Reuse long-Holter window logic (non-special-case records)
        eh, basal, approach = get_event_and_windows("external", np.array([0.0, total_h]), event_hr)
        return eh, basal, approach, "holter_ge6h"
    if pre < float(min_pre_event_h):
        # Degenerate but defined windows for bookkeeping
        b0, b1 = 0.0, max(pre * 0.3, 1e-6)
        a0, a1 = max(pre * 0.7, 0.0), pre
        return event_hr, (b0, b1), (a0, a1), "too_short"
    # thirds of pre-event
    third = pre / 3.0
    basal = (0.0, third)
    approach = (2.0 * third, pre)
    if pre >= 1.0:
        stratum = "short_ge60min"
    else:
        stratum = "short_15_60min"
    return event_hr, basal, approach, stratum


def count_alarm_episodes(
    metric: np.ndarray,
    t_hr: np.ndarray,
    basal: Tuple[float, float],
    *,
    search_start_hr: Optional[float] = None,
    search_end_hr: Optional[float] = None,
    z_threshold: float = FROZEN_Z_THRESHOLD,
    min_consecutive: int = FROZEN_MIN_CONSECUTIVE,
    use_abs: bool = True,
    refractory_h: float = 0.5,
) -> Dict[str, float]:
    """
    Count sustained abs-z alarm episodes on a control (or any) series.

    Basal stats from `basal`; search in (search_start, search_end). After each
    episode a refractory period suppresses re-triggers. Used for FAR on
    negative-control Holters.

    Returns n_episodes, first_alarm_hr, search_hours, alarmed (0/1).
    """
    metric = np.asarray(metric, dtype=float)
    t_hr = np.asarray(t_hr, dtype=float)
    mu, sd = basal_stats(metric, t_hr, basal)
    eps = 1e-12
    if not np.isfinite(mu):
        return {
            "n_episodes": 0.0,
            "first_alarm_hr": float("nan"),
            "search_hours": 0.0,
            "alarmed": 0.0,
            "basal_mean": float("nan"),
            "basal_std": float("nan"),
        }
    b_end = basal[1]
    t0 = b_end if search_start_hr is None else float(search_start_hr)
    t1 = float(np.nanmax(t_hr)) if search_end_hr is None else float(search_end_hr)
    mask = (t_hr > t0) & (t_hr < t1) & ~np.isnan(metric)
    t_s = t_hr[mask]
    m_s = metric[mask]
    if len(m_s) == 0:
        return {
            "n_episodes": 0.0,
            "first_alarm_hr": float("nan"),
            "search_hours": max(0.0, t1 - t0),
            "alarmed": 0.0,
            "basal_mean": mu,
            "basal_std": sd,
        }
    z = (m_s - mu) / (sd + eps)
    alarm = np.abs(z) >= z_threshold if use_abs else z >= z_threshold

    n_episodes = 0
    first_alarm_hr = float("nan")
    run = 0
    refractory_until = -np.inf
    for i, a in enumerate(alarm):
        if t_s[i] < refractory_until:
            run = 0
            continue
        if a:
            run += 1
            if run >= min_consecutive:
                n_episodes += 1
                if not np.isfinite(first_alarm_hr):
                    first_alarm_hr = float(t_s[i - min_consecutive + 1])
                refractory_until = float(t_s[i]) + float(refractory_h)
                run = 0
        else:
            run = 0
    return {
        "n_episodes": float(n_episodes),
        "first_alarm_hr": first_alarm_hr,
        "search_hours": float(max(0.0, t1 - t0)),
        "alarmed": 1.0 if n_episodes > 0 else 0.0,
        "basal_mean": mu,
        "basal_std": sd,
    }


def false_alarm_rate(
    control_rows: Sequence[Dict[str, float]],
    *,
    episode_key: str = "n_episodes",
    hours_key: str = "search_hours",
) -> Dict[str, float]:
    """
    Aggregate FAR from per-control episode counts.

    FAR (alarms / 24 h) = total_episodes / total_search_hours * 24.
    If total_search_hours == 0, returns nan (honest fallback).
    """
    if not control_rows:
        return {
            "n_controls": 0.0,
            "total_episodes": 0.0,
            "total_search_hours": 0.0,
            "far_per_24h": float("nan"),
            "fraction_alarmed": float("nan"),
            "reason": "no_controls",
        }
    eps = [float(r.get(episode_key, 0.0)) for r in control_rows]
    hrs = [float(r.get(hours_key, 0.0)) for r in control_rows]
    total_ep = float(np.nansum(eps))
    total_h = float(np.nansum(hrs))
    alarmed = [float(r.get("alarmed", 1.0 if e > 0 else 0.0)) for e, r in zip(eps, control_rows)]
    if total_h <= 0:
        return {
            "n_controls": float(len(control_rows)),
            "total_episodes": total_ep,
            "total_search_hours": 0.0,
            "far_per_24h": float("nan"),
            "fraction_alarmed": float(np.mean(alarmed)) if alarmed else float("nan"),
            "reason": "insufficient_control_hours",
        }
    return {
        "n_controls": float(len(control_rows)),
        "total_episodes": total_ep,
        "total_search_hours": total_h,
        "far_per_24h": total_ep / total_h * 24.0,
        "fraction_alarmed": float(np.mean(alarmed)),
        "reason": "ok",
    }
