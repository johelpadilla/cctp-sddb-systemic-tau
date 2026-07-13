"""
Exploratory light cascade fusion: SDD candidate → OPC L=50 confirmation.

Pure post-processing of two binary alarm series on a shared timebase.
Does NOT retune or rewrite opc_detect / sdd_detect / abs-z.

Rule
----
1. SDD is the high-sensitivity candidate filter.
2. An SDD alarm at time t becomes a cascade alarm only if OPC also alarms
   within a closed confirmation window [t − W, t + W].
3. Default W = ±5 minutes (W = 5/60 h). Local confirmation only — not a
   second independent long-horizon detector.
4. **Causality (event evaluation):** only OPC samples with t_opc < event_hr
   may confirm; the cascade *decision time* is
   decision_time = max(t_SDD, t_OPC_confirm) and must be < event_hr.
   Post-event OPC must never credit a pre-event cascade detection (no look-ahead).

Justification for ±5 min (documented for the exploratory experiment):
- VFDB pre-event horizons are often short (<15–60 min); a tight window avoids
  pairing SDD candidates with unrelated distant OPC spikes.
- Multi-hour SDDB leads still only need *local* co-occurrence to confirm a
  structural (SDD) candidate with an ordinal-persistence (OPC) signature.
- ±5 min is interpretable as “same short-term instability epoch,” not
  independent dual detection.

Exploratory only — no clinical / superiority / S5 claims.
"""

from __future__ import annotations

from typing import Any, Dict, Optional, Tuple

import numpy as np

# Primary confirmation half-window (hours). Closed interval [t-W, t+W].
CONFIRM_WINDOW_MIN = 5.0
CONFIRM_WINDOW_H = CONFIRM_WINDOW_MIN / 60.0  # ±5 minutes


def cascade_sdd_confirm_opc(
    sdd_alarm: np.ndarray,
    opc_alarm: np.ndarray,
    t_hr: np.ndarray,
    *,
    confirm_window_h: float = CONFIRM_WINDOW_H,
    opc_max_hr: Optional[float] = None,
) -> Dict[str, Any]:
    """
    Light cascade: SDD candidates confirmed by OPC within ±confirm_window_h.

    Parameters
    ----------
    sdd_alarm, opc_alarm : array-like of shape (T,)
        Binary (or nonzero-as-true) alarm streams from sdd_detect / opc_detect.
        Must share the same length and timebase as t_hr.
    t_hr : array-like of shape (T,)
        Sample times in hours (symbol-endpoint timebase).
    confirm_window_h : float
        Half-window in hours (default 5/60 ≈ ±5 min). Closed interval.
    opc_max_hr : float, optional
        If set, only OPC samples with t_opc < opc_max_hr are eligible to confirm
        (strict pre-horizon; use event_hr for causal event evaluation).

    Returns
    -------
    dict with keys:
        alarm : (T,) int8 cascade alarms (1 only at SDD-on indices confirmed by OPC)
        confirmed : same as alarm (alias for readability)
        decision_time_hr : (T,) float — max(t_SDD, earliest eligible OPC in window);
            NaN where no cascade
        sdd_candidate_count : int
        opc_alarm_count : int  (eligible OPC count after opc_max_hr filter)
        cascade_alarm_count : int
        confirm_window_h : float
        opc_max_hr : float or None
    """
    if confirm_window_h < 0:
        raise ValueError("confirm_window_h must be >= 0")

    sdd = np.asarray(sdd_alarm)
    opc = np.asarray(opc_alarm)
    t = np.asarray(t_hr, dtype=float)

    if sdd.shape != opc.shape or sdd.shape != t.shape:
        raise ValueError("sdd_alarm, opc_alarm, and t_hr must have the same shape")
    if sdd.ndim != 1:
        raise ValueError("inputs must be 1-D")

    cascade = np.zeros(sdd.shape, dtype=np.int8)
    decision_time = np.full(sdd.shape, np.nan, dtype=float)
    sdd_on = _as_on(sdd)
    opc_on = _as_on(opc)

    sdd_idx = np.flatnonzero(sdd_on)
    opc_idx = np.flatnonzero(opc_on)
    n_sdd = int(sdd_idx.size)

    # Eligible OPC times (causal filter)
    if opc_idx.size:
        opc_times_all = t[opc_idx].astype(float)
        finite = np.isfinite(opc_times_all)
        opc_times_all = opc_times_all[finite]
        if opc_max_hr is not None:
            opc_times = np.sort(opc_times_all[opc_times_all < float(opc_max_hr)])
        else:
            opc_times = np.sort(opc_times_all)
    else:
        opc_times = np.array([], dtype=float)

    n_opc = int(opc_times.size)

    empty = {
        "alarm": cascade,
        "confirmed": cascade,
        "decision_time_hr": decision_time,
        "sdd_candidate_count": n_sdd,
        "opc_alarm_count": n_opc,
        "cascade_alarm_count": 0,
        "confirm_window_h": float(confirm_window_h),
        "opc_max_hr": opc_max_hr if opc_max_hr is None else float(opc_max_hr),
    }

    if n_sdd == 0 or n_opc == 0:
        return empty

    w = float(confirm_window_h)
    for i in sdd_idx:
        ti = float(t[i])
        if not np.isfinite(ti):
            continue
        # Closed interval: any eligible OPC with |t_opc - ti| <= w
        left = int(np.searchsorted(opc_times, ti - w, side="left"))
        right = int(np.searchsorted(opc_times, ti + w, side="right"))
        if right > left:
            t_opc_star = float(opc_times[left])  # earliest confirming OPC in window
            # Among all in [left, right), earliest is opc_times[left] since sorted
            decision = max(ti, t_opc_star)
            # If opc_max_hr set, decision must also be strictly before horizon
            if opc_max_hr is not None and not (decision < float(opc_max_hr)):
                continue
            cascade[i] = 1
            decision_time[i] = decision

    return {
        "alarm": cascade,
        "confirmed": cascade,
        "decision_time_hr": decision_time,
        "sdd_candidate_count": n_sdd,
        "opc_alarm_count": n_opc,
        "cascade_alarm_count": int(np.sum(cascade)),
        "confirm_window_h": float(confirm_window_h),
        "opc_max_hr": opc_max_hr if opc_max_hr is None else float(opc_max_hr),
    }


def cascade_first_causal_detection(
    sdd_alarm: np.ndarray,
    opc_alarm: np.ndarray,
    t_hr: np.ndarray,
    *,
    confirm_window_h: float = CONFIRM_WINDOW_H,
    search_start_hr: float,
    event_hr: float,
) -> Dict[str, Any]:
    """
    First *causal* cascade detection before an event.

    - Only OPC with t_opc < event_hr may confirm (no post-event look-ahead).
    - SDD candidates considered only if t_sdd is in (search_start_hr, event_hr)
      (same exclusive bounds style as count_binary_alarm_episodes search).
    - decision_time = max(t_sdd, earliest confirming pre-event OPC in ±W).
    - Detection requires decision_time < event_hr and decision_time > search_start_hr.

    Returns dict:
        alarmed (0/1), detection_hr (decision time), lead_time_h,
        sdd_idx, sdd_hr, opc_confirm_hr, cascade_out (full merger dict)
    """
    t = np.asarray(t_hr, dtype=float)
    out = cascade_sdd_confirm_opc(
        sdd_alarm,
        opc_alarm,
        t,
        confirm_window_h=confirm_window_h,
        opc_max_hr=float(event_hr),
    )
    alarm = out["alarm"]
    dec = out["decision_time_hr"]
    t0 = float(search_start_hr)
    t_ev = float(event_hr)

    best_i: Optional[int] = None
    best_dec = float("inf")
    for i in range(len(alarm)):
        if not alarm[i]:
            continue
        ti = float(t[i])
        di = float(dec[i])
        if not np.isfinite(ti) or not np.isfinite(di):
            continue
        # SDD sample must lie in pre-event search (exclusive bounds)
        if not (ti > t0 and ti < t_ev):
            continue
        # Decision must be causal and inside search horizon
        if not (di > t0 and di < t_ev):
            continue
        if di < best_dec:
            best_dec = di
            best_i = i

    if best_i is None:
        return {
            "alarmed": 0,
            "detection_hr": float("nan"),
            "lead_time_h": float("nan"),
            "sdd_idx": None,
            "sdd_hr": float("nan"),
            "opc_confirm_hr": float("nan"),
            "cascade_out": out,
        }

    sdd_hr = float(t[best_i])
    # opc_confirm_hr = decision if OPC after SDD else decision equals sdd when OPC earlier
    # Recover earliest OPC in window from decision definition: decision = max(sdd, opc_star)
    # so opc_star = decision if decision > sdd else <= sdd; report decision - for clarity
    # store opc as decision when OPC is later, else sdd_hr - (decision - sdd) ... simpler:
    opc_confirm_hr = best_dec if best_dec > sdd_hr + 1e-15 else sdd_hr
    # If OPC was before or at SDD, decision == sdd_hr and opc <= sdd; leave opc_confirm as decision
    # when equal. For reporting, recompute earliest eligible OPC in window:
    opc_confirm_hr = _earliest_opc_in_window(
        opc_alarm, t, sdd_hr, confirm_window_h, event_hr
    )

    return {
        "alarmed": 1,
        "detection_hr": best_dec,
        "lead_time_h": float(t_ev - best_dec),
        "sdd_idx": int(best_i),
        "sdd_hr": sdd_hr,
        "opc_confirm_hr": opc_confirm_hr,
        "cascade_out": out,
    }


def cascade_first_alarm_index(
    sdd_alarm: np.ndarray,
    opc_alarm: np.ndarray,
    t_hr: np.ndarray,
    *,
    confirm_window_h: float = CONFIRM_WINDOW_H,
    search_start: int = 0,
    search_end: Optional[int] = None,
    opc_max_hr: Optional[float] = None,
) -> Tuple[Optional[int], Dict[str, Any]]:
    """
    First cascade alarm index in [search_start, search_end).

    Prefer cascade_first_causal_detection for event sensitivity (uses decision time).
    """
    out = cascade_sdd_confirm_opc(
        sdd_alarm,
        opc_alarm,
        t_hr,
        confirm_window_h=confirm_window_h,
        opc_max_hr=opc_max_hr,
    )
    alarm = out["alarm"]
    T = len(alarm)
    end = T if search_end is None else min(int(search_end), T)
    start = max(0, int(search_start))
    for idx in range(start, end):
        if alarm[idx]:
            return idx, out
    return None, out


def _earliest_opc_in_window(
    opc_alarm: np.ndarray,
    t_hr: np.ndarray,
    t_sdd: float,
    confirm_window_h: float,
    event_hr: float,
) -> float:
    opc_on = _as_on(opc_alarm)
    t = np.asarray(t_hr, dtype=float)
    w = float(confirm_window_h)
    best = float("nan")
    for j in np.flatnonzero(opc_on):
        to = float(t[j])
        if not np.isfinite(to) or not (to < float(event_hr)):
            continue
        if abs(to - float(t_sdd)) <= w:
            if not np.isfinite(best) or to < best:
                best = to
    return best


def _as_on(a: np.ndarray) -> np.ndarray:
    """Nonzero finite → True; NaN → False."""
    a = np.asarray(a, dtype=float)
    return np.isfinite(a) & (a != 0)
