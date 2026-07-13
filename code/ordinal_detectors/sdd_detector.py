"""
Option 2 — Symbolic Distribution Divergence (SDD) detector.

Pure reference implementation on a discrete ordinal symbol stream.
Uses Total Variation between basal and current empirical distributions.
Does NOT use mean/variance of any continuous series.
Does NOT call or wrap abs-z / detect_lead_time.
Does NOT fuse with Ordinal Persistence Collapse (Option 1).

See docs/ORDINAL_ALARM_DETECTORS.md §2.
"""

from __future__ import annotations

from typing import Dict, Optional, Sequence, Tuple

import numpy as np


def _validate_symbols(sigma: np.ndarray, K: int) -> np.ndarray:
    sigma = np.asarray(sigma)
    if sigma.ndim != 1:
        raise ValueError("sigma must be a 1-D integer symbol stream")
    if K < 2:
        raise ValueError("K (alphabet size) must be >= 2")
    if sigma.size and (np.min(sigma) < 0 or np.max(sigma) >= K):
        raise ValueError(f"symbols must lie in {{0,...,{K-1}}}")
    return sigma.astype(np.int64, copy=False)


def empirical_distribution(
    sigma: np.ndarray,
    start: int,
    end: int,
    K: int,
) -> np.ndarray:
    """
    Empirical distribution on symbols sigma[start:end] (end exclusive).
    Returns probability vector of length K.
    """
    if end <= start:
        raise ValueError("empty segment for empirical distribution")
    seg = sigma[start:end]
    counts = np.bincount(seg, minlength=K).astype(float)
    return counts / counts.sum()


def total_variation(p: np.ndarray, q: np.ndarray) -> float:
    """TV(p,q) = 0.5 * ||p-q||_1 on the probability simplex."""
    p = np.asarray(p, dtype=float)
    q = np.asarray(q, dtype=float)
    if p.shape != q.shape:
        raise ValueError("p and q must have the same shape")
    return float(0.5 * np.sum(np.abs(p - q)))


def kl_divergence_smoothed(
    p: np.ndarray,
    q: np.ndarray,
    eps: float = 1e-9,
) -> float:
    """
    Secondary diagnostic only: smoothed KL(p || q).
    Not used by the recommended SDD alarm rule.
    """
    p = np.asarray(p, dtype=float)
    q = np.asarray(q, dtype=float)
    K = p.size
    p_s = (p + eps) / (1.0 + K * eps)
    q_s = (q + eps) / (1.0 + K * eps)
    return float(np.sum(p_s * np.log(p_s / q_s)))


def sdd_detect(
    sigma: np.ndarray,
    basal: Tuple[int, int],
    *,
    L_c: int = 50,
    theta_TV: float = 0.35,
    theta_S: int = 1,
    K: Optional[int] = None,
    mask_basal: bool = True,
) -> Dict[str, np.ndarray]:
    """
    Symbolic Distribution Divergence (TV) on symbol stream sigma.

    Parameters
    ----------
    sigma : array of int
        Ordinal symbols in {0,...,K-1}.
    basal : (start, end)
        Inclusive-exclusive basal index range [start, end) used for P_basal.
    L_c : int
        Current window length.
    theta_TV : float
        Alarm if TV(P_t, P_basal) >= theta_TV.
    theta_S : int
        Optional sustainment (consecutive high-TV endpoints); default 1.
    K : int, optional
        Alphabet size.
    mask_basal : bool
        If True, do not alarm for t whose current window overlaps basal.

    Returns
    -------
    dict with keys:
        alarm : (T,) int {0,1}
        tv : (T,) float NaN where undefined
        p_basal : (K,) float
        high_tv : (T,) int indicator TV >= theta_TV
    """
    if L_c < 2:
        raise ValueError("L_c must be >= 2")
    if not (0.0 < theta_TV <= 1.0):
        raise ValueError("theta_TV must be in (0, 1]")
    if theta_S < 1:
        raise ValueError("theta_S must be >= 1")

    b0, b1 = int(basal[0]), int(basal[1])
    if b1 <= b0:
        raise ValueError("basal range must be non-empty [start, end)")

    if K is None:
        K = int(np.max(np.asarray(sigma))) + 1 if len(np.asarray(sigma)) else 2
        K = max(K, 2)

    sigma = _validate_symbols(sigma, K)
    T = sigma.size
    if b0 < 0 or b1 > T:
        raise ValueError("basal range out of bounds for sigma")

    p_basal = empirical_distribution(sigma, b0, b1, K)

    alarm = np.zeros(T, dtype=np.int8)
    tv = np.full(T, np.nan, dtype=float)
    high_tv = np.zeros(T, dtype=np.int8)

    run = 0
    for t in range(T):
        if t < L_c - 1:
            continue
        w0 = t - L_c + 1
        if mask_basal and not (w0 >= b1 or t < b0):
            # window overlaps basal: leave undefined / no alarm
            run = 0
            continue
        p_cur = empirical_distribution(sigma, w0, t + 1, K)
        d = total_variation(p_cur, p_basal)
        tv[t] = d
        if d >= theta_TV:
            high_tv[t] = 1
            run += 1
        else:
            run = 0
        if run >= theta_S:
            alarm[t] = 1

    return {
        "alarm": alarm,
        "tv": tv,
        "p_basal": p_basal,
        "high_tv": high_tv,
    }


def sdd_alarm_at(
    sigma: np.ndarray,
    t: int,
    basal: Tuple[int, int],
    *,
    L_c: int = 50,
    theta_TV: float = 0.35,
    theta_S: int = 1,
    K: Optional[int] = None,
    mask_basal: bool = True,
) -> bool:
    """Boolean SDD alarm at index t."""
    out = sdd_detect(
        sigma,
        basal,
        L_c=L_c,
        theta_TV=theta_TV,
        theta_S=theta_S,
        K=K,
        mask_basal=mask_basal,
    )
    if t < 0 or t >= len(out["alarm"]):
        return False
    return bool(out["alarm"][t])


def sdd_first_alarm_index(
    sigma: np.ndarray,
    basal: Tuple[int, int],
    *,
    L_c: int = 50,
    theta_TV: float = 0.35,
    theta_S: int = 1,
    K: Optional[int] = None,
    search_start: int = 0,
    search_end: Optional[int] = None,
    mask_basal: bool = True,
) -> Tuple[Optional[int], Dict[str, np.ndarray]]:
    """First index with A_SDD=1 in [search_start, search_end)."""
    out = sdd_detect(
        sigma,
        basal,
        L_c=L_c,
        theta_TV=theta_TV,
        theta_S=theta_S,
        K=K,
        mask_basal=mask_basal,
    )
    T = len(out["alarm"])
    end = T if search_end is None else min(search_end, T)
    start = max(0, search_start)
    for t in range(start, end):
        if out["alarm"][t]:
            return t, out
    return None, out
