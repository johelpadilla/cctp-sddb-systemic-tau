"""
Option 1 — Ordinal Persistence Collapse (OPC) detector.

Pure reference implementation on a discrete ordinal symbol stream.
Does NOT use mean/variance of any continuous series.
Does NOT call or wrap abs-z / detect_lead_time.
Does NOT fuse with Symbolic Distribution Divergence (Option 2).

See docs/ORDINAL_ALARM_DETECTORS.md §1.
"""

from __future__ import annotations

from typing import Dict, Optional, Tuple

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


def opc_detect(
    sigma: np.ndarray,
    *,
    L: int = 8,
    theta_D: float = 0.35,
    theta_R: int = 5,
    K: Optional[int] = None,
) -> Dict[str, np.ndarray]:
    """
    Ordinal Persistence Collapse on symbol stream sigma.

    Parameters
    ----------
    sigma : array of int
        Ordinal symbols in {0,...,K-1}.
    L : int
        Observation window length (symbols).
    theta_D : float
        Diversity collapse threshold in (0, 1].
    theta_R : int
        Minimum consecutive low-diversity windows for alarm.
    K : int, optional
        Alphabet size. If None, uses max(sigma)+1 (must be >= 2).

    Returns
    -------
    dict with keys:
        alarm : (T,) int {0,1}
        diversity : (T,) float, NaN where undefined
        persistence : (T,) int run length of low diversity
        low_div : (T,) int indicator D_t <= theta_D
    """
    if L < 2:
        raise ValueError("L must be >= 2")
    if not (0.0 < theta_D <= 1.0):
        raise ValueError("theta_D must be in (0, 1]")
    if theta_R < 1:
        raise ValueError("theta_R must be >= 1")

    if K is None:
        if sigma is None or len(np.asarray(sigma)) == 0:
            K = 2
        else:
            K = int(np.max(np.asarray(sigma))) + 1
            K = max(K, 2)

    sigma = _validate_symbols(sigma, K)
    T = sigma.size
    alarm = np.zeros(T, dtype=np.int8)
    diversity = np.full(T, np.nan, dtype=float)
    persistence = np.zeros(T, dtype=np.int64)
    low_div = np.zeros(T, dtype=np.int8)

    run = 0
    for t in range(T):
        if t < L - 1:
            continue
        window = sigma[t - L + 1 : t + 1]
        n_unique = len(np.unique(window))
        D = n_unique / float(K)
        diversity[t] = D
        if D <= theta_D:
            low_div[t] = 1
            run += 1
        else:
            run = 0
        persistence[t] = run
        if low_div[t] == 1 and run >= theta_R:
            alarm[t] = 1

    return {
        "alarm": alarm,
        "diversity": diversity,
        "persistence": persistence,
        "low_div": low_div,
    }


def opc_alarm_at(
    sigma: np.ndarray,
    t: int,
    *,
    L: int = 8,
    theta_D: float = 0.35,
    theta_R: int = 5,
    K: Optional[int] = None,
) -> bool:
    """Boolean OPC alarm at a single index t (drives full stream for correctness)."""
    out = opc_detect(sigma, L=L, theta_D=theta_D, theta_R=theta_R, K=K)
    if t < 0 or t >= len(out["alarm"]):
        return False
    return bool(out["alarm"][t])


def opc_first_alarm_index(
    sigma: np.ndarray,
    *,
    L: int = 8,
    theta_D: float = 0.35,
    theta_R: int = 5,
    K: Optional[int] = None,
    search_start: int = 0,
    search_end: Optional[int] = None,
) -> Tuple[Optional[int], Dict[str, np.ndarray]]:
    """First index with A_OPC=1 in [search_start, search_end)."""
    out = opc_detect(sigma, L=L, theta_D=theta_D, theta_R=theta_R, K=K)
    T = len(out["alarm"])
    end = T if search_end is None else min(search_end, T)
    start = max(0, search_start)
    for t in range(start, end):
        if out["alarm"][t]:
            return t, out
    return None, out
