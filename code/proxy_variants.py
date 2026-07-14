#!/usr/bin/env python3
"""
Bivariate proxy variants for Phase-α surplus bake-off (diagnostic only).

Does NOT mutate `build_bivariate_proxy` in cctp_metrics_core.py — production
τ_s / OPC / shipped detectors keep the frozen R0 path.

Phase-α frozen set (Bandt–Pompe m=3, L=50, TV surplus shared elsewhere):

  R0  [z(RR), z(|ΔRR|)]                         — current baseline
  P1  [z(RR), z(e_AR1)]                         — global AR(1) residual of RR
  P2  [z(RR_t), z(|ΔRR|_{t-ℓ})]  ℓ=10           — lagged absolute increment
  P3  [z(RR), z(ΔRR)]                           — signed increment (control)
  P4  [z(RR), z(r)]  r = |ΔRR| − OLS(|ΔRR|~RR) — residualized |ΔRR| on RR

Observational research only — no clinical claims.
"""
from __future__ import annotations

from typing import Callable, Dict, List, Sequence, Tuple

import numpy as np

from cctp_metrics_core import build_bivariate_proxy

# Frozen hyperparameters for Phase α (chosen a priori; do not retune on events)
P2_LAG = 10
P1_MIN_N = 20  # minimum samples to fit AR(1)

# Phase-β mini screen: discrete lags only (not a continuous optimization grid)
P2_LAG_SCREEN: Tuple[int, ...] = (5, 10, 20)


def _zscore(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    return (x - np.mean(x)) / (np.std(x) + 1e-12)


def _abs_drr(rr: np.ndarray) -> np.ndarray:
    rr = np.asarray(rr, dtype=float)
    return np.abs(np.diff(rr, prepend=rr[0]))


def _signed_drr(rr: np.ndarray) -> np.ndarray:
    rr = np.asarray(rr, dtype=float)
    return np.diff(rr, prepend=rr[0])


def build_proxy_r0(rr: np.ndarray) -> np.ndarray:
    """R0 — current production bivariate proxy (bitwise via core helper)."""
    return build_bivariate_proxy(np.asarray(rr, dtype=float))


def build_proxy_p1_ar1_residual(rr: np.ndarray) -> np.ndarray:
    """
    P1 — [z(RR), z(e)] where e is the global AR(1) residual of RR.

    Fit RR_t = c + φ RR_{t-1} + e_t by OLS on the full series (P1a).
    e_0 = 0. If the design matrix is degenerate, fall back to first difference.
    """
    rr = np.asarray(rr, dtype=float).ravel()
    n = rr.size
    e = np.zeros(n, dtype=float)
    if n < P1_MIN_N:
        e[1:] = rr[1:] - rr[:-1]
    else:
        y = rr[1:]
        x = rr[:-1]
        # OLS: y = c + φ x
        X = np.column_stack([np.ones(n - 1), x])
        try:
            beta, *_ = np.linalg.lstsq(X, y, rcond=None)
            e[1:] = y - X @ beta
            e[0] = 0.0
        except np.linalg.LinAlgError:
            e[1:] = y - x
            e[0] = 0.0
    return np.column_stack([_zscore(rr), _zscore(e)])


def build_proxy_p2_lagged(
    rr: np.ndarray,
    *,
    lag: int = P2_LAG,
) -> np.ndarray:
    """
    P2 — [z(RR_t), z(|ΔRR|_{t-lag})].

    lag is frozen at P2_LAG=10 for Phase α (screening arm).
    Leading samples use the earliest available |ΔRR| (causal pad).
    """
    rr = np.asarray(rr, dtype=float).ravel()
    if lag < 1:
        raise ValueError("lag must be >= 1")
    drr = _abs_drr(rr)
    drr_lag = np.empty_like(drr)
    drr_lag[:lag] = drr[0]
    drr_lag[lag:] = drr[:-lag]
    return np.column_stack([_zscore(rr), _zscore(drr_lag)])


def build_proxy_p3_signed(rr: np.ndarray) -> np.ndarray:
    """P3 — [z(RR), z(ΔRR)] signed increment (no absolute value)."""
    rr = np.asarray(rr, dtype=float).ravel()
    drr = _signed_drr(rr)
    return np.column_stack([_zscore(rr), _zscore(drr)])


def build_proxy_p4_drr_residual(rr: np.ndarray) -> np.ndarray:
    """
    P4 — [z(RR), z(r)] where r is the residual of |ΔRR| after OLS on RR.

    Fit |ΔRR|_t = a + b RR_t + r_t (global). Directly attacks structural
    coupling of the R0 pair (redesign P1c).
    """
    rr = np.asarray(rr, dtype=float).ravel()
    drr = _abs_drr(rr)
    n = rr.size
    if n < 10:
        r = drr.copy()
    else:
        X = np.column_stack([np.ones(n), rr])
        try:
            beta, *_ = np.linalg.lstsq(X, drr, rcond=None)
            r = drr - X @ beta
        except np.linalg.LinAlgError:
            r = drr - np.mean(drr)
    return np.column_stack([_zscore(rr), _zscore(r)])


# Canonical Phase-α registry (order fixed)
PROXY_BUILDERS: Dict[str, Callable[[np.ndarray], np.ndarray]] = {
    "R0": build_proxy_r0,
    "P1": build_proxy_p1_ar1_residual,
    "P2": build_proxy_p2_lagged,
    "P3": build_proxy_p3_signed,
    "P4": build_proxy_p4_drr_residual,
}

PROXY_DESCRIPTIONS: Dict[str, str] = {
    "R0": "[z(RR), z(|ΔRR|)] current baseline",
    "P1": f"[z(RR), z(e_AR1)] global AR(1) residual of RR",
    "P2": f"[z(RR_t), z(|ΔRR|_{{t-{P2_LAG}}})] lagged absolute increment",
    "P3": "[z(RR), z(ΔRR)] signed increment (control arm)",
    "P4": "[z(RR), z(r)] residual of |ΔRR| after OLS on RR",
}

PROXY_IDS: Tuple[str, ...] = ("R0", "P1", "P2", "P3", "P4")


def _p2_lag_id(lag: int) -> str:
    return f"P2-{int(lag)}"


def _make_p2_lag_builder(lag: int) -> Callable[[np.ndarray], np.ndarray]:
    def _builder(rr: np.ndarray) -> np.ndarray:
        return build_proxy_p2_lagged(rr, lag=lag)

    _builder.__name__ = f"build_proxy_p2_lag{lag}"
    _builder.__doc__ = f"P2 with fixed lag ℓ={lag}."
    return _builder


# Discrete lag screen (Phase β mini): P2-5, P2-10, P2-20
P2_LAG_BUILDERS: Dict[str, Callable[[np.ndarray], np.ndarray]] = {
    _p2_lag_id(lag): _make_p2_lag_builder(lag) for lag in P2_LAG_SCREEN
}
P2_LAG_DESCRIPTIONS: Dict[str, str] = {
    _p2_lag_id(lag): (
        f"[z(RR_t), z(|ΔRR|_{{t-{lag}}})] lagged |ΔRR| (ℓ={lag})"
    )
    for lag in P2_LAG_SCREEN
}
P2_LAG_IDS: Tuple[str, ...] = tuple(_p2_lag_id(lag) for lag in P2_LAG_SCREEN)

# Unified dispatch (Phase-α + lag screen). Production path still uses R0 only.
ALL_PROXY_BUILDERS: Dict[str, Callable[[np.ndarray], np.ndarray]] = {
    **PROXY_BUILDERS,
    **P2_LAG_BUILDERS,
}
ALL_PROXY_DESCRIPTIONS: Dict[str, str] = {
    **PROXY_DESCRIPTIONS,
    **P2_LAG_DESCRIPTIONS,
}


def build_proxy(rr: np.ndarray, kind: str) -> np.ndarray:
    """Dispatch by proxy id. Raises KeyError for unknown kind."""
    kind = str(kind).upper()
    if kind not in ALL_PROXY_BUILDERS:
        known = list(PROXY_IDS) + list(P2_LAG_IDS)
        raise KeyError(f"Unknown proxy kind {kind!r}; choose from {known}")
    return ALL_PROXY_BUILDERS[kind](np.asarray(rr, dtype=float))


def list_proxies() -> List[Dict[str, str]]:
    return [
        {"id": pid, "description": PROXY_DESCRIPTIONS[pid]}
        for pid in PROXY_IDS
    ]


def list_p2_lag_proxies() -> List[Dict[str, str]]:
    return [
        {"id": pid, "description": P2_LAG_DESCRIPTIONS[pid]}
        for pid in P2_LAG_IDS
    ]


def assert_proxy_shape(X: np.ndarray, n: int) -> None:
    X = np.asarray(X)
    if X.ndim != 2 or X.shape != (n, 2):
        raise ValueError(f"proxy must be shape (n,2); got {X.shape}")
    if not np.all(np.isfinite(X)):
        raise ValueError("proxy contains non-finite values")
