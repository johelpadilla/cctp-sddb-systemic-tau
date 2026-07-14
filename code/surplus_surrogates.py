#!/usr/bin/env python3
"""
Surplus-over-surrogate helpers for Phase β (diagnostic only).

Null model: destroy *cross*-dependence by phase-randomizing (or IAAFT) the
second continuous channel of a bivariate proxy while keeping channel 1 fixed,
then recompute Bandt–Pompe TV surplus under the same (m, L) settings.

Does NOT mutate production `build_bivariate_proxy` or any detector path.

Observational research only — no clinical claims.
"""
from __future__ import annotations

from typing import Any, Dict, Optional, Tuple

import numpy as np

from ordinal_detectors.opc_refinements import ordinal_synergy_surplus
from proxy_variants import build_proxy
from recd_ordinal_levels import generate_multivariate_symbols

# Frozen defaults (match Phase α / lag screen)
M_EMB = 3
DELAY = 1
BP_ALPHABET = 6
L_SURPLUS = 50
N_SURR_DEFAULT = 19  # odd → robust pointwise median
PROXY_DEFAULT = "P2-10"


def phase_randomize_1d(x: np.ndarray, seed: Optional[int] = None) -> np.ndarray:
    """
    Fourier phase randomization: preserve amplitude spectrum, randomize phases.

    Destroys nonlinear / temporal ordinal structure while keeping power spectrum
    (and approximately linear autocovariance).
    """
    rng = np.random.default_rng(seed)
    x = np.asarray(x, dtype=float).ravel()
    n = x.size
    if n < 2:
        return x.copy()
    xf = np.fft.rfft(x)
    amp = np.abs(xf)
    phase = rng.uniform(0.0, 2.0 * np.pi, size=amp.size)
    phase[0] = 0.0
    if amp.size > 1:
        # Nyquist bin for even n: keep real (phase 0 or π only) — set 0
        phase[-1] = 0.0
    s = np.fft.irfft(amp * np.exp(1j * phase), n=n)
    return np.asarray(s, dtype=float)


def iaaft_1d(
    x: np.ndarray,
    n_iters: int = 20,
    seed: Optional[int] = None,
) -> np.ndarray:
    """
    Iterative Amplitude Adjusted Fourier Transform (IAAFT) 1D surrogate.

    Preserves the amplitude distribution of x and approximately the power
    spectrum. Adapted from run_cctp_surrogates.py (local copy for isolation).
    """
    rng = np.random.default_rng(seed)
    x = np.asarray(x, dtype=float).ravel()
    n = x.size
    if n < 2:
        return x.copy()
    xf = np.fft.rfft(x)
    amp = np.abs(xf)
    phase = rng.uniform(0.0, 2.0 * np.pi, size=amp.size)
    phase[0] = 0.0
    if amp.size > 1:
        phase[-1] = 0.0
    s = np.fft.irfft(amp * np.exp(1j * phase), n=n)
    x_sorted = np.sort(x)
    for _ in range(int(n_iters)):
        sf = np.fft.rfft(s)
        s = np.fft.irfft(amp * np.exp(1j * np.angle(sf)), n=n)
        ranks = np.argsort(np.argsort(s))
        s = x_sorted[ranks]
    return np.asarray(s, dtype=float)


def surrogate_channel2(
    X: np.ndarray,
    *,
    method: str = "phase",
    seed: Optional[int] = None,
    iaaft_iters: int = 20,
) -> np.ndarray:
    """
    Return a copy of bivariate proxy X with channel 2 replaced by a surrogate.

    method:
      'phase' — Fourier phase randomization (default; fast)
      'iaaft' — IAAFT (slower; matches amplitude histogram)
    """
    X = np.asarray(X, dtype=float)
    if X.ndim != 2 or X.shape[1] != 2:
        raise ValueError(f"X must be (n,2); got {X.shape}")
    method = str(method).lower().strip()
    Xs = X.copy()
    if method in ("phase", "phase_random", "fourier"):
        Xs[:, 1] = phase_randomize_1d(X[:, 1], seed=seed)
    elif method in ("iaaft", "iAAFT"):
        Xs[:, 1] = iaaft_1d(X[:, 1], n_iters=iaaft_iters, seed=seed)
    else:
        raise ValueError(f"Unknown surrogate method {method!r}; use 'phase' or 'iaaft'")
    return Xs


def surplus_from_proxy_matrix(
    X: np.ndarray,
    *,
    L: int = L_SURPLUS,
    m: int = M_EMB,
    delay: int = DELAY,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, int]:
    """Bandt–Pompe factors + TV synergy surplus from continuous proxy (n,2)."""
    X = np.asarray(X, dtype=float)
    S_sym = generate_multivariate_symbols(X, m=m, delay=delay)
    offset = (m - 1) * delay
    if S_sym.size == 0 or S_sym.ndim < 2 or S_sym.shape[1] < 2:
        empty_i = np.array([], dtype=np.int64)
        empty_f = np.array([], dtype=float)
        return empty_i, empty_i, empty_f, offset
    pi1 = S_sym[:, 0].astype(np.int64)
    pi2 = S_sym[:, 1].astype(np.int64)
    S = ordinal_synergy_surplus(pi1, pi2, L=L, k1=BP_ALPHABET, k2=BP_ALPHABET)
    return pi1, pi2, S, offset


def surplus_raw_and_excess(
    rr: np.ndarray,
    *,
    proxy_id: str = PROXY_DEFAULT,
    L: int = L_SURPLUS,
    n_surr: int = N_SURR_DEFAULT,
    method: str = "phase",
    base_seed: int = 20260714,
    iaaft_iters: int = 20,
) -> Dict[str, Any]:
    """
    Compute raw S on proxy and pointwise surplus-over-surrogate.

    Returns dict with:
      S          — raw surplus series
      S_ex       — S - median_k(S_surr_k)   (raw excess)
      S_z        — (S - mean_k) / std_k      (per-time z over surrogates)
      S_surr_med — pointwise median of surrogates
      S_surr_mean, S_surr_std
      pi1, pi2, offset
      n_surr, method, proxy_id
    """
    if n_surr < 1:
        raise ValueError("n_surr must be >= 1")
    rr = np.asarray(rr, dtype=float).ravel()
    X = build_proxy(rr, proxy_id)
    pi1, pi2, S, offset = surplus_from_proxy_matrix(X, L=L)
    T = int(S.size)
    if T == 0:
        empty = np.array([], dtype=float)
        return {
            "proxy_id": proxy_id,
            "method": method,
            "n_surr": int(n_surr),
            "L": int(L),
            "offset": int(offset),
            "S": empty,
            "S_ex": empty,
            "S_z": empty,
            "S_surr_med": empty,
            "S_surr_mean": empty,
            "S_surr_std": empty,
            "pi1": pi1,
            "pi2": pi2,
            "X": X,
        }

    stack = np.empty((n_surr, T), dtype=float)
    for k in range(n_surr):
        seed = int(base_seed) + 10007 * k
        Xs = surrogate_channel2(X, method=method, seed=seed, iaaft_iters=iaaft_iters)
        _, _, Sk, _ = surplus_from_proxy_matrix(Xs, L=L)
        # Align lengths if embedding edge cases differ (should not)
        n = min(T, Sk.size)
        row = np.full(T, np.nan, dtype=float)
        row[:n] = Sk[:n]
        stack[k] = row

    # Leading samples of TV surplus are NaN (window warm-up); suppress all-NaN warns
    with np.errstate(all="ignore"):
        surr_med = np.nanmedian(stack, axis=0)
        surr_mean = np.nanmean(stack, axis=0)
        surr_std = np.nanstd(stack, axis=0, ddof=0)
    S_ex = S - surr_med
    S_z = (S - surr_mean) / (surr_std + 1e-12)

    return {
        "proxy_id": proxy_id,
        "method": method,
        "n_surr": int(n_surr),
        "L": int(L),
        "offset": int(offset),
        "S": S,
        "S_ex": S_ex,
        "S_z": S_z,
        "S_surr_med": surr_med,
        "S_surr_mean": surr_mean,
        "S_surr_std": surr_std,
        "pi1": pi1,
        "pi2": pi2,
        "X": X,
        "S_surr_stack_mean_over_t": float(np.nanmean(stack)),
    }
