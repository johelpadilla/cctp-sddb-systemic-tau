"""Smoke tests for Phase-α proxy variants (shapes + R0 identity)."""
from __future__ import annotations

import numpy as np
import pytest

from cctp_metrics_core import build_bivariate_proxy
from proxy_variants import (
    P2_LAG,
    P2_LAG_IDS,
    PROXY_IDS,
    build_proxy,
    build_proxy_p1_ar1_residual,
    build_proxy_p2_lagged,
    build_proxy_p3_signed,
    build_proxy_p4_drr_residual,
)


@pytest.fixture
def rr() -> np.ndarray:
    rng = np.random.default_rng(42)
    t = np.linspace(0, 40, 800)
    return 800.0 + 40.0 * np.sin(t) + rng.normal(0, 12.0, size=t.size)


def test_r0_matches_core(rr: np.ndarray) -> None:
    assert np.allclose(build_proxy(rr, "R0"), build_bivariate_proxy(rr))


@pytest.mark.parametrize("pid", list(PROXY_IDS))
def test_proxy_shape_and_finite(rr: np.ndarray, pid: str) -> None:
    X = build_proxy(rr, pid)
    assert X.shape == (rr.size, 2)
    assert np.all(np.isfinite(X))
    # global z-score → unit std per channel (within float noise)
    assert np.allclose(np.std(X, axis=0), 1.0, atol=1e-6)


def test_p2_lag_changes_channel2(rr: np.ndarray) -> None:
    X0 = build_proxy(rr, "R0")
    X2 = build_proxy_p2_lagged(rr, lag=P2_LAG)
    # lag should alter second channel vs contemporaneous |ΔRR|
    assert not np.allclose(X0[:, 1], X2[:, 1])


@pytest.mark.parametrize("pid", list(P2_LAG_IDS))
def test_p2_lag_screen_ids(rr: np.ndarray, pid: str) -> None:
    X = build_proxy(rr, pid)
    assert X.shape == (rr.size, 2)
    assert np.all(np.isfinite(X))


def test_p2_lags_differ(rr: np.ndarray) -> None:
    X5 = build_proxy(rr, "P2-5")
    X10 = build_proxy(rr, "P2-10")
    X20 = build_proxy(rr, "P2-20")
    # different lags → different second channel (not identical)
    assert not np.allclose(X5[:, 1], X10[:, 1])
    assert not np.allclose(X10[:, 1], X20[:, 1])
    # first channel is always z(RR) → same across lags
    assert np.allclose(X5[:, 0], X10[:, 0])
    assert np.allclose(X10[:, 0], X20[:, 0])
    # P2 and P2-10 must match (same ℓ=10)
    assert np.allclose(build_proxy(rr, "P2"), X10)


def test_p3_signed_differs_from_r0(rr: np.ndarray) -> None:
    X0 = build_proxy(rr, "R0")
    X3 = build_proxy_p3_signed(rr)
    assert not np.allclose(X0[:, 1], X3[:, 1])


def test_p1_residual_mean_near_zero(rr: np.ndarray) -> None:
    X = build_proxy_p1_ar1_residual(rr)
    # residual before z-score is mean≈0; after z-score still mean≈0
    assert abs(float(np.mean(X[:, 1]))) < 1e-9


def test_p4_residual_uncorrelated_linearly(rr: np.ndarray) -> None:
    """OLS residual of |ΔRR| on RR should be near-orthogonal to RR (Pearson)."""
    X = build_proxy_p4_drr_residual(rr)
    # after z-score, linear corr of residual channel with RR channel ≈ 0
    r = float(np.corrcoef(X[:, 0], X[:, 1])[0, 1])
    assert abs(r) < 0.05


def test_unknown_proxy_raises(rr: np.ndarray) -> None:
    with pytest.raises(KeyError):
        build_proxy(rr, "P99")
