"""Tests for Phase β surplus-over-surrogate helpers (P2-10)."""
from __future__ import annotations

import numpy as np
import pytest

from proxy_variants import build_proxy
from run_surplus_proxy_bakeoff import factors_and_surplus
from surplus_surrogates import (
    PROXY_DEFAULT,
    iaaft_1d,
    phase_randomize_1d,
    surplus_from_proxy_matrix,
    surplus_raw_and_excess,
    surrogate_channel2,
)


@pytest.fixture
def rr() -> np.ndarray:
    rng = np.random.default_rng(7)
    t = np.linspace(0, 80, 1600)
    return 820.0 + 35.0 * np.sin(0.3 * t) + rng.normal(0, 15.0, size=t.size)


def test_phase_randomize_shape_and_finite(rr: np.ndarray) -> None:
    s = phase_randomize_1d(rr, seed=1)
    assert s.shape == rr.shape
    assert np.all(np.isfinite(s))
    # not identical to original (almost surely)
    assert not np.allclose(s, rr)


def test_iaaft_preserves_sorted_values(rr: np.ndarray) -> None:
    s = iaaft_1d(rr, n_iters=5, seed=2)
    assert s.shape == rr.shape
    assert np.allclose(np.sort(s), np.sort(rr))


def test_surrogate_channel2_keeps_ch1(rr: np.ndarray) -> None:
    X = build_proxy(rr, "P2-10")
    Xs = surrogate_channel2(X, method="phase", seed=3)
    assert Xs.shape == X.shape
    assert np.allclose(Xs[:, 0], X[:, 0])
    assert not np.allclose(Xs[:, 1], X[:, 1])


def test_surplus_shapes(rr: np.ndarray) -> None:
    X = build_proxy(rr, PROXY_DEFAULT)
    pi1, pi2, S, offset = surplus_from_proxy_matrix(X, L=50)
    assert pi1.shape == pi2.shape
    assert S.shape == pi1.shape
    assert offset == 2  # (m-1)*delay for m=3,delay=1
    assert np.all(np.isfinite(S[np.isfinite(S)]))  # no infs


def test_p2_10_raw_matches_bakeoff_path(rr: np.ndarray) -> None:
    """Raw S from surplus_surrogates must match bake-off factors_and_surplus(P2-10)."""
    _, _, S_ref, off_ref = factors_and_surplus(rr, "P2-10", L=50)
    out = surplus_raw_and_excess(
        rr, proxy_id="P2-10", L=50, n_surr=3, method="phase", base_seed=99
    )
    assert out["offset"] == off_ref
    assert out["S"].shape == S_ref.shape
    assert np.allclose(out["S"], S_ref, equal_nan=True)


def test_excess_shapes_and_relation(rr: np.ndarray) -> None:
    out = surplus_raw_and_excess(
        rr, proxy_id="P2-10", L=50, n_surr=5, method="phase", base_seed=11
    )
    S, S_ex, S_z, med = out["S"], out["S_ex"], out["S_z"], out["S_surr_med"]
    assert S.shape == S_ex.shape == S_z.shape == med.shape
    assert np.allclose(S_ex, S - med, equal_nan=True)
    # With n_surr>=1, std>0 almost everywhere after warm-up
    finite = np.isfinite(S_z)
    assert finite.sum() > 100


def test_unknown_method_raises(rr: np.ndarray) -> None:
    X = build_proxy(rr, "P2-10")
    with pytest.raises(ValueError):
        surrogate_channel2(X, method="not_a_method", seed=0)
