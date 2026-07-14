"""Unit tests for structural I0 surplus arms (basal modes, hop, scorer)."""

from __future__ import annotations

import numpy as np
import pytest

from ordinal_detectors.opc_refinements import (
    compute_surplus_threshold,
    ops_detect,
    opsp_integrated_detect,
)
from run_i0_structural_arms import build_structural_arms, score_arm


def test_compute_surplus_threshold_mean_delta():
    vals = np.array([0.10, 0.12, 0.08, 0.11])
    s_ref, thr, mode = compute_surplus_threshold(
        vals, basal_mode="mean_delta", theta_delta_S=0.08
    )
    assert mode == "mean_delta"
    assert abs(s_ref - float(np.mean(vals))) < 1e-12
    assert abs(thr - (s_ref + 0.08)) < 1e-12


def test_compute_surplus_threshold_percentile_above_median():
    vals = np.linspace(0.05, 0.20, 40)
    s_ref, thr, mode = compute_surplus_threshold(
        vals, basal_mode="percentile", basal_q=90.0
    )
    assert mode == "percentile"
    assert thr >= s_ref
    assert thr >= float(np.percentile(vals, 90)) - 1e-12


def test_compute_surplus_threshold_mad_scales():
    vals = np.array([0.10] * 20 + [0.11, 0.09, 0.10, 0.105])
    s_ref, thr, mode = compute_surplus_threshold(
        vals, basal_mode="mad", mad_kappa=2.5
    )
    assert mode == "mad"
    assert thr > s_ref


def test_percentile_basal_quieter_than_mean_delta_on_noisy_basal():
    """High basal tail: q90 thr should suppress post-basal noise better than mean+0.08."""
    rng = np.random.default_rng(0)
    L, k1, k2 = 30, 4, 4
    n_b, n_app = 100, 80
    # Mild dependence throughout (high basal surplus noise)
    pi1 = rng.integers(0, k1, size=n_b + n_app)
    pi2 = pi1.copy()
    flip = rng.random(n_b + n_app) < 0.55
    pi2[flip] = rng.integers(0, k2, size=int(flip.sum()))
    basal_end = n_b

    mean_d = ops_detect(
        pi1, pi2, L=L, theta_R=3, k1=k1, k2=k2,
        basal_end=basal_end, theta_delta_S=0.02, basal_mode="mean_delta", hop=1,
    )
    q90 = ops_detect(
        pi1, pi2, L=L, theta_R=3, k1=k1, k2=k2,
        basal_end=basal_end, basal_mode="percentile", basal_q=90.0, hop=1,
    )
    # q90 threshold is higher → not more alarms than loose mean_delta
    assert int(q90["alarm"][basal_end:].sum()) <= int(mean_d["alarm"][basal_end:].sum())
    assert float(q90["S_thr"][0]) >= float(mean_d["S_thr"][0]) - 1e-9


def test_hop_reduces_alarm_mass_vs_hop1():
    """Semi-independent hop credits fewer pseudo-runs than hop=1."""
    rng = np.random.default_rng(3)
    L, k1, k2 = 20, 4, 4
    n_b, n_app = 60, 120
    pi1 = rng.integers(0, k1, size=n_b + n_app)
    pi2 = rng.integers(0, k2, size=n_b + n_app)
    # short intermittent dependence bursts (autocorrelated under hop=1)
    for start in range(n_b, n_b + n_app, 15):
        end = min(start + 8, n_b + n_app)
        pi2[start:end] = pi1[start:end]
    basal_end = n_b

    h1 = ops_detect(
        pi1, pi2, L=L, theta_R=3, k1=k1, k2=k2,
        basal_end=basal_end, theta_delta_S=0.05, hop=1,
    )
    h10 = ops_detect(
        pi1, pi2, L=L, theta_R=3, k1=k1, k2=k2,
        basal_end=basal_end, theta_delta_S=0.05, hop=10,
    )
    assert int(h10["alarm"][basal_end:].sum()) <= int(h1["alarm"][basal_end:].sum())
    assert int(h10["hop"][0]) == 10


def test_hop_still_alarms_on_long_sustained_dependence():
    rng = np.random.default_rng(5)
    L, k1, k2 = 25, 4, 4
    n_b, n_app = 70, 150
    pi1 = rng.integers(0, k1, size=n_b + n_app)
    pi2 = rng.integers(0, k2, size=n_b + n_app)
    pi2[n_b:] = pi1[n_b:]
    basal_end = n_b
    hop = 12
    out = opsp_integrated_detect(
        pi1, pi2, L=L, theta_R=3, k1=k1, k2=k2,
        basal_end=basal_end, theta_delta_S=0.06,
        collapse_role="none", hop=hop,
    )
    assert int(out["alarm"][basal_end:].sum()) > 0
    assert int(out["core_alarm"][basal_end:].sum()) > 0


def test_legacy_mean_delta_hop1_matches_between_ops_and_integrated():
    rng = np.random.default_rng(11)
    L, theta_R = 30, 4
    k1 = k2 = 4
    n_b, n_app = 60, 80
    pi1 = rng.integers(0, k1, size=n_b + n_app)
    pi2 = rng.integers(0, k2, size=n_b + n_app)
    pi2[n_b:] = pi1[n_b:]
    basal_end = n_b
    ops = ops_detect(
        pi1, pi2, L=L, theta_R=theta_R, k1=k1, k2=k2,
        basal_end=basal_end, theta_delta_S=0.08, hop=1,
    )
    integ = opsp_integrated_detect(
        pi1, pi2, L=L, theta_R=theta_R, k1=k1, k2=k2,
        basal_end=basal_end, theta_delta_S=0.08,
        collapse_role="none", hop=1,
    )
    assert np.array_equal(integ["high_syn"], ops["high_syn"])
    assert np.array_equal(integ["core_alarm"], ops["alarm"])


def test_build_structural_arms_six_fixed():
    arms = build_structural_arms()
    assert len(arms) == 6
    names = [a["name"] for a in arms]
    assert names[0] == "R0_I0_ref"
    assert names[-1] == "R5_combo_Ls100"
    assert arms[1]["hop"] == 25 and arms[1]["theta_R"] == 3
    assert arms[2]["basal_mode"] == "percentile"
    assert arms[3]["basal_mode"] == "mad"
    assert arms[4]["hop"] == 25 and arms[4]["basal_mode"] == "percentile"
    assert arms[5]["L"] == 100 and arms[5]["hop"] == 50
    # no Cartesian explosion
    assert len({(a["L"], a["basal_mode"], a["hop"], a["theta_R"]) for a in arms}) == 6


def test_score_arm_structural_and_clear():
    # FAR just under 2×OPC (~7.47), sens 0.55 → structural only
    s = score_arm(0.55, 7.0)
    assert s["structural_win"] is True
    assert s["clear_advance_original"] is False
    # jackpot
    s2 = score_arm(0.70, 7.0)
    assert s2["clear_advance_original"] is True
    # fail FAR
    s3 = score_arm(0.80, 40.0)
    assert s3["structural_win"] is False
    assert s3["approaches"] is False


def test_percentile_requires_no_theta_delta():
    rng = np.random.default_rng(9)
    pi1 = rng.integers(0, 4, size=200)
    pi2 = pi1.copy()
    out = opsp_integrated_detect(
        pi1, pi2, L=40, theta_R=3, k1=4, k2=4,
        basal_end=80, basal_mode="percentile", basal_q=90.0,
        collapse_role="none",
    )
    assert "S_thr" in out and np.isfinite(out["S_thr"][0])
    assert str(out["basal_mode"][0]) == "percentile"
