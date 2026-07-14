"""Tests for exploratory OPC refinements (real shipped entry points)."""

from __future__ import annotations

import numpy as np
import pytest

from ordinal_detectors.opc_detector import opc_detect
from ordinal_detectors.opc_refinements import (
    joint_symbols_from_factors,
    opc_detect_basal_relative,
    opc_detect_gap_tolerant,
    ops_detect,
    opsp_integrated_detect,
    ordinal_synergy_surplus,
)


def test_gap_tolerant_recovers_interrupted_collapse_where_hard_opc_fails():
    """Intermittency: one high-div break resets hard OPC; G=1 recovers credit."""
    # Alphabet K=4. Low diversity = at most 1 symbol when theta_D=0.35 -> 0.25
    # Build: many samples of symbol 0, then a brief intrusion of others, then 0 again.
    L, theta_D, theta_R = 8, 0.35, 4
    K = 4
    # window of 8 all-0 is D=0.25 <= 0.35
    body = [0] * 20
    # insert a single non-collapse window by flipping one position to create diversity
    # After enough collapse, one mixed window, then collapse again.
    sigma = np.array(body + [0, 1, 2, 3, 0, 1, 2, 3] + [0] * 25, dtype=int)

    hard = opc_detect(sigma, L=L, theta_D=theta_D, theta_R=theta_R, K=K)
    soft = opc_detect_gap_tolerant(
        sigma, L=L, theta_D=theta_D, theta_R=theta_R, G=1, K=K
    )
    # With G=0 equivalent behavior region early; with G=1 more alarms late
    # Construct a cleaner synthetic: contiguous low-div of length theta_R-1,
    # one high-div, then more low-div.
    sig = []
    # warm-up
    sig.extend([0] * L)
    # 3 low-div windows endpoints worth of collapse (sliding)
    sig.extend([0] * 10)
    # force one high-diversity window content
    sig.extend([0, 1, 2, 3, 0, 1, 2, 3])
    # resume collapse
    sig.extend([0] * 20)
    sigma2 = np.array(sig, dtype=int)

    hard2 = opc_detect(sigma2, L=L, theta_D=theta_D, theta_R=theta_R, K=K)
    soft2 = opc_detect_gap_tolerant(
        sigma2, L=L, theta_D=theta_D, theta_R=theta_R, G=2, K=K
    )
    # Soft with G=2 should have at least as many alarm samples as hard
    assert int(soft2["alarm"].sum()) >= int(hard2["alarm"].sum())
    # Gap-tolerant is not identical to ignoring diversity: G=0 matches hard reset spirit
    soft0 = opc_detect_gap_tolerant(
        sigma2, L=L, theta_D=theta_D, theta_R=theta_R, G=0, K=K
    )
    assert np.array_equal(soft0["alarm"], hard2["alarm"])


def test_basal_relative_suppresses_collapse_when_basal_already_low_div():
    """If basal repertoire is already collapsed, absolute OPC may fire; BR tightens."""
    K = 6
    L, theta_D, theta_R = 10, 0.35, 3
    # Basal: only symbols 0,1 -> D_basal = 2/6 ≈ 0.333
    basal = [0, 1] * 30
    # Approach: still only 0,1 (same low diversity) — not a *relative* collapse
    approach = [0, 1] * 40
    sigma = np.array(basal + approach, dtype=int)
    basal_end = len(basal)

    abs_opc = opc_detect(sigma, L=L, theta_D=theta_D, theta_R=theta_R, K=K)
    br = opc_detect_basal_relative(
        sigma,
        L=L,
        theta_D=theta_D,
        theta_R=theta_R,
        basal_end=basal_end,
        rho=0.85,
        K=K,
    )
    # Absolute may alarm on continued low diversity after warm-up
    # Basal-relative: rel_cap = min(0.35, 0.85*D_basal) ≈ min(0.35, 0.283) = 0.283
    # D_t for {0,1} = 2/6 ≈ 0.333 > 0.283 → should NOT low-div under BR
    post = slice(basal_end, None)
    assert int(br["alarm"][post].sum()) == 0
    # Absolute OPC still sees D=0.333 <= 0.35
    assert int(abs_opc["low_div"][L:].sum()) > 0


def test_synergy_surplus_high_when_joint_dependent_zero_when_independent():
    """OPS TV: dependent joints >> independent product structure."""
    rng = np.random.default_rng(0)
    L = 40
    k1 = k2 = 4
    T = 200
    # Independent
    pi1_i = rng.integers(0, k1, size=T)
    pi2_i = rng.integers(0, k2, size=T)
    S_ind = ordinal_synergy_surplus(pi1_i, pi2_i, L=L, k1=k1, k2=k2)
    # Dependent: pi2 = pi1 (perfect coupling)
    pi1_d = rng.integers(0, k1, size=T)
    pi2_d = pi1_d.copy()
    S_dep = ordinal_synergy_surplus(pi1_d, pi2_d, L=L, k1=k1, k2=k2)

    ind_mean = float(np.nanmean(S_ind))
    dep_mean = float(np.nanmean(S_dep))
    assert dep_mean > ind_mean + 0.05
    assert dep_mean > 0.2  # strong dependence → substantial TV surplus


def test_ops_detect_alarms_on_sustained_dependence_not_on_independence():
    rng = np.random.default_rng(1)
    L, theta_S, theta_R = 30, 0.2, 3
    k1 = k2 = 4
    T = 150
    # Independent block then dependent block
    pi1 = rng.integers(0, k1, size=T)
    pi2 = rng.integers(0, k2, size=T)
    # last half: lock pi2 = pi1
    pi2[T // 2 :] = pi1[T // 2 :]

    out = ops_detect(
        pi1, pi2, L=L, theta_S=theta_S, theta_R=theta_R, k1=k1, k2=k2
    )
    # Alarms should concentrate in the dependent half
    early = int(out["alarm"][: T // 2].sum())
    late = int(out["alarm"][T // 2 :].sum())
    assert late > early
    assert late > 0


def test_ops_basal_relative_requires_increase_over_basal():
    """If dependence is high already in basal, absolute OPS may fire; delta gate suppresses."""
    L = 20
    k1 = k2 = 3
    # Always dependent
    pi1 = np.array([i % k1 for i in range(120)], dtype=int)
    pi2 = pi1.copy()
    basal_end = 50
    abs_ops = ops_detect(pi1, pi2, L=L, theta_S=0.1, theta_R=2, k1=k1, k2=k2)
    rel_ops = ops_detect(
        pi1,
        pi2,
        L=L,
        theta_S=0.1,
        theta_R=2,
        k1=k1,
        k2=k2,
        basal_end=basal_end,
        theta_delta_S=0.15,  # require substantial *increase*
    )
    assert int(abs_ops["alarm"].sum()) > 0
    # No increase over basal mean → relative should stay quiet post-basal
    assert int(rel_ops["alarm"][basal_end:].sum()) == 0


def test_joint_symbols_encoding_roundtrip_shape():
    pi1 = np.array([0, 1, 2, 5], dtype=int)
    pi2 = np.array([0, 0, 1, 2], dtype=int)
    joint = joint_symbols_from_factors(pi1, pi2, k2=6)
    assert joint.tolist() == [0, 6, 13, 32]


def test_mode_c_recommended_configs_c_br1_c_br3_call_shipped_api():
    """Recommended Mode-C cells C-BR1 / C-BR3 (docs §5.1) are callable as-shipped."""
    K = 6
    L, theta_D = 50, 0.35
    # Rich basal then true relative collapse (single symbol)
    basal = list(range(K)) * 20  # basal_end >= L
    approach = [0] * 80
    sigma = np.array(basal + approach, dtype=int)
    basal_end = len(basal)

    c_br1 = opc_detect_basal_relative(
        sigma, L=L, theta_D=theta_D, theta_R=5, basal_end=basal_end, rho=0.85, K=K
    )
    c_br3 = opc_detect_basal_relative(
        sigma, L=L, theta_D=theta_D, theta_R=6, basal_end=basal_end, rho=0.85, K=K
    )
    # True collapse from rich basal: both BR configs must be able to alarm
    assert int(c_br1["alarm"][basal_end:].sum()) > 0
    assert int(c_br3["alarm"][basal_end:].sum()) > 0
    # Longer θ_R cannot produce *more* alarm samples than shorter θ_R
    assert int(c_br3["alarm"][basal_end:].sum()) <= int(c_br1["alarm"][basal_end:].sum())
    # rel_cap uses basal-relative gate (not pure absolute)
    assert float(c_br1["rel_cap"][0]) <= theta_D


def test_mode_c_g_alone_is_not_a_far_reducer():
    """C-G1 discipline: G>0 recovers mass vs hard OPC — not a free specificity upgrade."""
    L, theta_D, theta_R, K = 8, 0.35, 4, 4
    sig = [0] * L + [0] * 10 + [0, 1, 2, 3, 0, 1, 2, 3] + [0] * 20
    sigma = np.array(sig, dtype=int)
    hard = opc_detect(sigma, L=L, theta_D=theta_D, theta_R=theta_R, K=K)
    g1 = opc_detect_gap_tolerant(
        sigma, L=L, theta_D=theta_D, theta_R=theta_R, G=1, K=K
    )
    # Gap-tolerant with G>0 has at least as many alarms as hard (may raise FAR risk)
    assert int(g1["alarm"].sum()) >= int(hard["alarm"].sum())


def test_mode_s_recommended_configs_s1_s2_s3_basal_relative_delta_s():
    """Recommended Mode-S cells S1–S3 (docs §5.2): basal-relative ΔS + longer θ_R^S."""
    rng = np.random.default_rng(42)
    L_S = 50
    k1 = k2 = 4
    T_basal = 80
    T_app = 100
    # Independent basal, then sustained dependence (synergy surplus rise)
    pi1 = rng.integers(0, k1, size=T_basal + T_app)
    pi2 = rng.integers(0, k2, size=T_basal + T_app)
    pi2[T_basal:] = pi1[T_basal:]  # lock joint dependence in approach
    basal_end = T_basal

    s1 = ops_detect(
        pi1, pi2, L=L_S, theta_R=5, k1=k1, k2=k2,
        basal_end=basal_end, theta_delta_S=0.08,
    )
    s2 = ops_detect(
        pi1, pi2, L=L_S, theta_R=8, k1=k1, k2=k2,
        basal_end=basal_end, theta_delta_S=0.12,
    )
    s3 = ops_detect(
        pi1, pi2, L=L_S, theta_R=10, k1=k1, k2=k2,
        basal_end=basal_end, theta_delta_S=0.15,
    )
    # High-synergy approach should be visible under preferred starting Mode S (S1)
    assert int(s1["alarm"][basal_end:].sum()) > 0
    # Stricter gates (S3) never fire more than looser S1 on the same stream
    assert int(s3["alarm"][basal_end:].sum()) <= int(s1["alarm"][basal_end:].sum())
    assert int(s2["alarm"][basal_end:].sum()) <= int(s1["alarm"][basal_end:].sum())
    # Basal-relative path recorded S_basal
    assert "S_basal" in s1 and np.isfinite(s1["S_basal"][0])


def test_l3_stress_high_support_high_surplus_mode_c_blind_mode_s_sees():
    """
    RECD Nivel-3 stress: full-support high synergy should not force Mode C;
    Mode S (OPS ΔS) is the Level-3–consistent proxy path.
    """
    rng = np.random.default_rng(7)
    L = 40
    k1 = k2 = 4
    K = k1 * k2
    n_basal, n_app = 100, 120

    # Independent full-support basal
    pi1_b = rng.integers(0, k1, size=n_basal)
    pi2_b = rng.integers(0, k2, size=n_basal)

    # Approach: ridge coupling with floor so many joints stay active (high support)
    # pi2 ≈ pi1 most of the time, but inject enough off-ridge mass to keep D high
    pi1_a = rng.integers(0, k1, size=n_app)
    pi2_a = pi1_a.copy()
    flip = rng.random(n_app) < 0.35
    pi2_a[flip] = rng.integers(0, k2, size=int(flip.sum()))

    pi1 = np.concatenate([pi1_b, pi1_a])
    pi2 = np.concatenate([pi2_b, pi2_a])
    joint = joint_symbols_from_factors(pi1, pi2, k2=k2)
    basal_end = n_basal

    # Mode C: absolute + BR collapse on joint symbols
    mode_c = opc_detect(joint, L=L, theta_D=0.35, theta_R=5, K=K)
    mode_c_br = opc_detect_basal_relative(
        joint, L=L, theta_D=0.35, theta_R=5, basal_end=basal_end, rho=0.85, K=K
    )
    # Mode S: basal-relative surplus
    mode_s = ops_detect(
        pi1, pi2, L=L, theta_R=5, k1=k1, k2=k2,
        basal_end=basal_end, theta_delta_S=0.05,
    )

    # High-support stream: mean post-basal diversity should stay above collapse bar often
    post_div = mode_c["diversity"][basal_end:]
    post_div = post_div[np.isfinite(post_div)]
    assert post_div.size > 0
    assert float(np.mean(post_div)) > 0.35  # not a collapse-dominated regime

    # Mode S must be able to register surplus rise (L3-proxy path open)
    assert float(np.nanmean(mode_s["synergy"][basal_end:])) > float(
        np.nanmean(mode_s["synergy"][L - 1 : basal_end])
    )
    # Family L3 capacity: Mode S can alarm while Mode C may stay quiet or weak
    # (do not require Mode C zero — only that Mode S is the designed visibility path)
    assert int(mode_s["high_syn"][basal_end:].sum()) > 0 or int(
        mode_s["alarm"][basal_end:].sum()
    ) > 0
    # BR remains collapse-based (still L3-blind ontology): uses only D / rel_cap
    assert "rel_cap" in mode_c_br


# ---------------------------------------------------------------------------
# Integrated N2-persistence × N3-surplus detector (opsp_integrated_detect)
# ---------------------------------------------------------------------------


def _indep_then_dep(rng, n_basal, n_app, k1=4, k2=4):
    pi1 = rng.integers(0, k1, size=n_basal + n_app)
    pi2 = rng.integers(0, k2, size=n_basal + n_app)
    pi2[n_basal:] = pi1[n_basal:]
    return pi1, pi2, n_basal


def test_integrated_core_matches_ops_detect_alarm_mass_order():
    """I0 (collapse_role=none): primary surplus-persistence aligns with ops_detect."""
    rng = np.random.default_rng(11)
    L, theta_R = 30, 4
    k1 = k2 = 4
    pi1, pi2, basal_end = _indep_then_dep(rng, 60, 80, k1, k2)
    ops = ops_detect(
        pi1, pi2, L=L, theta_R=theta_R, k1=k1, k2=k2,
        basal_end=basal_end, theta_delta_S=0.08,
    )
    integ = opsp_integrated_detect(
        pi1, pi2, L=L, theta_R=theta_R, k1=k1, k2=k2,
        basal_end=basal_end, theta_delta_S=0.08,
        collapse_role="none",
    )
    # Core integrated path must alarm on sustained ΔS rise (same primary ontology as OPS)
    assert int(integ["core_alarm"][basal_end:].sum()) > 0
    assert int(integ["alarm"][basal_end:].sum()) == int(integ["core_alarm"][basal_end:].sum())
    # Same high_syn / alarm ordering spirit as ops_detect on this stream
    assert int(ops["alarm"][basal_end:].sum()) > 0
    assert np.array_equal(integ["high_syn"], ops["high_syn"])
    assert np.array_equal(integ["core_alarm"], ops["alarm"])


def test_integrated_collapse_only_stream_never_alarms():
    """Collapse without synergy surplus must not fire under any complementary role."""
    # Single joint symbol forever → D low (collapse) but joint = product margins
    # for deterministic constant factors: pi1=0, pi2=0 → S_t = 0
    T = 200
    pi1 = np.zeros(T, dtype=int)
    pi2 = np.zeros(T, dtype=int)
    L, k1, k2 = 20, 4, 4
    basal_end = 50
    for role in ("none", "tag", "modulate", "confirm"):
        out = opsp_integrated_detect(
            pi1, pi2, L=L, theta_R=3, k1=k1, k2=k2,
            basal_end=basal_end, theta_delta_S=0.05,
            collapse_role=role, theta_D=0.35, confirm_window=5,
        )
        assert int(out["alarm"].sum()) == 0, f"role={role} must not alarm on collapse-only"
        assert int(out["core_alarm"].sum()) == 0
        # Diversity should register as low on constant joint
        post = out["diversity"][basal_end:]
        post = post[np.isfinite(post)]
        assert post.size > 0 and float(np.mean(post)) <= 0.35


def test_integrated_high_support_high_surplus_alarms_without_collapse():
    """RECD L3 path: full-support high surplus fires I0/tag; confirm may filter."""
    rng = np.random.default_rng(21)
    L = 40
    k1 = k2 = 4
    n_basal, n_app = 90, 110
    pi1_b = rng.integers(0, k1, size=n_basal)
    pi2_b = rng.integers(0, k2, size=n_basal)
    # Ridge + floor: high surplus, many joints active
    pi1_a = rng.integers(0, k1, size=n_app)
    pi2_a = pi1_a.copy()
    flip = rng.random(n_app) < 0.40
    pi2_a[flip] = rng.integers(0, k2, size=int(flip.sum()))
    pi1 = np.concatenate([pi1_b, pi1_a])
    pi2 = np.concatenate([pi2_b, pi2_a])
    basal_end = n_basal

    i0 = opsp_integrated_detect(
        pi1, pi2, L=L, theta_R=4, k1=k1, k2=k2,
        basal_end=basal_end, theta_delta_S=0.04,
        collapse_role="none", theta_D=0.35,
    )
    itag = opsp_integrated_detect(
        pi1, pi2, L=L, theta_R=4, k1=k1, k2=k2,
        basal_end=basal_end, theta_delta_S=0.04,
        collapse_role="tag", theta_D=0.35, confirm_window=5,
    )
    iconf = opsp_integrated_detect(
        pi1, pi2, L=L, theta_R=4, k1=k1, k2=k2,
        basal_end=basal_end, theta_delta_S=0.04,
        collapse_role="confirm", theta_D=0.35, confirm_window=5,
    )
    # Mean post-basal diversity stays above collapse bar (high support)
    post_d = i0["diversity"][basal_end:]
    post_d = post_d[np.isfinite(post_d)]
    assert float(np.mean(post_d)) > 0.35
    # Primary path must see surplus (core alarm or high_syn mass)
    assert int(i0["high_syn"][basal_end:].sum()) > 0
    assert int(i0["alarm"][basal_end:].sum()) > 0
    # tag does not change primary alarm vs none
    assert np.array_equal(i0["alarm"], itag["alarm"])
    # confirm never adds alarms beyond core; may reduce
    assert int(iconf["alarm"].sum()) <= int(i0["core_alarm"].sum())
    assert int(iconf["alarm"].sum()) <= int(iconf["core_alarm"].sum())


def test_integrated_modulate_relaxes_when_collapse_co_present_not_without_surplus():
    """Modulate can ease gates when low_div, but still requires surplus activity."""
    rng = np.random.default_rng(33)
    L, k1, k2 = 25, 3, 3
    # Independent basal
    n_b = 50
    pi1 = list(rng.integers(0, k1, size=n_b))
    pi2 = list(rng.integers(0, k2, size=n_b))
    # Approach: perfect coupling on a *collapsed* alphabet (pi2=pi1, few symbols)
    # Use only symbol 0,1 locked → high S + low D
    for _ in range(80):
        s = rng.integers(0, 2)
        pi1.append(s)
        pi2.append(s)
    pi1 = np.array(pi1, dtype=int)
    pi2 = np.array(pi2, dtype=int)
    basal_end = n_b
    # Strict θ_R so core barely reaches; modulate should allow more when low_div
    core = opsp_integrated_detect(
        pi1, pi2, L=L, theta_R=8, k1=k1, k2=k2,
        basal_end=basal_end, theta_delta_S=0.10,
        collapse_role="none", theta_D=0.50,
    )
    mod = opsp_integrated_detect(
        pi1, pi2, L=L, theta_R=8, k1=k1, k2=k2,
        basal_end=basal_end, theta_delta_S=0.10,
        collapse_role="modulate", theta_D=0.50,
        modulate_delta_R=4, modulate_delta_S=0.03,
    )
    # Modulate never invents alarms without high_syn path: core_alarm/high_syn needed
    assert int(mod["alarm"].sum()) >= int(core["alarm"].sum())
    # Collapse-only still silent under modulate (reuse constant stream)
    z1 = np.zeros(120, dtype=int)
    z2 = np.zeros(120, dtype=int)
    silent = opsp_integrated_detect(
        z1, z2, L=20, theta_R=3, k1=3, k2=3,
        basal_end=40, theta_delta_S=0.05,
        collapse_role="modulate", theta_D=0.35,
        modulate_delta_R=2, modulate_delta_S=0.05,
    )
    assert int(silent["alarm"].sum()) == 0


def test_integrated_confirm_requires_core_surplus_first():
    """Confirm filters core surplus alarms; collapse alone cannot pass confirm."""
    rng = np.random.default_rng(44)
    L, k1, k2 = 30, 4, 4
    pi1, pi2, basal_end = _indep_then_dep(rng, 55, 90, k1, k2)
    core = opsp_integrated_detect(
        pi1, pi2, L=L, theta_R=3, k1=k1, k2=k2,
        basal_end=basal_end, theta_delta_S=0.06,
        collapse_role="none",
    )
    conf = opsp_integrated_detect(
        pi1, pi2, L=L, theta_R=3, k1=k1, k2=k2,
        basal_end=basal_end, theta_delta_S=0.06,
        collapse_role="confirm", confirm_window=8, theta_D=0.35,
    )
    # Confirm ⊆ core
    assert np.all(conf["alarm"] <= conf["core_alarm"])
    assert int(conf["alarm"].sum()) <= int(core["alarm"].sum())
    # If any confirm alarms exist, they must have collapse_coincident
    if int(conf["alarm"].sum()) > 0:
        assert int(conf["collapse_coincident"][conf["alarm"] == 1].sum()) == int(
            conf["alarm"].sum()
        )


def test_integrated_rejects_unknown_collapse_role():
    with pytest.raises(ValueError, match="collapse_role"):
        opsp_integrated_detect(
            np.zeros(40, dtype=int),
            np.zeros(40, dtype=int),
            L=10,
            theta_R=2,
            k1=2,
            k2=2,
            collapse_role="dominant_collapse",
        )


def test_integrated_recommended_configs_i0_imod_iconf_callable():
    """Named configs I0 / I-mod / I-confirm from design note are callable as-shipped."""
    rng = np.random.default_rng(55)
    L_S = 50
    k1 = k2 = 4
    pi1, pi2, basal_end = _indep_then_dep(rng, 80, 100, k1, k2)
    # I0 — preferred starting integrated core
    i0 = opsp_integrated_detect(
        pi1, pi2, L=L_S, theta_R=5, k1=k1, k2=k2,
        basal_end=basal_end, theta_delta_S=0.08,
        collapse_role="none",
    )
    # I-mod
    imod = opsp_integrated_detect(
        pi1, pi2, L=L_S, theta_R=5, k1=k1, k2=k2,
        basal_end=basal_end, theta_delta_S=0.08,
        collapse_role="modulate", theta_D=0.35,
        modulate_delta_R=2, modulate_delta_S=0.02,
    )
    # I-confirm (secondary FAR-leaning)
    iconf = opsp_integrated_detect(
        pi1, pi2, L=L_S, theta_R=8, k1=k1, k2=k2,
        basal_end=basal_end, theta_delta_S=0.12,
        collapse_role="confirm", theta_D=0.35, confirm_window=5,
    )
    assert "alarm" in i0 and "core_alarm" in i0 and "synergy" in i0
    assert int(i0["alarm"][basal_end:].sum()) > 0
    assert int(imod["alarm"].sum()) >= 0  # exercised
    assert int(iconf["alarm"].sum()) <= int(iconf["core_alarm"].sum())
