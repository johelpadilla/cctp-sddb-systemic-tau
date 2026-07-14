"""
Unit + structural tests for the pre-registered I0 surplus-persist param grid.

Locks success rules without Holter data; drives shipped detector on high-surplus
synthetic streams (including strictest grid corner). Full-cohort numbers come
from real runs under results/ when present.
"""
from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
CODE = ROOT / "code"
sys.path.insert(0, str(CODE))

from ordinal_detectors.opc_refinements import opsp_integrated_detect  # noqa: E402


def _load_grid():
    path = CODE / "run_i0_surplus_param_grid.py"
    spec = importlib.util.spec_from_file_location("run_i0_surplus_param_grid", path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_grid_has_16_cells_fixed_l_collapse_none():
    mod = _load_grid()
    grid = mod.build_i0_grid()
    assert len(grid) == 16
    assert [c["theta_delta_S"] for c in grid[:4]] == [0.08, 0.08, 0.08, 0.08]
    assert [c["theta_R"] for c in grid[:4]] == [5, 8, 10, 12]
    assert {c["theta_delta_S"] for c in grid} == {0.08, 0.10, 0.12, 0.15}
    assert {c["theta_R"] for c in grid} == {5, 8, 10, 12}
    for c in grid:
        assert c["L"] == 50
        assert c["collapse_role"] == "none"
        assert c["name"].startswith("I0_")


def test_pre_registered_rules_constants():
    mod = _load_grid()
    rules = mod.PRE_REGISTERED_RULES
    assert rules["require_far_le_multiplier_times_opc"] == 2.0
    assert rules["require_sens_all_ge"] == 0.65
    assert rules["both_required"] is True
    assert abs(rules["opc_sens_all"] - 0.4242424242) < 1e-6
    assert abs(rules["opc_far_per_24h"] - 3.733376994) < 1e-6
    assert abs(rules["far_cap_per_24h"] - 2.0 * rules["opc_far_per_24h"]) < 1e-9
    assert "I-confirm" in " ".join(rules["non_success"]) or "OR" in " ".join(
        rules["non_success"]
    )


def test_score_clear_advance_both_bars_required():
    mod = _load_grid()
    # Clear advance: low FAR and high sens
    ok = mod.score_clear_advance(0.70, 5.0)
    assert ok["clear_advance"] is True
    assert ok["far_ok"] is True
    assert ok["sens_ok"] is True

    # High sens but FAR ~40 (I0 default regime) — NOT success
    high_far = mod.score_clear_advance(0.788, 40.134)
    assert high_far["clear_advance"] is False
    assert high_far["sens_ok"] is True
    assert high_far["far_ok"] is False
    assert "one-sided" in high_far["reason"] or "exceeds" in high_far["reason"]

    # FAR under cap but sens only at OPC level — NOT success
    low_sens = mod.score_clear_advance(0.424, 5.0)
    assert low_sens["clear_advance"] is False
    assert low_sens["far_ok"] is True
    assert low_sens["sens_ok"] is False

    # Border: FAR exactly 2×OPC, sens exactly 0.65
    far_cap = 2.0 * mod.OPC_FAR
    border = mod.score_clear_advance(0.65, far_cap)
    assert border["clear_advance"] is True

    # Just over FAR cap
    over = mod.score_clear_advance(0.70, far_cap + 0.01)
    assert over["clear_advance"] is False


def test_score_grid_rows_attaches_booleans():
    mod = _load_grid()
    rows = [
        {"detector": "a", "sens_all": 0.80, "far_per_24h": 40.0},
        {"detector": "b", "sens_all": 0.70, "far_per_24h": 5.0},
    ]
    scored = mod.score_grid_rows(rows)
    assert scored[0]["clear_advance"] is False
    assert scored[1]["clear_advance"] is True
    assert "score_reason" in scored[0]


def test_strict_corner_still_alarms_on_high_surplus_stream():
    """Strictest grid corner (θ_ΔS=0.15, θ_R=12) must not silence pure high-synergy."""
    rng = np.random.default_rng(11)
    T = 600
    pi1 = rng.integers(0, 6, T)
    pi2 = rng.integers(0, 6, T)
    # Strong dependence after basal: high surplus without needing collapse
    pi1[250:] = rng.integers(0, 6, T - 250)
    pi2[250:] = pi1[250:]
    out = opsp_integrated_detect(
        pi1,
        pi2,
        L=50,
        theta_R=12,
        basal_end=200,
        theta_delta_S=0.15,
        collapse_role="none",
        k1=6,
        k2=6,
    )
    assert int(out["alarm"].sum()) > 0
    assert np.array_equal(out["alarm"], out["core_alarm"])
    # Collapse-only should not be required: low_div may be sparse
    assert "high_syn" in out


def test_loose_and_strict_i0_via_grid_configs():
    mod = _load_grid()
    grid = mod.build_i0_grid()
    loose = grid[0]  # 0.08, 5
    strict = grid[-1]  # 0.15, 12
    assert loose["theta_delta_S"] == 0.08 and loose["theta_R"] == 5
    assert strict["theta_delta_S"] == 0.15 and strict["theta_R"] == 12

    rng = np.random.default_rng(21)
    T = 500
    pi1 = rng.integers(0, 6, T)
    pi2 = pi1.copy()
    pi2[:80] = rng.integers(0, 6, 80)
    out_loose = opsp_integrated_detect(
        pi1, pi2, L=50, theta_R=loose["theta_R"], basal_end=100,
        theta_delta_S=loose["theta_delta_S"], collapse_role="none", k1=6, k2=6,
    )
    out_strict = opsp_integrated_detect(
        pi1, pi2, L=50, theta_R=strict["theta_R"], basal_end=100,
        theta_delta_S=strict["theta_delta_S"], collapse_role="none", k1=6, k2=6,
    )
    # Strict gates ⊆ loose (monotone non-increase of alarms)
    assert int(out_strict["alarm"].sum()) <= int(out_loose["alarm"].sum())
    assert int(out_loose["alarm"].sum()) > 0


def test_full_grid_artifacts_when_present():
    summary = ROOT / "results" / "i0_surplus_param_grid_summary.json"
    if not summary.exists():
        pytest.skip("no i0 grid summary yet")
    data = json.loads(summary.read_text(encoding="utf-8"))
    if data.get("smoke"):
        pytest.skip("summary is smoke-only")
    assert data["phase"] == "i0_surplus_param_grid"
    rules = data["pre_registered_rules"]
    assert rules["require_sens_all_ge"] == 0.65
    assert rules["require_far_le_multiplier_times_opc"] == 2.0
    grid = data["grid"]
    assert len(grid) == 16
    for cell in grid:
        assert "clear_advance" in cell
        assert "sens_all" in cell
        assert "far_per_24h" in cell
        assert cell["collapse_role"] == "none"
        assert 0.0 <= float(cell["sens_all"]) <= 1.0
        assert np.isfinite(float(cell["far_per_24h"]))
    assert "n_clear_advance" in data
    assert data["n_clear_advance"] == sum(1 for c in grid if c["clear_advance"])
    csv_path = ROOT / "results" / "i0_surplus_param_grid.csv"
    assert csv_path.exists() and csv_path.stat().st_size > 100
