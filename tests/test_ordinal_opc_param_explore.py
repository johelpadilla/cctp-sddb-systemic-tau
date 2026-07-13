#!/usr/bin/env python3
"""
Tests for modest OPC (L, θ_D, θ_R) parameter exploration.

Drives shipped helpers: grid construction, sensitivity aggregation,
FAR pooling, baseline flagging, recommendation, and real opc_detect
on synthetic symbol streams. Does not mock the unit under test.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
CODE = ROOT / "code"
RES = ROOT / "results"
sys.path.insert(0, str(CODE))

from cctp_metrics_core import (  # noqa: E402
    FROZEN_MIN_CONSECUTIVE,
    FROZEN_Z_THRESHOLD,
    count_binary_alarm_episodes,
)
from ordinal_detectors.opc_detector import opc_detect  # noqa: E402
from run_ordinal_opc_param_explore import (  # noqa: E402
    BASELINE_L,
    BASELINE_THETA_D,
    BASELINE_THETA_R,
    DEFAULT_L_VALUES,
    DEFAULT_THETA_D_VALUES,
    DEFAULT_THETA_R_VALUES,
    FAR_SLACK_ABS,
    FAR_SLACK_FACTOR,
    aggregate_far_from_episodes,
    aggregate_sensitivity_from_flags,
    baseline_deltas,
    build_param_grid,
    build_recommendation,
    classify_cell,
    evaluate_grid_cell,
    evaluate_opc_on_control,
    evaluate_opc_on_event,
    far_within_slack,
    finalize_cell_rows,
    is_baseline,
    param_key,
    qualitative_parameter_effects,
)


# ---------------------------------------------------------------------------
# Grid / baseline
# ---------------------------------------------------------------------------


def test_default_grid_has_36_cells_and_includes_baseline():
    grid = build_param_grid()
    assert len(grid) == (
        len(DEFAULT_L_VALUES)
        * len(DEFAULT_THETA_D_VALUES)
        * len(DEFAULT_THETA_R_VALUES)
    )
    assert len(grid) == 36
    assert (BASELINE_L, BASELINE_THETA_D, BASELINE_THETA_R) in grid
    # axes as declared
    assert DEFAULT_L_VALUES == (40, 50, 60, 70)
    assert DEFAULT_THETA_D_VALUES == (0.30, 0.35, 0.40)
    assert DEFAULT_THETA_R_VALUES == (4, 5, 6)
    # unique
    assert len(set(grid)) == 36


def test_is_baseline_and_param_key():
    assert is_baseline(50, 0.35, 5)
    assert not is_baseline(40, 0.35, 5)
    assert not is_baseline(50, 0.30, 5)
    assert not is_baseline(50, 0.35, 4)
    assert param_key(50, 0.35, 5) == "L50_D0.35_R5"


def test_custom_grid_subset():
    g = build_param_grid(
        L_values=(50,), theta_D_values=(0.35,), theta_R_values=(4, 5)
    )
    assert g == [(50, 0.35, 4), (50, 0.35, 5)]


# ---------------------------------------------------------------------------
# Aggregation helpers (real math)
# ---------------------------------------------------------------------------


def test_aggregate_sensitivity_from_flags():
    flags = [
        {"source": "sddb", "alarmed": 1},
        {"source": "sddb", "alarmed": 0},
        {"source": "vfdb", "alarmed": 1},
        {"source": "vfdb", "alarmed": 1},
        {"source": "vfdb", "alarmed": 0},
    ]
    out = aggregate_sensitivity_from_flags(flags)
    assert out["sddb"]["n"] == 2
    assert out["sddb"]["n_detected"] == 1
    assert abs(out["sddb"]["sensitivity"] - 0.5) < 1e-12
    assert out["vfdb"]["n"] == 3
    assert out["vfdb"]["n_detected"] == 2
    assert abs(out["vfdb"]["sensitivity"] - 2 / 3) < 1e-12
    assert out["all_events"]["n"] == 5
    assert out["all_events"]["n_detected"] == 3
    assert abs(out["all_events"]["sensitivity"] - 0.6) < 1e-12


def test_aggregate_far_from_episodes_uses_shipped_formula():
    rows = [
        {"n_episodes": 4.0, "search_hours": 10.0, "alarmed": 1.0},
        {"n_episodes": 2.0, "search_hours": 14.0, "alarmed": 1.0},
    ]
    far = aggregate_far_from_episodes(rows)
    # (4+2)/(10+14)*24 = 6
    assert abs(far["far_per_24h"] - 6.0) < 1e-9
    assert far["total_episodes"] == 6.0
    assert far["n_controls"] == 2


def test_baseline_deltas_and_slack():
    d = baseline_deltas(0.5, 5.0, 0.4, 4.0)
    assert abs(d["delta_sens_all"] - 0.1) < 1e-12
    assert abs(d["delta_far_per_24h"] - 1.0) < 1e-12
    assert far_within_slack(5.0, 4.0)  # within max(1.5*4, 4+2)=6
    assert not far_within_slack(10.0, 4.0)
    assert far_within_slack(3.0, 4.0)  # better FAR


def test_classify_and_recommendation_keep_baseline_when_no_gain():
    rows = [
        {
            "param_key": "L50_D0.35_R5",
            "L": 50,
            "theta_D": 0.35,
            "theta_R": 5,
            "is_baseline": True,
            "sens_all": 0.42,
            "far_per_24h": 3.73,
            "delta_sens_all": 0.0,
            "delta_far_per_24h": 0.0,
            "class": "baseline",
        },
        {
            "param_key": "L40_D0.40_R4",
            "L": 40,
            "theta_D": 0.40,
            "theta_R": 4,
            "is_baseline": False,
            "sens_all": 0.50,
            "far_per_24h": 20.0,  # far worse
            "delta_sens_all": 0.08,
            "delta_far_per_24h": 16.27,
            "class": "higher_sens_far_worse",
        },
        {
            "param_key": "L70_D0.30_R6",
            "L": 70,
            "theta_D": 0.30,
            "theta_R": 6,
            "is_baseline": False,
            "sens_all": 0.30,
            "far_per_24h": 1.0,
            "delta_sens_all": -0.12,
            "delta_far_per_24h": -2.73,
            "class": "lower_sens_lower_far",
        },
    ]
    # fix class via shipped classifier for middle row
    rows[1]["class"] = classify_cell(0.50, 20.0, 0.42, 3.73)
    assert rows[1]["class"] == "higher_sens_far_worse"
    rec = build_recommendation(rows)
    assert rec["decision"] == "keep_baseline"
    assert rec["exploratory_only"] is True
    assert rec["clinical_claim"] is False
    assert rec["candidates"] == []


def test_recommendation_consider_adopt_when_higher_sens_far_ok():
    rows = [
        {
            "param_key": "L50_D0.35_R5",
            "L": 50,
            "theta_D": 0.35,
            "theta_R": 5,
            "is_baseline": True,
            "sens_all": 0.42,
            "far_per_24h": 3.73,
            "delta_sens_all": 0.0,
            "delta_far_per_24h": 0.0,
            "class": "baseline",
        },
        {
            "param_key": "L40_D0.35_R4",
            "L": 40,
            "theta_D": 0.35,
            "theta_R": 4,
            "is_baseline": False,
            "sens_all": 0.48,
            "far_per_24h": 4.0,
            "delta_sens_all": 0.06,
            "delta_far_per_24h": 0.27,
            "class": classify_cell(0.48, 4.0, 0.42, 3.73),
        },
    ]
    assert rows[1]["class"] == "higher_sens_far_ok"
    rec = build_recommendation(rows)
    assert rec["decision"] == "consider_adopt_exploratory"
    assert rec["adopt_params"]["L"] == 40
    assert rec["candidates"][0]["param_key"] == "L40_D0.35_R4"
    assert rec["clinical_claim"] is False


def test_finalize_cell_rows_flags_baseline():
    raw = [
        {
            "L": 50,
            "theta_D": 0.35,
            "theta_R": 5,
            "param_key": "L50_D0.35_R5",
            "is_baseline": True,
            "n_sddb": 11,
            "n_detected_sddb": 6,
            "sens_sddb": 6 / 11,
            "n_vfdb": 22,
            "n_detected_vfdb": 8,
            "sens_vfdb": 8 / 22,
            "n_all": 33,
            "n_detected_all": 14,
            "sens_all": 14 / 33,
            "n_controls": 18,
            "total_episodes": 28.0,
            "total_search_hours": 180.0,
            "far_per_24h": 28 / 180 * 24,
            "fraction_alarmed": 0.4,
        },
        {
            "L": 40,
            "theta_D": 0.40,
            "theta_R": 4,
            "param_key": "L40_D0.40_R4",
            "is_baseline": False,
            "n_sddb": 11,
            "n_detected_sddb": 7,
            "sens_sddb": 7 / 11,
            "n_vfdb": 22,
            "n_detected_vfdb": 9,
            "sens_vfdb": 9 / 22,
            "n_all": 33,
            "n_detected_all": 16,
            "sens_all": 16 / 33,
            "n_controls": 18,
            "total_episodes": 40.0,
            "total_search_hours": 180.0,
            "far_per_24h": 40 / 180 * 24,
            "fraction_alarmed": 0.5,
        },
    ]
    rows = finalize_cell_rows(raw)
    assert rows[0]["is_baseline"] is True
    assert rows[0]["class"] == "baseline"
    assert abs(rows[0]["delta_sens_all"]) < 1e-12
    assert rows[1]["delta_sens_all"] > 0


# ---------------------------------------------------------------------------
# Real shipped OPC on synthetic streams (smoke path)
# ---------------------------------------------------------------------------


def _collapse_stream(T: int = 200, K: int = 36, collapse_start: int = 100) -> np.ndarray:
    """Diverse then collapsed to few symbols (should trip OPC at moderate params)."""
    rng = np.random.default_rng(0)
    sigma = rng.integers(0, K, size=T)
    # force collapse: only symbols 0,1
    sigma[collapse_start:] = rng.integers(0, 2, size=T - collapse_start)
    return sigma.astype(np.int64)


def test_opc_detect_real_path_collapse_vs_diverse():
    K = 36
    sigma = _collapse_stream()
    out_strict = opc_detect(sigma, L=50, theta_D=0.30, theta_R=6, K=K)
    out_loose = opc_detect(sigma, L=50, theta_D=0.40, theta_R=4, K=K)
    # loose should fire at least as often
    assert int(out_loose["alarm"].sum()) >= int(out_strict["alarm"].sum())
    assert out_loose["alarm"].sum() > 0


def test_evaluate_opc_on_event_and_control_synthetic():
    K = 36
    sigma = _collapse_stream(T=300)
    t_hr = np.linspace(0, 5, len(sigma))
    # event at t=4.5; basal ends at 1.0; collapse starts mid
    event_pack = {
        "source": "sddb",
        "record": "syn",
        "sigma": sigma,
        "K": K,
        "t_hr_sym": t_hr,
        "event_hr": 4.5,
        "basal": (0.25, 1.0),
        "approach": (1.0, 4.5),
        "stratum": "test",
        "search_start": int(np.searchsorted(t_hr, 1.0)),
        "search_end": int(np.searchsorted(t_hr, 4.5)),
    }
    hit = evaluate_opc_on_event(event_pack, L=50, theta_D=0.35, theta_R=5)
    assert hit["source"] == "sddb"
    assert hit["alarmed"] in (0, 1)
    if hit["alarmed"]:
        assert hit["lead_time_h"] > 0

    control_pack = {
        "record_id": "syn_ctrl",
        "sigma": sigma,
        "K": K,
        "t_hr_sym": t_hr,
        "basal": (0.25, 1.0),
        "search_start": 1.0,
        "search_end": 5.0,
        "total_hours_used": 5.0,
    }
    ep = evaluate_opc_on_control(
        control_pack, L=50, theta_D=0.35, theta_R=5, refractory_h=0.5
    )
    assert "n_episodes" in ep
    assert "search_hours" in ep
    assert abs(ep["search_hours"] - 4.0) < 1e-6
    # episode count consistent with shipped binary counter on same alarm
    alarm = opc_detect(sigma, L=50, theta_D=0.35, theta_R=5, K=K)["alarm"]
    direct = count_binary_alarm_episodes(
        alarm, t_hr, search_start_hr=1.0, search_end_hr=5.0, refractory_h=0.5
    )
    assert ep["n_episodes"] == direct["n_episodes"]


def test_evaluate_grid_cell_two_params_synthetic():
    K = 36
    sigma = _collapse_stream(T=250)
    t_hr = np.linspace(0, 4, len(sigma))
    event_packs = [
        {
            "source": "sddb",
            "record": "a",
            "sigma": sigma,
            "K": K,
            "t_hr_sym": t_hr,
            "event_hr": 3.5,
            "basal": (0.25, 0.8),
            "approach": (0.8, 3.5),
            "stratum": "t",
            "search_start": int(np.searchsorted(t_hr, 0.8)),
            "search_end": int(np.searchsorted(t_hr, 3.5)),
        },
        {
            "source": "vfdb",
            "record": "b",
            "sigma": sigma,
            "K": K,
            "t_hr_sym": t_hr,
            "event_hr": 3.5,
            "basal": (0.25, 0.8),
            "approach": (0.8, 3.5),
            "stratum": "t",
            "search_start": int(np.searchsorted(t_hr, 0.8)),
            "search_end": int(np.searchsorted(t_hr, 3.5)),
        },
    ]
    control_packs = [
        {
            "record_id": "c1",
            "sigma": sigma,
            "K": K,
            "t_hr_sym": t_hr,
            "basal": (0.25, 0.8),
            "search_start": 0.8,
            "search_end": 4.0,
            "total_hours_used": 4.0,
        }
    ]
    cell_a = evaluate_grid_cell(
        event_packs, control_packs, L=50, theta_D=0.35, theta_R=5
    )
    cell_b = evaluate_grid_cell(
        event_packs, control_packs, L=50, theta_D=0.40, theta_R=4
    )
    assert cell_a["n_all"] == 2
    assert cell_a["n_controls"] == 1
    assert cell_a["is_baseline"] is True
    assert cell_b["is_baseline"] is False
    # looser params should not have lower sens on this collapse stream
    assert cell_b["sens_all"] >= cell_a["sens_all"] - 1e-12


def test_qualitative_notes_nonempty():
    rows = finalize_cell_rows(
        [
            {
                "L": 50,
                "theta_D": 0.35,
                "theta_R": 5,
                "param_key": "L50_D0.35_R5",
                "is_baseline": True,
                "n_sddb": 1,
                "n_detected_sddb": 1,
                "sens_sddb": 1.0,
                "n_vfdb": 1,
                "n_detected_vfdb": 0,
                "sens_vfdb": 0.0,
                "n_all": 2,
                "n_detected_all": 1,
                "sens_all": 0.5,
                "n_controls": 1,
                "total_episodes": 1.0,
                "total_search_hours": 10.0,
                "far_per_24h": 2.4,
                "fraction_alarmed": 1.0,
            },
            {
                "L": 60,
                "theta_D": 0.30,
                "theta_R": 6,
                "param_key": "L60_D0.30_R6",
                "is_baseline": False,
                "n_sddb": 1,
                "n_detected_sddb": 0,
                "sens_sddb": 0.0,
                "n_vfdb": 1,
                "n_detected_vfdb": 0,
                "sens_vfdb": 0.0,
                "n_all": 2,
                "n_detected_all": 0,
                "sens_all": 0.0,
                "n_controls": 1,
                "total_episodes": 0.0,
                "total_search_hours": 10.0,
                "far_per_24h": 0.0,
                "fraction_alarmed": 0.0,
            },
        ]
    )
    notes = qualitative_parameter_effects(rows)
    assert len(notes) >= 3
    assert any("θ_D" in n or "theta" in n.lower() or "θ_D" in n for n in notes)


def test_absz_frozen_not_retuned():
    """Exploration must not touch abs-z frozen constants."""
    assert FROZEN_Z_THRESHOLD == 2.0
    assert FROZEN_MIN_CONSECUTIVE == 3
    # runner source must not redefine abs-z thresholds as tunable
    src = (CODE / "run_ordinal_opc_param_explore.py").read_text(encoding="utf-8")
    assert "FROZEN_Z_THRESHOLD" in src
    assert "retuned" in src.lower() or "retune" in src.lower()
    assert "opc_detect" in src
    assert "count_binary_alarm_episodes" in src
    # must not call abs-z lead detector as primary exploration path
    assert "detect_lead_time" not in src


def test_far_slack_constants_documented():
    assert FAR_SLACK_FACTOR == 1.5
    assert FAR_SLACK_ABS == 2.0


def test_existing_baseline_artifacts_still_match_expected():
    """Sanity: prior bake-off / FAR baseline references still on disk."""
    exp = RES / "ordinal_exploratory_summary.json"
    far = RES / "ordinal_nsrdb_far_summary.json"
    if not exp.is_file() or not far.is_file():
        pytest.skip("prior ordinal artifacts missing")
    e = json.loads(exp.read_text(encoding="utf-8"))
    f = json.loads(far.read_text(encoding="utf-8"))
    opc = e["opc_L50_companion"]["all_events"]
    assert opc["n"] == 33
    assert opc["n_detected"] == 14
    assert abs(opc["sensitivity"] - 14 / 33) < 1e-12
    pf = f["pooled_far"]["opc_L50"]
    assert pf["n_controls"] == 18
    assert abs(pf["far_per_24h"] - 3.7333769938233505) < 1e-9
