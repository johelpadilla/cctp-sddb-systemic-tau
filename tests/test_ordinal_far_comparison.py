#!/usr/bin/env python3
"""
Tests for NSRDB ordinal FAR comparison — drive shipped episode counting and runner path.

No reimplementation of OPC/SDD math; no fusion paths.
"""
from __future__ import annotations

import ast
import json
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
CODE = ROOT / "code"
sys.path.insert(0, str(CODE))

from cctp_metrics_core import (  # noqa: E402
    FROZEN_MIN_CONSECUTIVE,
    FROZEN_Z_THRESHOLD,
    count_alarm_episodes,
    count_binary_alarm_episodes,
    false_alarm_rate,
)
from ordinal_detectors.opc_detector import opc_detect  # noqa: E402
from ordinal_detectors.sdd_detector import sdd_detect  # noqa: E402
from run_ordinal_far_comparison import (  # noqa: E402
    CONTROL_MAX_HOURS,
    K_JOINT,
    OPC_L,
    OPC_THETA_D,
    OPC_THETA_R,
    REFRACTORY_H,
    SDD_L_C,
    SDD_THETA_S,
    SDD_THETA_TV,
    joint_bivariate_symbols,
    list_nsrdb_npz,
    process_control,
)


# ---------------------------------------------------------------------------
# Episode counting / FAR math (shipped functions)
# ---------------------------------------------------------------------------


def test_binary_zero_alarms_zero_episodes():
    t = np.linspace(0, 10, 1001)
    alarm = np.zeros(len(t), dtype=np.int8)
    out = count_binary_alarm_episodes(
        alarm, t, search_start_hr=2.0, search_end_hr=10.0, refractory_h=0.5
    )
    assert out["n_episodes"] == 0.0
    assert out["alarmed"] == 0.0
    assert abs(out["search_hours"] - 8.0) < 1e-9
    assert not np.isfinite(out["first_alarm_hr"])


def test_binary_sustained_run_with_refractory_hand_checked():
    """
    Search (1.0, 5.0). Alarms True at t=1.1, 1.2, 1.3 (same episode after first),
    then True again at t=2.0 (still in 0.5h refractory from 1.1 → until 1.6),
    then True at t=1.7 (wait — 1.7 > 1.6) and t=3.0.

    Timeline (refractory 0.5 h from episode start time):
      t=1.1 ON → episode 1, refractory until 1.6
      t=1.2, 1.3 ON but < 1.6 → suppressed
      t=1.7 ON → episode 2, refractory until 2.2
      t=2.0 ON but < 2.2 → suppressed
      t=3.0 ON → episode 3, refractory until 3.5
    """
    times = np.array([0.5, 1.1, 1.2, 1.3, 1.7, 2.0, 3.0, 4.0, 5.5])
    alarm = np.array([0, 1, 1, 1, 1, 1, 1, 0, 0], dtype=np.int8)
    out = count_binary_alarm_episodes(
        alarm, times, search_start_hr=1.0, search_end_hr=5.0, refractory_h=0.5
    )
    # Points in (1.0, 5.0): 1.1, 1.2, 1.3, 1.7, 2.0, 3.0, 4.0
    assert out["n_episodes"] == 3.0
    assert out["alarmed"] == 1.0
    assert abs(out["first_alarm_hr"] - 1.1) < 1e-9
    assert abs(out["search_hours"] - 4.0) < 1e-9


def test_binary_dense_grid_one_episode_per_refractory_slot():
    """Dense True alarms for 2 h after basal → with 0.5 h refractory ≈ 4 episodes."""
    t = np.linspace(0, 4, 4001)
    alarm = np.zeros(len(t), dtype=np.int8)
    # alarm ON entire search (1.0, 3.0)
    alarm[(t > 1.0) & (t < 3.0)] = 1
    out = count_binary_alarm_episodes(
        alarm, t, search_start_hr=1.0, search_end_hr=3.0, refractory_h=0.5
    )
    # first at ~1.0+, then +0.5 each → ~1.0, 1.5, 2.0, 2.5 → 4 episodes
    assert out["n_episodes"] == 4.0
    assert abs(out["search_hours"] - 2.0) < 1e-9


def test_false_alarm_rate_two_rows_known():
    rows = [
        {"n_episodes": 4.0, "search_hours": 10.0, "alarmed": 1.0},
        {"n_episodes": 2.0, "search_hours": 14.0, "alarmed": 1.0},
    ]
    fr = false_alarm_rate(rows)
    # (4+2)/(10+14)*24 = 6/24*24 = 6.0
    assert abs(fr["far_per_24h"] - 6.0) < 1e-9
    assert fr["total_episodes"] == 6.0
    assert fr["total_search_hours"] == 24.0
    assert fr["reason"] == "ok"


def test_absz_frozen_constants_untouched():
    assert FROZEN_Z_THRESHOLD == 2.0
    assert FROZEN_MIN_CONSECUTIVE == 3
    assert REFRACTORY_H == 0.5
    assert OPC_L == 50
    assert OPC_THETA_D == 0.35
    assert OPC_THETA_R == 5
    assert SDD_L_C == 50
    assert SDD_THETA_TV == 0.35
    assert SDD_THETA_S == 1


def test_count_alarm_episodes_still_uses_frozen_defaults():
    """Regression: continuous path defaults remain frozen."""
    n = 500
    t = np.linspace(0, 5, n)
    y = np.zeros(n)
    y[300:320] = 10.0
    out = count_alarm_episodes(y, t, basal=(0.2, 1.0), use_abs=True, refractory_h=0.5)
    assert out["alarmed"] == 1.0
    assert out["n_episodes"] >= 1.0


# ---------------------------------------------------------------------------
# Real NSRDB path smoke
# ---------------------------------------------------------------------------


def _first_nsrdb() -> Path:
    paths = list_nsrdb_npz()
    if not paths:
        pytest.skip("no NSRDB clean npz")
    return paths[0]


def test_list_nsrdb_expect_18():
    paths = list_nsrdb_npz()
    if not paths:
        pytest.skip("no NSRDB")
    assert len(paths) == 18


def test_real_nsrdb_joint_symbols_and_shipped_detectors():
    path = _first_nsrdb()
    d = np.load(path)
    rr = d["rr_ms"].astype(float)
    # Cap like Phase 2 for speed
    t_hr = d["t_sec"].astype(float) / 3600.0
    keep = t_hr <= 3.0
    rr = rr[keep]
    sigma, K, offset = joint_bivariate_symbols(rr)
    assert K == K_JOINT == 36
    assert len(sigma) > 200
    assert int(sigma.min()) >= 0 and int(sigma.max()) < K

    opc = opc_detect(sigma, L=OPC_L, theta_D=OPC_THETA_D, theta_R=OPC_THETA_R, K=K)
    assert opc["alarm"].shape == sigma.shape

    basal = (0, min(200, len(sigma) // 5))
    sdd = sdd_detect(
        sigma,
        basal,
        L_c=SDD_L_C,
        theta_TV=SDD_THETA_TV,
        theta_S=SDD_THETA_S,
        K=K,
        mask_basal=True,
    )
    assert sdd["alarm"].shape == sigma.shape


def test_process_control_smoke_finite_fields():
    path = _first_nsrdb()
    row = process_control(path, max_hours=3.0, basal_hours=1.0, refractory_h=0.5)
    if row.get("skipped"):
        pytest.skip(f"control skipped: {row.get('skip_reason')}")
    assert row["fusion"] is False
    for prefix in ("opc", "sdd", "absz"):
        assert float(row[f"{prefix}_search_hours"]) >= 0.0
        assert float(row[f"{prefix}_n_episodes"]) >= 0.0
        far = float(row[f"{prefix}_far_per_24h"])
        assert np.isfinite(far) or float(row[f"{prefix}_search_hours"]) == 0.0
    assert row["absz_z_threshold"] == FROZEN_Z_THRESHOLD
    assert row["absz_min_consecutive"] == FROZEN_MIN_CONSECUTIVE


def test_runner_source_no_fusion():
    src_path = CODE / "run_ordinal_far_comparison.py"
    src = src_path.read_text(encoding="utf-8")
    tree = ast.parse(src)
    assert isinstance(tree, ast.Module)
    assert "opc_detect(" in src
    assert "sdd_detect(" in src
    assert "count_binary_alarm_episodes" in src
    assert "count_alarm_episodes" in src
    assert "false_alarm_rate" in src
    # Separate episode fields (not a single fused count)
    assert "opc_n_episodes" in src and "sdd_n_episodes" in src
    assert "absz_n_episodes" in src
    assert "FUSION_FORBIDDEN" in src or "fusion" in src.lower()


def test_control_max_hours_matches_phase2_default():
    assert CONTROL_MAX_HOURS == 12.0
