#!/usr/bin/env python3
"""
Unit tests for Phase-1 FAR / short-window helpers in cctp_metrics_core.

Drives the real shipped functions (no reimplementation of detector math).
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

CODE = Path(__file__).resolve().parents[1] / "code"
sys.path.insert(0, str(CODE))

from cctp_metrics_core import (  # noqa: E402
    FROZEN_MIN_CONSECUTIVE,
    FROZEN_THETA3,
    FROZEN_Z_THRESHOLD,
    W_TAU,
    count_alarm_episodes,
    detect_lead_time,
    false_alarm_rate,
    short_db_windows,
)


def test_frozen_constants_match_jul12_discovery():
    assert FROZEN_THETA3 == 0.08
    assert FROZEN_Z_THRESHOLD == 2.0
    assert FROZEN_MIN_CONSECUTIVE == 3
    assert W_TAU == 101


def test_short_db_windows_thirds_for_sub_6h():
    event_hr, basal, approach, stratum = short_db_windows(0.5, 0.58)
    assert stratum in ("short_15_60min", "short_ge60min")
    assert basal[0] < basal[1]
    assert approach[0] < approach[1]
    assert basal[1] <= approach[0] + 1e-9 or basal[1] < approach[0] + 0.2
    # pre=0.5h → short_15_60min
    assert stratum == "short_15_60min"
    assert approach[1] == event_hr


def test_short_db_windows_holter_ge6h():
    event_hr, basal, approach, stratum = short_db_windows(10.0, 20.0)
    assert stratum == "holter_ge6h"
    assert approach[1] == event_hr
    assert basal[1] < approach[0]


def test_short_db_windows_too_short():
    _, _, _, stratum = short_db_windows(0.1, 0.2, min_pre_event_h=0.25)
    assert stratum == "too_short"


def test_count_alarm_episodes_on_flat_control():
    rng = np.random.default_rng(0)
    n = 2000
    t_hr = np.linspace(0, 10, n)
    y = rng.normal(0, 0.2, n)
    out = count_alarm_episodes(
        y, t_hr, basal=(0.5, 2.5), z_threshold=5.0, min_consecutive=5, use_abs=True
    )
    assert out["n_episodes"] == 0.0
    assert out["alarmed"] == 0.0
    assert out["search_hours"] > 0


def test_count_alarm_episodes_finds_spike_run():
    n = 1000
    t_hr = np.linspace(0, 5, n)
    y = np.zeros(n)
    # large run after basal
    y[600:620] = 10.0
    out = count_alarm_episodes(
        y,
        t_hr,
        basal=(0.2, 1.0),
        z_threshold=2.0,
        min_consecutive=3,
        use_abs=True,
        refractory_h=0.5,
    )
    assert out["alarmed"] == 1.0
    assert out["n_episodes"] >= 1.0
    assert np.isfinite(out["first_alarm_hr"])


def test_false_alarm_rate_nan_without_controls():
    fr = false_alarm_rate([])
    assert fr["n_controls"] == 0.0
    assert np.isnan(fr["far_per_24h"])
    assert fr["reason"] == "no_controls"


def test_false_alarm_rate_from_episodes():
    rows = [
        {"n_episodes": 2.0, "search_hours": 12.0, "alarmed": 1.0},
        {"n_episodes": 0.0, "search_hours": 12.0, "alarmed": 0.0},
    ]
    fr = false_alarm_rate(rows)
    # 2 episodes / 24 h * 24 = 2.0 alarms/24h
    assert abs(fr["far_per_24h"] - 2.0) < 1e-9
    assert fr["reason"] == "ok"
    assert fr["n_controls"] == 2.0


def test_detect_lead_time_still_compatible_with_short_windows():
    """Smoke: short basal/approach still drives real detect_lead_time."""
    n = 800
    t_hr = np.linspace(0, 0.5, n)
    y = np.zeros(n)
    y[500:520] = 5.0
    event_hr, basal, approach, _ = short_db_windows(0.5, 0.5)
    out = detect_lead_time(
        y, t_hr, event_hr, basal, z_threshold=2.0, min_consecutive=3, use_abs=True
    )
    # may or may not alarm depending on basal variance of zeros — ensure no crash
    assert "alarmed" in out
    assert "lead_time_h" in out
