#!/usr/bin/env python3
"""
Unit tests for exploratory SDD→OPC cascade merger (shipped pure function).

Drives cascade_sdd_confirm_opc / cascade_first_alarm_index — not a reimplementation.
Does not retune opc_detect / sdd_detect.
"""
from __future__ import annotations

import ast
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
CODE = ROOT / "code"
sys.path.insert(0, str(CODE))

from ordinal_detectors.cascade_fusion import (  # noqa: E402
    CONFIRM_WINDOW_H,
    CONFIRM_WINDOW_MIN,
    cascade_first_alarm_index,
    cascade_first_causal_detection,
    cascade_sdd_confirm_opc,
)
from ordinal_detectors.opc_detector import opc_detect  # noqa: E402
from ordinal_detectors.sdd_detector import sdd_detect  # noqa: E402


def _oracle_cascade(sdd, opc, t, w, opc_max_hr=None):
    """Independent oracle: closed ±W window; optional causal opc_max_hr (t_opc < max)."""
    sdd = np.asarray(sdd)
    opc = np.asarray(opc)
    t = np.asarray(t, dtype=float)
    out = np.zeros(len(t), dtype=np.int8)
    opc_times = t[np.asarray(opc) != 0]
    if opc_max_hr is not None:
        opc_times = opc_times[opc_times < float(opc_max_hr)]
    for i in range(len(t)):
        if sdd[i] == 0:
            continue
        if opc_times.size and np.any(np.abs(opc_times - t[i]) <= w):
            # decision_time = max(t_sdd, earliest opc in window)
            conf = opc_times[np.abs(opc_times - t[i]) <= w]
            dec = max(float(t[i]), float(np.min(conf)))
            if opc_max_hr is not None and not (dec < float(opc_max_hr)):
                continue
            out[i] = 1
    return out


# ---------------------------------------------------------------------------
# Pure merger edge cases
# ---------------------------------------------------------------------------


def test_confirm_window_default_is_five_minutes():
    assert CONFIRM_WINDOW_MIN == 5.0
    assert abs(CONFIRM_WINDOW_H - 5.0 / 60.0) < 1e-12


def test_sdd_only_no_cascade():
    """SDD spike without OPC in window → no cascade alarm."""
    t = np.array([0.0, 0.1, 0.2, 0.3, 0.4])  # hours
    sdd = np.array([0, 0, 1, 0, 0], dtype=np.int8)
    opc = np.array([0, 0, 0, 0, 0], dtype=np.int8)
    out = cascade_sdd_confirm_opc(sdd, opc, t)
    assert int(out["cascade_alarm_count"]) == 0
    assert np.all(out["alarm"] == 0)
    assert int(out["sdd_candidate_count"]) == 1


def test_opc_only_no_cascade():
    """OPC-only → cascade never fires (SDD is the candidate filter)."""
    t = np.linspace(0, 1, 11)
    sdd = np.zeros(11, dtype=np.int8)
    opc = np.ones(11, dtype=np.int8)
    out = cascade_sdd_confirm_opc(sdd, opc, t)
    assert np.all(out["alarm"] == 0)
    assert int(out["cascade_alarm_count"]) == 0
    assert int(out["opc_alarm_count"]) == 11


def test_sdd_plus_opc_within_window_cascade_on():
    """SDD + OPC within ±5 min → cascade at the SDD sample."""
    # 1 sample per minute
    t = np.arange(0, 30) / 60.0  # 0..29 min in hours
    sdd = np.zeros(30, dtype=np.int8)
    opc = np.zeros(30, dtype=np.int8)
    sdd[15] = 1  # t = 15 min
    opc[18] = 1  # t = 18 min, |Δ| = 3 min < 5
    out = cascade_sdd_confirm_opc(sdd, opc, t, confirm_window_h=CONFIRM_WINDOW_H)
    assert out["alarm"][15] == 1
    assert int(out["cascade_alarm_count"]) == 1
    # Cascade only marks SDD-on indices
    assert out["alarm"][18] == 0


def test_sdd_opc_outside_window_no_cascade():
    """SDD + OPC more than 5 min apart → no cascade."""
    t = np.arange(0, 30) / 60.0
    sdd = np.zeros(30, dtype=np.int8)
    opc = np.zeros(30, dtype=np.int8)
    sdd[10] = 1  # 10 min
    opc[20] = 1  # 20 min, |Δ| = 10 min > 5
    out = cascade_sdd_confirm_opc(sdd, opc, t, confirm_window_h=CONFIRM_WINDOW_H)
    assert np.all(out["alarm"] == 0)


def test_boundary_exactly_five_minutes_included():
    """Closed interval: OPC exactly ±5 min from SDD → cascade ON."""
    t = np.arange(0, 20) / 60.0  # minutes 0..19 as hours
    sdd = np.zeros(20, dtype=np.int8)
    opc = np.zeros(20, dtype=np.int8)
    sdd[10] = 1  # 10 min
    opc[15] = 1  # 15 min → exactly +5 min
    out = cascade_sdd_confirm_opc(sdd, opc, t, confirm_window_h=CONFIRM_WINDOW_H)
    assert out["alarm"][10] == 1

    # Exactly −5 min
    opc2 = np.zeros(20, dtype=np.int8)
    opc2[5] = 1
    out2 = cascade_sdd_confirm_opc(sdd, opc2, t, confirm_window_h=CONFIRM_WINDOW_H)
    assert out2["alarm"][10] == 1


def test_boundary_just_outside_five_minutes_excluded():
    """Just beyond ±5 min → no cascade."""
    # Use fractional minutes so we are clearly outside
    # SDD at 0.0 h; OPC at 5 min + 1 second
    t = np.array([0.0, (5.0 + 1.0 / 60.0) / 60.0])  # 0 and ~5.0167 min
    sdd = np.array([1, 0], dtype=np.int8)
    opc = np.array([0, 1], dtype=np.int8)
    out = cascade_sdd_confirm_opc(sdd, opc, t, confirm_window_h=CONFIRM_WINDOW_H)
    assert out["alarm"][0] == 0


def test_matches_closed_window_oracle():
    """Random streams: shipped merger matches independent closed-window oracle."""
    rng = np.random.default_rng(42)
    n = 200
    t = np.cumsum(rng.uniform(0.001, 0.02, size=n))  # irregular hours
    sdd = (rng.random(n) < 0.15).astype(np.int8)
    opc = (rng.random(n) < 0.12).astype(np.int8)
    w = CONFIRM_WINDOW_H
    out = cascade_sdd_confirm_opc(sdd, opc, t, confirm_window_h=w)
    expected = _oracle_cascade(sdd, opc, t, w)
    np.testing.assert_array_equal(out["alarm"], expected)


def test_post_event_opc_cannot_confirm_pre_event_sdd():
    """Look-ahead bug fix: OPC only after event_hr must not confirm SDD before event."""
    # Timeline in minutes → hours: SDD at 10 min, event at 12 min, OPC at 14 min
    # |14-10|=4 min < 5 → non-causal would confirm; causal must reject
    t = np.array([10.0, 12.0, 14.0]) / 60.0
    sdd = np.array([1, 0, 0], dtype=np.int8)
    opc = np.array([0, 0, 1], dtype=np.int8)
    event_hr = 12.0 / 60.0
    # Non-causal (no opc_max): would fire
    noncausal = cascade_sdd_confirm_opc(sdd, opc, t, confirm_window_h=CONFIRM_WINDOW_H)
    assert noncausal["alarm"][0] == 1
    # Causal: opc_max_hr = event_hr
    causal = cascade_sdd_confirm_opc(
        sdd, opc, t, confirm_window_h=CONFIRM_WINDOW_H, opc_max_hr=event_hr
    )
    assert causal["alarm"][0] == 0
    assert int(causal["cascade_alarm_count"]) == 0


def test_pre_event_opc_confirms_with_decision_time():
    """SDD then OPC both before event → cascade ON; decision_time = max(t_sdd, t_opc)."""
    t = np.array([10.0, 12.0, 14.0, 20.0]) / 60.0  # min
    sdd = np.array([1, 0, 0, 0], dtype=np.int8)
    opc = np.array([0, 1, 0, 0], dtype=np.int8)  # OPC 2 min after SDD, both pre-event
    event_hr = 20.0 / 60.0
    out = cascade_sdd_confirm_opc(
        sdd, opc, t, confirm_window_h=CONFIRM_WINDOW_H, opc_max_hr=event_hr
    )
    assert out["alarm"][0] == 1
    # decision = max(10, 12) min
    assert abs(out["decision_time_hr"][0] - 12.0 / 60.0) < 1e-12


def test_causal_first_detection_rejects_post_event_only_confirm():
    """cascade_first_causal_detection: post-event-only OPC → not alarmed."""
    t = np.array([0.10, 0.20, 0.30, 0.40])  # hours
    sdd = np.array([0, 1, 0, 0], dtype=np.int8)  # SDD at 0.20
    opc = np.array([0, 0, 0, 1], dtype=np.int8)  # OPC at 0.40
    # event between SDD and OPC; |0.40-0.20|=0.2h=12min > 5min anyway
    # tighter: OPC just after event within 5 min of SDD
    t = np.array([0.10, 0.15, 0.18, 0.20])
    sdd = np.array([0, 1, 0, 0], dtype=np.int8)  # 0.15
    opc = np.array([0, 0, 0, 1], dtype=np.int8)  # 0.20; event at 0.18
    # |0.20-0.15|=0.05h=3min < 5min; OPC after event 0.18
    det = cascade_first_causal_detection(
        sdd,
        opc,
        t,
        confirm_window_h=CONFIRM_WINDOW_H,
        search_start_hr=0.0,
        event_hr=0.18,
    )
    assert det["alarmed"] == 0


def test_causal_first_detection_accepts_pre_event_pair():
    t = np.array([0.10, 0.15, 0.16, 0.30])
    sdd = np.array([0, 1, 0, 0], dtype=np.int8)
    opc = np.array([0, 0, 1, 0], dtype=np.int8)
    det = cascade_first_causal_detection(
        sdd,
        opc,
        t,
        confirm_window_h=CONFIRM_WINDOW_H,
        search_start_hr=0.05,
        event_hr=0.30,
    )
    assert det["alarmed"] == 1
    # decision = max(0.15, 0.16) = 0.16
    assert abs(det["detection_hr"] - 0.16) < 1e-12
    assert abs(det["lead_time_h"] - (0.30 - 0.16)) < 1e-12


def test_oracle_with_opc_max_matches_shipped():
    rng = np.random.default_rng(7)
    n = 150
    t = np.cumsum(rng.uniform(0.001, 0.015, size=n))
    sdd = (rng.random(n) < 0.2).astype(np.int8)
    opc = (rng.random(n) < 0.15).astype(np.int8)
    event = float(t[n // 2])
    w = CONFIRM_WINDOW_H
    out = cascade_sdd_confirm_opc(
        sdd, opc, t, confirm_window_h=w, opc_max_hr=event
    )
    expected = _oracle_cascade(sdd, opc, t, w, opc_max_hr=event)
    np.testing.assert_array_equal(out["alarm"], expected)


def test_runner_documents_causal_confirm():
    runner = CODE / "run_ordinal_cascade_fusion.py"
    src = runner.read_text(encoding="utf-8")
    assert "cascade_first_causal_detection" in src
    assert "opc_max_hr" in src or "event_hr" in src
    assert "look-ahead" in src.lower() or "causal" in src.lower()


def test_first_alarm_index_search_window():
    t = np.arange(0, 20) / 60.0
    sdd = np.zeros(20, dtype=np.int8)
    opc = np.zeros(20, dtype=np.int8)
    sdd[5] = 1
    opc[6] = 1
    sdd[12] = 1
    opc[12] = 1
    idx, out = cascade_first_alarm_index(
        sdd, opc, t, search_start=0, search_end=10
    )
    assert idx == 5
    idx2, _ = cascade_first_alarm_index(
        sdd, opc, t, search_start=10, search_end=20
    )
    assert idx2 == 12
    idx3, _ = cascade_first_alarm_index(
        sdd, opc, t, search_start=6, search_end=10
    )
    assert idx3 is None  # cascade only at SDD index 5, outside search


def test_shape_mismatch_raises():
    with pytest.raises(ValueError):
        cascade_sdd_confirm_opc([1, 0], [1], [0.0, 0.1])


def test_negative_window_raises():
    with pytest.raises(ValueError):
        cascade_sdd_confirm_opc([1], [1], [0.0], confirm_window_h=-0.01)


def test_cascade_does_not_modify_singleton_detectors():
    """opc_detect / sdd_detect remain independent; cascade is post-process only."""
    sigma = np.array([0, 1, 2, 0, 1, 2, 0, 1, 2, 3] * 20, dtype=np.int64)
    opc = opc_detect(sigma, L=8, theta_D=0.35, theta_R=5, K=6)
    sdd = sdd_detect(sigma, (0, 30), L_c=10, theta_TV=0.2, theta_S=1, K=6)
    t = np.arange(len(sigma), dtype=float) / 3600.0
    before_opc = opc["alarm"].copy()
    before_sdd = sdd["alarm"].copy()
    _ = cascade_sdd_confirm_opc(sdd["alarm"], opc["alarm"], t)
    np.testing.assert_array_equal(opc["alarm"], before_opc)
    np.testing.assert_array_equal(sdd["alarm"], before_sdd)


def test_runner_source_documents_cascade_as_exploratory():
    """Structural: evaluation entry exists and declares exploratory cascade params."""
    runner = CODE / "run_ordinal_cascade_fusion.py"
    assert runner.is_file(), "missing cascade evaluation entry point"
    src = runner.read_text(encoding="utf-8")
    assert "cascade_sdd_confirm_opc" in src
    assert "CONFIRM_WINDOW" in src or "confirm_window" in src
    assert "exploratory" in src.lower()
    # Must still call singleton detectors (not reimplement)
    assert "opc_detect" in src and "sdd_detect" in src
    # Parses as Python
    ast.parse(src)


def test_singleton_detectors_untouched_source():
    """opc_detector / sdd_detector source still declare no fusion internals."""
    opc_src = (CODE / "ordinal_detectors" / "opc_detector.py").read_text()
    sdd_src = (CODE / "ordinal_detectors" / "sdd_detector.py").read_text()
    assert "Does NOT fuse" in opc_src
    assert "Does NOT fuse" in sdd_src
    # Cascade module is separate
    assert (CODE / "ordinal_detectors" / "cascade_fusion.py").is_file()
