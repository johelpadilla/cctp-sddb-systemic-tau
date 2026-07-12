#!/usr/bin/env python3
"""
Unit tests for pure lead-time / detector functions in cctp_metrics_core.

Drives the real shipped functions (no reimplementation, no hard-coded cohort answers).
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

CODE = Path(__file__).resolve().parents[1] / "code"
sys.path.insert(0, str(CODE))

from cctp_metrics_core import (  # noqa: E402
    cumulative_detection_curve,
    detect_lead_time,
    detector_performance,
    event_type_label,
    get_event_and_windows,
    regime_delta,
    resolve_event_timing_from_npz,
    rolling_ar1,
    rolling_var,
    sign_concordance,
    substrate_label,
)


def _synthetic_series(
    n: int = 2000,
    event_i: int = 1800,
    basal_level: float = 0.0,
    approach_level: float = 3.0,
    noise: float = 0.15,
    seed: int = 0,
):
    """Metric series flat in basal, rises near event; t_hr linear 0..10h."""
    rng = np.random.default_rng(seed)
    t_hr = np.linspace(0.0, 10.0, n)
    event_hr = float(t_hr[event_i])
    y = np.full(n, basal_level, dtype=float) + rng.normal(0, noise, n)
    # ramp last ~1.5h before event
    ramp_start = event_i - 300
    for i in range(ramp_start, event_i):
        frac = (i - ramp_start) / max(1, event_i - ramp_start)
        y[i] = basal_level + frac * approach_level + rng.normal(0, noise)
    basal = (0.5, 3.5)
    return y, t_hr, event_hr, basal


def test_detect_lead_time_finds_pre_event_rise():
    y, t_hr, event_hr, basal = _synthetic_series()
    out = detect_lead_time(
        y, t_hr, event_hr, basal, z_threshold=2.0, min_consecutive=3, use_abs=True
    )
    assert out["alarmed"] == 1.0
    assert np.isfinite(out["lead_time_h"])
    assert out["lead_time_h"] > 0.0
    assert out["detection_hr"] < event_hr
    assert out["detection_hr"] > basal[1]


def test_detect_lead_time_no_alarm_on_flat_noise():
    rng = np.random.default_rng(1)
    n = 1500
    t_hr = np.linspace(0, 8, n)
    y = rng.normal(0, 0.2, n)
    event_hr = 7.5
    basal = (0.5, 3.0)
    out = detect_lead_time(
        y, t_hr, event_hr, basal, z_threshold=5.0, min_consecutive=5, use_abs=True
    )
    assert out["alarmed"] == 0.0
    assert not np.isfinite(out["lead_time_h"])


def test_detect_lead_time_respects_min_consecutive():
    t_hr = np.linspace(0, 6, 600)
    y = np.zeros(600)
    # single spike after basal — should not alarm if min_consecutive=3
    y[400] = 10.0
    out = detect_lead_time(
        y, t_hr, event_hr=5.5, basal=(0.2, 2.0), z_threshold=2.0, min_consecutive=3
    )
    assert out["alarmed"] == 0.0
    # three consecutive large deviations
    y[450:453] = 10.0
    out2 = detect_lead_time(
        y, t_hr, event_hr=5.5, basal=(0.2, 2.0), z_threshold=2.0, min_consecutive=3
    )
    assert out2["alarmed"] == 1.0
    assert out2["lead_time_h"] > 0


def test_detector_performance_sensitivity():
    rows = [
        {"alarmed": 1.0, "lead_time_h": 2.0},
        {"alarmed": 1.0, "lead_time_h": 1.0},
        {"alarmed": 0.0, "lead_time_h": float("nan")},
        {"alarmed": 1.0, "lead_time_h": 3.0},
    ]
    perf = detector_performance(rows)
    assert perf["n_records"] == 4.0
    assert perf["n_detected"] == 3.0
    assert abs(perf["sensitivity"] - 0.75) < 1e-9
    assert abs(perf["median_lead_time_h"] - 2.0) < 1e-9
    assert np.isnan(perf["false_alarm_rate"])


def test_cumulative_detection_curve_horizons():
    rows = [
        {"alarmed": 1.0, "lead_time_h": 4.0},
        {"alarmed": 1.0, "lead_time_h": 1.5},
        {"alarmed": 0.0, "lead_time_h": float("nan")},
    ]
    curve = cumulative_detection_curve(rows, horizons_h=(1.0, 2.0, 3.0, 6.0))
    assert len(curve) == 4
    by_h = {c["horizon_h"]: c["detection_rate"] for c in curve}
    assert abs(by_h[1.0] - 2 / 3) < 1e-9
    assert abs(by_h[2.0] - 1 / 3) < 1e-9
    assert abs(by_h[6.0] - 0.0) < 1e-9


def test_sign_concordance():
    c = sign_concordance([1.0, -1.0, 0.5], [2.0, -0.1, -0.5])
    assert c["n"] == 3.0
    assert c["n_concordant"] == 2.0
    assert abs(c["concordance"] - 2 / 3) < 1e-9


def test_regime_delta_and_rolling():
    x = np.concatenate([np.ones(100), np.ones(100) * 3.0])
    t = np.linspace(0, 4, 200)
    d = regime_delta(x, t, basal=(0.0, 1.5), approach=(2.0, 4.0))
    assert abs(d["delta"] - 2.0) < 0.2
    v = rolling_var(np.arange(200, dtype=float), w=21, stride=5)
    assert np.isfinite(v[20::5]).all()
    a = rolling_ar1(np.sin(np.linspace(0, 20, 300)), w=51, stride=5)
    assert np.nanmax(np.abs(a[np.isfinite(a)])) <= 1.0 + 1e-6


def test_get_event_and_windows_record_special_cases():
    t = np.linspace(0, 24, 1000)
    eh, basal, app = get_event_and_windows("35", t, vfon_hr=25.0)
    assert basal == (6.0, 16.0)
    assert app[1] == eh
    eh30, basal30, _ = get_event_and_windows("30", t, vfon_hr=7.91)
    assert basal30 == (0.5, 3.5)
    assert abs(eh30 - 7.91) < 1e-6


def test_substrate_label():
    assert substrate_label("Sinus", "none") == "sinus"
    assert substrate_label("Atrial fibrillation", "none") == "AF"
    assert substrate_label("Sinus", "intermittent") == "paced"


def test_event_type_label_gap_rules():
    """Duration must not be used as event_hr; gap>3h → intermediate."""
    assert event_type_label("99", 10.0, 24.0) == "intermediate"
    assert event_type_label("99", 23.5, 24.0) == "terminal"
    assert event_type_label("99", 22.5, 24.0) == "terminal"  # 1.5h gap band
    assert event_type_label("30", 23.9, 24.0) == "intermediate"  # manuscript override
    # Never treat missing event as duration proxy inside the pure labeler
    assert event_type_label("99", None, 24.0) == "unknown"


def test_resolve_event_timing_from_real_npz_cohort():
    """Drive shipped resolve_event_timing_from_npz on analytic RR npz files."""
    data_dir = Path(__file__).resolve().parents[1] / "data"
    expected = {
        "30": "intermediate",
        "31": "terminal",
        "32": "intermediate",
        "35": "terminal",
        "36": "terminal",
        "38": "intermediate",
        "45": "intermediate",
        "47": "intermediate",
        "50": "intermediate",
        "51": "terminal",
    }
    types = {}
    for rec, exp in expected.items():
        npz = data_dir / f"rr_{rec}_clean.npz"
        assert npz.exists(), f"missing fixture {npz}"
        t = resolve_event_timing_from_npz(rec, str(npz))
        assert np.isfinite(t["event_hr"]), rec
        assert np.isfinite(t["duration_h"]), rec
        # Must not silently set event_hr == duration for intermediate cases with vfon mid-recording
        assert t["event_type"] == exp, (rec, t["event_type"], t["event_hr"], t["duration_h"], t["span_h"])
        types[t["event_type"]] = types.get(t["event_type"], 0) + 1
    assert types.get("intermediate") == 6
    assert types.get("terminal") == 4


def test_load_key_metrics_populates_event_hr():
    """Batch aggregator must write finite event_hr from npz (not empty)."""
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
    from run_cctp_batch import load_key_metrics

    m = load_key_metrics("38")
    assert m.get("event_hr") is not None
    assert float(m["event_hr"]) > 0
    # Record 38 VF is mid-recording (~8h), not at duration (~18h)
    assert float(m["event_hr"]) < float(m["duration_h"]) - 3.0
    assert m.get("event_type") == "intermediate"


if __name__ == "__main__":
    # allow `python tests/test_leadtime_detector.py` without pytest
    import traceback

    tests = [v for k, v in list(globals().items()) if k.startswith("test_") and callable(v)]
    failed = 0
    for fn in tests:
        try:
            fn()
            print(f"PASS {fn.__name__}")
        except Exception as e:
            failed += 1
            print(f"FAIL {fn.__name__}: {e}")
            traceback.print_exc()
    if failed:
        raise SystemExit(1)
    print(f"All {len(tests)} tests passed.")
