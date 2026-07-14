"""
Structural + smoke checks for the OPSP Holter evaluation entry path.

Does not hard-code multi-cohort sens/FAR (those come from real runs under results/).
Drives the shipped runner helpers and asserts non-trivial detector output.
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
CODE = ROOT / "code"
sys.path.insert(0, str(CODE))

from ordinal_detectors.opc_refinements import opsp_integrated_detect  # noqa: E402


def _load_runner():
    path = CODE / "run_opsp_integrated_holter_eval.py"
    spec = importlib.util.spec_from_file_location("run_opsp_integrated_holter_eval", path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_runner_module_exports_i0_config_and_helpers():
    mod = _load_runner()
    assert mod.I0["collapse_role"] == "none"
    assert mod.I0["theta_delta_S"] == 0.08
    assert mod.I0["theta_R"] == 5
    assert mod.I0["L"] == 50
    assert mod.I_MOD["collapse_role"] == "modulate"
    assert mod.I_CONFIRM["collapse_role"] == "confirm"
    assert callable(mod.process_event_record)
    assert callable(mod.process_control)
    assert callable(mod.run_opsp)
    assert callable(mod.list_event_records)
    assert callable(mod.list_nsrdb)


def test_consumer_import_opsp_i0_nontrivial_on_dependent_stream():
    """Fresh consumer path: high-synergy stream must produce alarms under I0."""
    rng = np.random.default_rng(7)
    T = 400
    # Independent basal-like segment
    pi1 = rng.integers(0, 6, T)
    pi2 = rng.integers(0, 6, T)
    # Dependent surplus segment
    pi1[200:] = rng.integers(0, 6, T - 200)
    pi2[200:] = pi1[200:]
    out = opsp_integrated_detect(
        pi1,
        pi2,
        L=50,
        theta_R=5,
        basal_end=150,
        theta_delta_S=0.08,
        collapse_role="none",
        k1=6,
        k2=6,
    )
    assert "alarm" in out and "core_alarm" in out
    assert int(out["alarm"].sum()) > 0
    assert np.array_equal(out["alarm"], out["core_alarm"])


def test_i0_identical_to_ops_on_shared_helper_path():
    mod = _load_runner()
    rng = np.random.default_rng(3)
    T = 500
    pi1 = rng.integers(0, 6, T)
    pi2 = pi1.copy()
    pi2[:100] = rng.integers(0, 6, 100)
    basal_end = 120
    i0 = mod.run_opsp(pi1, pi2, mod.I0, basal_end)
    ops = mod.run_ops_s1(pi1, pi2, basal_end)
    assert np.array_equal(i0["alarm"], ops["alarm"])


def test_evaluation_artifacts_exist_when_full_run_present():
    """
    If a full (non-smoke) summary exists, require numeric real-run fields.
    Skips cleanly if only smoke artifacts remain.
    """
    summary = ROOT / "results" / "opsp_integrated_summary.json"
    if not summary.exists():
        pytest.skip("no opsp_integrated_summary.json yet")
    import json

    data = json.loads(summary.read_text(encoding="utf-8"))
    if data.get("smoke"):
        pytest.skip("summary is smoke-only")
    i0 = data["sensitivity"]["I0"]["all_events"]
    far = data["far"]["I0"]
    assert i0["n"] >= 30
    assert 0.0 <= float(i0["sensitivity"]) <= 1.0
    assert far["n_controls"] >= 10
    assert np.isfinite(float(far["far_per_24h"]))
    # Primary remains surplus-persist twin of OPS-S1
    assert data["sensitivity"]["I0"]["all_events"]["n_detected"] == data["sensitivity"]["OPS-S1"]["all_events"]["n_detected"]
    sens_csv = ROOT / "results" / "opsp_integrated_sens_per_record.csv"
    far_csv = ROOT / "results" / "opsp_integrated_nsrdb_far_per_record.csv"
    assert sens_csv.exists() and sens_csv.stat().st_size > 100
    assert far_csv.exists() and far_csv.stat().st_size > 100


def test_report_doc_exists_with_recommendation_language():
    report = ROOT / "docs" / "OPSP_INTEGRATED_HOLTER_EVAL.md"
    assert report.exists()
    text = report.read_text(encoding="utf-8")
    assert "I0" in text
    assert "FAR" in text
    assert "Promote I0" in text or "primary exploratory" in text
    assert "No" in text  # recommendation section uses No for promote
    assert "Level-3" in text or "Nivel-3" in text or "proxy" in text
