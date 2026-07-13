#!/usr/bin/env python3
"""
Tests for ordinal sensitivity–specificity trade-off analysis.

Drives the shipped join runner against real existing result artifacts.
Does not reprocess .wfdb or reimplement detectors.
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
CODE = ROOT / "code"
RES = ROOT / "results"
sys.path.insert(0, str(CODE))

from run_ordinal_tradeoff_analysis import (  # noqa: E402
    DISPLAY_NAMES,
    FORBIDDEN_CLAIM_PHRASES,
    PRIMARY_ARMS,
    balance_metrics,
    build_cohort_summary_rows,
    build_tradeoff_rows,
    extract_far_block,
    extract_sensitivity_block,
    load_json,
    qualitative_notes,
    rank_detectors,
    run_analysis,
    sum_episodes_from_per_record,
)


REQUIRED_INPUTS = [
    RES / "ordinal_exploratory_summary.json",
    RES / "ordinal_opc_sdd_absz_comparison.csv",
    RES / "ordinal_nsrdb_far_summary.json",
    RES / "ordinal_nsrdb_far_per_record.csv",
]


@pytest.fixture(scope="module")
def exploratory():
    path = RES / "ordinal_exploratory_summary.json"
    if not path.is_file():
        pytest.skip(f"missing {path}")
    return load_json(path)


@pytest.fixture(scope="module")
def far_summary():
    path = RES / "ordinal_nsrdb_far_summary.json"
    if not path.is_file():
        pytest.skip(f"missing {path}")
    return load_json(path)


@pytest.fixture(scope="module")
def per_record_rows():
    path = RES / "ordinal_nsrdb_far_per_record.csv"
    if not path.is_file():
        pytest.skip(f"missing {path}")
    with open(path, encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f))


def test_required_inputs_exist():
    missing = [p for p in REQUIRED_INPUTS if not p.is_file()]
    assert not missing, f"Missing required inputs: {missing}"


def test_primary_arms_are_three_and_separate():
    assert len(PRIMARY_ARMS) == 3
    ids = [a[0] for a in PRIMARY_ARMS]
    assert ids == ["opc_L50", "sdd", "absz_tau_s"]
    # OPC L=8 must not be a primary arm
    assert "opc" not in ids or "opc_L8" not in ids
    assert "opc_default" not in ids
    assert "opc_L8" not in ids


def test_extract_sensitivity_matches_json(exploratory):
    opc = extract_sensitivity_block(exploratory, "opc_L50_companion")
    sdd = extract_sensitivity_block(exploratory, "sdd")
    absz = extract_sensitivity_block(exploratory, "absz_tau_s")

    assert opc["sddb"]["n"] == 11
    assert opc["vfdb"]["n"] == 22
    assert opc["all_events"]["n"] == 33
    assert opc["sddb"]["n_detected"] == 6
    assert opc["vfdb"]["n_detected"] == 8
    assert opc["all_events"]["n_detected"] == 14
    assert abs(opc["sddb"]["sensitivity"] - 6 / 11) < 1e-12
    assert abs(opc["vfdb"]["sensitivity"] - 8 / 22) < 1e-12
    assert abs(opc["all_events"]["sensitivity"] - 14 / 33) < 1e-12

    assert sdd["sddb"]["n_detected"] == 11
    assert sdd["vfdb"]["n_detected"] == 21
    assert sdd["all_events"]["n_detected"] == 32
    assert abs(sdd["all_events"]["sensitivity"] - 32 / 33) < 1e-12

    assert absz["sddb"]["n_detected"] == 11
    assert absz["vfdb"]["n_detected"] == 19
    assert absz["all_events"]["n_detected"] == 30
    assert abs(absz["all_events"]["sensitivity"] - 30 / 33) < 1e-12


def test_extract_far_matches_json(far_summary):
    opc = extract_far_block(far_summary, "opc_L50")
    sdd = extract_far_block(far_summary, "sdd")
    absz = extract_far_block(far_summary, "absz_tau_s")

    assert opc["n_controls"] == 18
    assert opc["total_episodes"] == 28.0
    assert abs(opc["far_per_24h"] - 3.7333769938233505) < 1e-9
    assert sdd["total_episodes"] == 347.0
    assert absz["total_episodes"] == 253.0
    # FAR ranking: OPC < abs-z < SDD
    assert opc["far_per_24h"] < absz["far_per_24h"] < sdd["far_per_24h"]


def test_episode_sums_match_pooled(per_record_rows, far_summary):
    mapping = {
        "opc_L50": "opc_n_episodes",
        "sdd": "sdd_n_episodes",
        "absz_tau_s": "absz_n_episodes",
    }
    for far_key, col in mapping.items():
        summed = sum_episodes_from_per_record(per_record_rows, col)
        pooled = float(far_summary["pooled_far"][far_key]["total_episodes"])
        assert abs(summed - pooled) < 1e-6, f"{far_key}: {summed} != {pooled}"


def test_build_tradeoff_rows_join(exploratory, far_summary):
    rows = build_tradeoff_rows(exploratory, far_summary)
    assert len(rows) == 3
    by_id = {r["detector"]: r for r in rows}

    # Cross-check join fidelity
    assert abs(by_id["opc_L50"]["sens_all_events"] - 14 / 33) < 1e-12
    assert by_id["opc_L50"]["total_control_episodes"] == 28.0
    assert abs(by_id["opc_L50"]["far_per_24h"] - 3.7333769938233505) < 1e-9

    assert abs(by_id["sdd"]["sens_all_events"] - 32 / 33) < 1e-12
    assert by_id["sdd"]["total_control_episodes"] == 347.0

    assert abs(by_id["absz_tau_s"]["sens_all_events"] - 30 / 33) < 1e-12
    assert by_id["absz_tau_s"]["total_control_episodes"] == 253.0

    # No fusion
    assert all(r["fusion"] is False for r in rows)
    assert all(r["primary_arm"] is True for r in rows)

    # Balance fields present
    for r in rows:
        assert "sens_per_far_unit" in r
        assert "geometric_balance_sens_x_spec_proxy" in r
        assert r["sens_per_far_unit"] > 0


def test_cohort_summary_has_sddb_vfdb_nsrdb(exploratory, far_summary):
    rows = build_tradeoff_rows(exploratory, far_summary)
    cohort = build_cohort_summary_rows(rows)
    # 3 detectors × 4 cohort rows (SDDB, VFDB, all_events, NSRDB)
    assert len(cohort) == 12
    cohorts = {c["cohort"] for c in cohort}
    assert cohorts == {"SDDB", "VFDB", "all_events", "NSRDB"}
    nsrdb = [c for c in cohort if c["cohort"] == "NSRDB"]
    assert all(c["role"] == "controls_far" for c in nsrdb)
    assert all(c["sensitivity"] is None for c in nsrdb)
    assert all(c["far_per_24h"] is not None for c in nsrdb)


def test_balance_metrics_unit():
    b = balance_metrics(sens_all=0.5, far_per_24h=10.0, sens_sddb=0.6, sens_vfdb=0.4)
    assert abs(b["sens_per_far_unit"] - 0.05) < 1e-12
    assert abs(b["far_per_sens_unit"] - 20.0) < 1e-12
    assert 0 < b["specificity_proxy_1_over_1_plus_far24"] < 1
    assert b["sens_sddb_minus_vfdb_gap"] == pytest.approx(0.2)


def test_rank_detectors_lowest_far(exploratory, far_summary):
    rows = build_tradeoff_rows(exploratory, far_summary)
    ranking = rank_detectors(rows, "far_per_24h", higher_is_better=False)
    assert ranking[0]["detector"] == "opc_L50"
    assert ranking[1]["detector"] == "absz_tau_s"
    assert ranking[2]["detector"] == "sdd"


def test_rank_detectors_highest_sens(exploratory, far_summary):
    rows = build_tradeoff_rows(exploratory, far_summary)
    ranking = rank_detectors(rows, "sens_all_events", higher_is_better=True)
    assert ranking[0]["detector"] == "sdd"
    assert ranking[1]["detector"] == "absz_tau_s"
    assert ranking[2]["detector"] == "opc_L50"


def test_qualitative_notes_no_forbidden_claims(exploratory, far_summary):
    rows = build_tradeoff_rows(exploratory, far_summary)
    narrative = qualitative_notes(rows)
    blob = json.dumps(narrative).lower()
    for phrase in FORBIDDEN_CLAIM_PHRASES:
        assert phrase not in blob, f"Forbidden claim phrase present: {phrase}"
    assert "executive_summary" in narrative
    assert "recommendations" in narrative
    assert len(narrative["recommendations"]) >= 4
    assert "strengths_weaknesses" in narrative
    for arm in ("opc_L50", "sdd", "absz_tau_s"):
        assert arm in narrative["strengths_weaknesses"]
        assert narrative["strengths_weaknesses"][arm]["strengths"]
        assert narrative["strengths_weaknesses"][arm]["weaknesses"]
    # Must state exploratory / no clinical superiority framing in caveats
    caveats_blob = " ".join(narrative["caveats"]).lower()
    assert "clinical" in caveats_blob or "no clinical" in narrative["executive_summary"].lower()
    assert "s5" in caveats_blob


def test_run_analysis_writes_artifacts(tmp_path):
    """Drive real run_analysis; copy inputs to temp so we don't depend on cwd pollution."""
    import shutil

    for name in (
        "ordinal_exploratory_summary.json",
        "ordinal_nsrdb_far_summary.json",
        "ordinal_nsrdb_far_per_record.csv",
        "ordinal_opc_sdd_absz_comparison.csv",
    ):
        src = RES / name
        if not src.is_file():
            pytest.skip(f"missing {src}")
        shutil.copy(src, tmp_path / name)

    summary = run_analysis(results_dir=tmp_path, write_doc=False)

    tradeoff_csv = tmp_path / "ordinal_sensitivity_specificity_tradeoff.csv"
    cohort_csv = tmp_path / "ordinal_tradeoff_by_cohort.csv"
    summary_json = tmp_path / "ordinal_tradeoff_summary.json"

    assert tradeoff_csv.is_file()
    assert cohort_csv.is_file()
    assert summary_json.is_file()

    with open(tradeoff_csv, encoding="utf-8", newline="") as f:
        trows = list(csv.DictReader(f))
    assert len(trows) == 3
    detectors = {r["detector"] for r in trows}
    assert detectors == {"opc_L50", "sdd", "absz_tau_s"}

    # Recompute sens from source JSON and match CSV
    explor = load_json(tmp_path / "ordinal_exploratory_summary.json")
    far = load_json(tmp_path / "ordinal_nsrdb_far_summary.json")
    by_csv = {r["detector"]: r for r in trows}

    # CSV uses %.10g; allow tiny float serialization error
    atol = 1e-8
    assert abs(float(by_csv["opc_L50"]["sens_sddb"]) - explor["opc_L50_companion"]["sddb"]["sensitivity"]) < atol
    assert abs(float(by_csv["opc_L50"]["sens_vfdb"]) - explor["opc_L50_companion"]["vfdb"]["sensitivity"]) < atol
    assert abs(float(by_csv["opc_L50"]["sens_all_events"]) - explor["opc_L50_companion"]["all_events"]["sensitivity"]) < atol
    assert abs(float(by_csv["opc_L50"]["far_per_24h"]) - far["pooled_far"]["opc_L50"]["far_per_24h"]) < atol

    assert abs(float(by_csv["sdd"]["sens_all_events"]) - explor["sdd"]["all_events"]["sensitivity"]) < atol
    assert abs(float(by_csv["sdd"]["far_per_24h"]) - far["pooled_far"]["sdd"]["far_per_24h"]) < atol

    assert abs(float(by_csv["absz_tau_s"]["sens_all_events"]) - explor["absz_tau_s"]["all_events"]["sensitivity"]) < atol
    assert abs(float(by_csv["absz_tau_s"]["far_per_24h"]) - far["pooled_far"]["absz_tau_s"]["far_per_24h"]) < atol

    # Summary JSON flags
    assert summary["clinical_claim"] is False
    assert summary["superiority_claim"] is False
    assert summary["s5_claim"] is False
    assert summary["fusion"] is False
    assert summary["methodology"]["reprocessed_wfdb"] is False
    assert summary["exploratory_only"] is True
    assert all(summary["episode_crosscheck_per_record_vs_pooled"][k]["match"] for k in ("opc_L50", "sdd", "absz_tau_s"))
    assert summary["forbidden_claim_scan"]["hits"] == []
    assert "executive_summary" in summary["narrative"]
    assert len(summary["narrative"]["recommendations"]) >= 4


def test_project_results_after_main_run_if_present():
    """If project results already contain trade-off outputs, validate them."""
    tradeoff_csv = RES / "ordinal_sensitivity_specificity_tradeoff.csv"
    summary_json = RES / "ordinal_tradeoff_summary.json"
    if not tradeoff_csv.is_file() or not summary_json.is_file():
        pytest.skip("project trade-off outputs not yet generated")

    with open(tradeoff_csv, encoding="utf-8", newline="") as f:
        rows = list(csv.DictReader(f))
    assert len(rows) == 3
    with open(summary_json, encoding="utf-8") as f:
        summary = json.load(f)
    assert summary["fusion"] is False
    assert summary["clinical_claim"] is False
    assert summary["s5_claim"] is False
    for r in rows:
        assert r["detector"] in DISPLAY_NAMES
        assert r.get("fusion", "False") in ("False", "false", "0", "")


def test_no_fusion_in_source_module():
    """Structural: runner module must not define fused alarm combination."""
    src = (CODE / "run_ordinal_tradeoff_analysis.py").read_text(encoding="utf-8")
    assert "PRIMARY_ARMS" in src
    # Must not implement a fused primary arm
    assert "opc_and_sdd" not in src.lower()
    assert "fused_alarm" not in src.lower()
    assert "no fusion" in src.lower() or "no_fusion" in src.lower() or "fusion: false" in src.lower() or '"fusion": false' in src.lower() or "fusion" in src
