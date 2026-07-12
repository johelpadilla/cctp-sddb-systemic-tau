"""Tests for Phase 2 interim FAR artifacts: frozen params, no S5 claim, inventory.

Drives the real shipped result JSON / CSV produced by
`code/run_external_validation_phase2_far.py` and inventory updates.
"""
from __future__ import annotations

import csv
import json
import re
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))

from cctp_metrics_core import (  # noqa: E402
    FROZEN_HIGH_THRESHOLD,
    FROZEN_MIN_CONSECUTIVE,
    FROZEN_THETA3,
    FROZEN_Z_THRESHOLD,
    W_TAU,
)


class TestPhase2FarJson(unittest.TestCase):
    def setUp(self):
        self.path = ROOT / "results" / "external_phase2_far.json"
        self.assertTrue(self.path.is_file(), f"missing {self.path}")
        self.data = json.loads(self.path.read_text(encoding="utf-8"))

    def test_frozen_params_match_core(self):
        p = self.data["params"]
        self.assertEqual(p["theta3"], FROZEN_THETA3)
        self.assertEqual(p["high_threshold"], FROZEN_HIGH_THRESHOLD)
        self.assertEqual(p["W_TAU"], W_TAU)
        self.assertEqual(p["z_threshold"], FROZEN_Z_THRESHOLD)
        self.assertEqual(p["min_consecutive"], FROZEN_MIN_CONSECUTIVE)
        self.assertEqual(p["theta3"], 0.08)
        self.assertEqual(p["high_threshold"], 0.65)
        self.assertEqual(p["W_TAU"], 101)
        self.assertEqual(p["z_threshold"], 2.0)
        self.assertEqual(p["min_consecutive"], 3)

    def test_n_controls_expanded(self):
        tau = self.data["by_metric"]["tau_s"]
        self.assertGreaterEqual(float(tau["n_controls"]), 6.0)
        # Full NSRDB target is 18 when all process
        self.assertGreaterEqual(float(tau["n_controls"]), 12.0)
        self.assertTrue(float(tau["total_search_hours"]) > 0)
        self.assertTrue(float(tau["far_per_24h"]) == float(tau["far_per_24h"]))  # finite

    def test_device_mismatch_and_no_s5_claim(self):
        self.assertTrue(self.data.get("device_mismatch") is True)
        note = (self.data.get("device_mismatch_note") or "") + " ".join(
            self.data.get("notes") or []
        )
        self.assertRegex(note, re.compile(r"device.?match|not device-matched", re.I))
        self.assertFalse(self.data.get("s5_claim", True))
        self.assertFalse(self.data.get("clinical_claim", True))
        self.assertFalse(self.data.get("deployability_claim", True))
        joined = note.lower()
        self.assertIn("not claimed", joined.replace("is not claimed", "not claimed"))


class TestPhase2Summary(unittest.TestCase):
    def setUp(self):
        self.path = ROOT / "results" / "external_phase2_summary.json"
        self.assertTrue(self.path.is_file())
        self.data = json.loads(self.path.read_text(encoding="utf-8"))

    def test_recommendation_and_limitations(self):
        self.assertFalse(self.data["s5_claim"])
        self.assertFalse(self.data["clinical_claim"])
        self.assertTrue(self.data["device_mismatch"])
        rec = (self.data.get("recommendation") or "") + (self.data.get("next_step") or "")
        self.assertRegex(rec, re.compile(r"institutional|IRB|device-matched", re.I))
        # Must not claim S5 success (allow honest "not met" / "not claimed")
        self.assertFalse(self.data.get("s5_claim"))
        blob = json.dumps(self.data).lower()
        self.assertNotRegex(
            blob,
            re.compile(r"s5.{0,40}\b(passed|achieved|satisfied)\b"),
        )
        self.assertRegex(blob, re.compile(r"not (met|claimed)"))


class TestPhase2InventoryUpdated(unittest.TestCase):
    def test_no_nsrdb_left_not_downloaded_if_processed(self):
        path = ROOT / "results" / "phase2_public_control_inventory.csv"
        self.assertTrue(path.is_file())
        rows = list(csv.DictReader(path.open(encoding="utf-8")))
        nsrdb = [r for r in rows if "NSRDB" in (r.get("source") or "")]
        self.assertGreaterEqual(len(nsrdb), 18)
        still = [r for r in nsrdb if r.get("local_status") == "not_downloaded"]
        # After Phase 2 interim, all 18 should be present if download succeeded
        self.assertEqual(
            still,
            [],
            f"NSRDB still not_downloaded: {[r.get('record_id') for r in still]}",
        )
        present = [r for r in nsrdb if "present" in (r.get("local_status") or "")]
        self.assertEqual(len(present), 18)


class TestPhase2ProgressDoc(unittest.TestCase):
    def test_progress_has_far_numbers_and_recommendation(self):
        path = ROOT / "docs" / "EXTERNAL_VALIDATION_PHASE2_PROGRESS.md"
        self.assertTrue(path.is_file())
        text = path.read_text(encoding="utf-8")
        self.assertRegex(text, re.compile(r"33\.|~33"))
        self.assertRegex(text, re.compile(r"32\.|~32"))
        self.assertIn("device", text.lower())
        self.assertIn("institutional", text.lower())
        self.assertRegex(text, re.compile(r"S5", re.I))
        self.assertRegex(text, re.compile(r"not (met|claimed)", re.I))
        self.assertIn("0.08", text)
        self.assertIn("0.65", text)
        self.assertIn("101", text)


class TestPhase2RunnerUsesFrozenImports(unittest.TestCase):
    """Structural: Phase-2 runner imports Phase-1 helpers (no forked detector)."""

    def test_runner_reuses_phase1_analyze_control(self):
        src = (ROOT / "code" / "run_external_validation_phase2_far.py").read_text(
            encoding="utf-8"
        )
        self.assertIn("from run_external_validation_phase1 import", src)
        self.assertIn("analyze_control_record", src)
        self.assertIn("false_alarm_rate", src)
        self.assertIn("FROZEN_THETA3", src)
        self.assertIn("device_mismatch", src)
        self.assertIn("s5_claim", src)


if __name__ == "__main__":
    unittest.main()
