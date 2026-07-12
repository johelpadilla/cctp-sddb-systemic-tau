"""Structural checks for Phase 2 planning artifacts + frozen primary params.

These tests drive the real in-repo planning deliverables and the shipped
frozen constants in cctp_metrics_core (no threshold retune).
"""
from __future__ import annotations

import csv
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


class TestFrozenPrimaryParams(unittest.TestCase):
    """Primary rule must remain discovery/Phase-1 frozen values."""

    def test_frozen_constants(self):
        self.assertEqual(FROZEN_THETA3, 0.08)
        self.assertEqual(FROZEN_HIGH_THRESHOLD, 0.65)
        self.assertEqual(W_TAU, 101)
        self.assertEqual(FROZEN_Z_THRESHOLD, 2.0)
        self.assertEqual(FROZEN_MIN_CONSECUTIVE, 3)


class TestPhase2PlanDocument(unittest.TestCase):
    def setUp(self):
        self.path = ROOT / "docs" / "EXTERNAL_VALIDATION_PHASE2_PLAN.md"
        self.text = self.path.read_text(encoding="utf-8")

    def test_plan_exists(self):
        self.assertTrue(self.path.is_file(), f"missing {self.path}")

    def test_scientific_goal_is_specificity_far(self):
        low = self.text.lower()
        self.assertIn("specificity", low)
        self.assertRegex(self.text, re.compile(r"\bFAR\b|false-alarm", re.I))
        # Explicitly frames Phase 2 as NOT more short-DB sensitivity hunting
        self.assertRegex(
            self.text,
            re.compile(r"more short-DB sensitivity hunting", re.I),
        )
        self.assertIn("phase 2 is **not**", low)

    def test_phase1_bottleneck_numbers(self):
        # Honest Phase 1 numbers that motivate Phase 2
        self.assertIn("n = 11", self.text.replace("n=11", "n = 11"))
        self.assertTrue("1.00" in self.text or "1.0" in self.text)
        self.assertTrue("0.82" in self.text or "≈ 0.82" in self.text or "~0.82" in self.text)
        self.assertTrue("34.4" in self.text or "~34" in self.text)
        self.assertTrue("28.8" in self.text or "~28" in self.text)
        self.assertIn("S5", self.text)
        self.assertRegex(self.text, re.compile(r"no clinical|NONE", re.I))

    def test_frozen_params_echoed(self):
        self.assertIn("0.08", self.text)
        self.assertIn("0.65", self.text)
        self.assertIn("101", self.text)
        self.assertRegex(self.text, re.compile(r"abs-z|abs.?z", re.I))
        self.assertIn("not retune", self.text.lower().replace("do not retune", "not retune"))

    def test_priority_sources_and_strategies(self):
        low = self.text.lower()
        self.assertIn("institutional", low)
        self.assertIn("device-matched", low)
        self.assertIn("irb", low)
        self.assertIn("exploratory", low)
        self.assertTrue("basal" in low or "refractory" in low or "fusion" in low)


class TestPhase2InventoryAndChecklist(unittest.TestCase):
    REQUIRED_INVENTORY_COLS = {
        "candidate_id",
        "source",
        "record_id",
        "tier",
        "role",
        "local_status",
        "device_match_to_vfdb",
        "include_for_primary_far",
        "exclusion_or_caveat",
        "phase1_used",
        "quality_notes",
        "next_action",
    }

    def test_inventory_csv_schema_and_quality_flags(self):
        path = ROOT / "results" / "phase2_public_control_inventory.csv"
        self.assertTrue(path.is_file(), f"missing {path}")
        with path.open(newline="", encoding="utf-8") as f:
            rows = list(csv.DictReader(f))
        self.assertGreaterEqual(len(rows), 6, "expect at least Phase-1 NSRDB controls inventoried")
        cols = set(rows[0].keys()) if rows else set()
        missing = self.REQUIRED_INVENTORY_COLS - cols
        self.assertFalse(missing, f"inventory missing columns: {missing}")

        # Quality-over-quantity: Phase-1 NSRDB present; device mismatch flagged
        phase1 = [r for r in rows if r.get("phase1_used", "").strip().lower() == "yes"]
        self.assertEqual(len(phase1), 6)
        for r in phase1:
            self.assertEqual(r["device_match_to_vfdb"].strip().lower(), "no")
            self.assertIn("mismatch", r["include_for_primary_far"].lower() + r["exclusion_or_caveat"].lower())

        # Institutional rows are planned only (not fabricated as available)
        inst = [r for r in rows if r.get("tier", "").startswith("A")]
        self.assertGreaterEqual(len(inst), 1)
        for r in inst:
            self.assertIn(
                r["local_status"],
                {"not_available", "planned", "not_downloaded"},
            )

        # No primary claim that S5 is met
        blob = path.read_text(encoding="utf-8").lower()
        self.assertNotIn("s5_met", blob)
        self.assertNotIn("clinical_claim=true", blob)

    def test_irb_checklist_exists(self):
        path = ROOT / "docs" / "PHASE2_IRB_DATA_CHECKLIST.md"
        self.assertTrue(path.is_file())
        text = path.read_text(encoding="utf-8").lower()
        for needle in ("irb", "de-identification", "device-matched", "frozen", "no clinical"):
            self.assertIn(needle, text)


class TestHandoffPhase2(unittest.TestCase):
    def test_handoff_mentions_phase2_and_next_prompt(self):
        path = ROOT / "HANDOFF.md"
        self.assertTrue(path.is_file())
        text = path.read_text(encoding="utf-8")
        self.assertIn("Phase 2", text)
        self.assertIn("EXTERNAL_VALIDATION_PHASE2_PLAN.md", text)
        self.assertIn("Next-session prompt", text)
        # Frozen params still present
        self.assertIn("0.08", text)
        self.assertIn("0.65", text)
        # Tier A path active after public interim
        self.assertIn("PHASE2_INSTITUTIONAL_DATA_REQUEST.md", text)
        self.assertIn("Tier A", text)


class TestPhase2InstitutionalDataRequest(unittest.TestCase):
    """Partner-facing Tier A request brief — sendable skeleton, no fabricated PHI."""

    def setUp(self):
        self.path = ROOT / "docs" / "PHASE2_INSTITUTIONAL_DATA_REQUEST.md"
        self.text = self.path.read_text(encoding="utf-8")
        self.low = self.text.lower()

    def test_request_exists(self):
        self.assertTrue(self.path.is_file(), f"missing {self.path}")

    def test_partner_archive_classes_not_fabricated_sites(self):
        # Generic Tier A classes present
        for needle in (
            "holter",
            "telemetry",
            "device-matched",
            "partner",
            "archive",
        ):
            self.assertIn(needle, self.low)
        # Must not invent real institutional PHI inventory IDs
        self.assertNotRegex(self.text, re.compile(r"\bMRN\s*[:=]\s*\d+", re.I))
        self.assertNotIn("clinical_claim: true", self.low)
        self.assertNotIn("s5_met", self.low.replace("not met", "x"))

    def test_control_inclusion_criteria(self):
        self.assertIn("adult", self.low)
        self.assertRegex(self.text, re.compile(r"≥\s*4|>=\s*4|ge.?4", re.I))
        self.assertRegex(self.text, re.compile(r"24\s*h|24h", re.I))
        self.assertIn("no", self.low)  # no sustained VT/VF language
        self.assertRegex(self.text, re.compile(r"sustained\s+VT|VT/VF", re.I))
        self.assertIn("device-matched", self.low)
        self.assertRegex(self.text, re.compile(r"10\s*[–-]\s*20|n\s*≈\s*10", re.I))

    def test_deidentification_scheme(self):
        for needle in (
            "de-identification",
            "research_id",
            "mrn",
            "relative",
            "rr",
            "annotation",
            "git",
        ):
            self.assertIn(needle, self.low)

    def test_optional_prevf_and_frozen_no_claim(self):
        self.assertRegex(self.text, re.compile(r"pre-?event|pre-?vf|spontaneous", re.I))
        self.assertIn("0.08", self.text)
        self.assertIn("0.65", self.text)
        self.assertIn("101", self.text)
        self.assertRegex(self.text, re.compile(r"abs-z|abs.?z", re.I))
        self.assertRegex(self.low, re.compile(r"not\s+retun|no\s+retun"))
        self.assertRegex(self.low, re.compile(r"no clinical|not claimed|none"))
        self.assertIn("decision support", self.low)

    def test_progress_records_tier_a_and_public_done(self):
        prog = (ROOT / "docs" / "EXTERNAL_VALIDATION_PHASE2_PROGRESS.md").read_text(
            encoding="utf-8"
        )
        self.assertIn("PHASE2_INSTITUTIONAL_DATA_REQUEST.md", prog)
        self.assertIn("n=18", prog.replace("n = 18", "n=18"))
        self.assertIn("Tier A", prog)
        low = prog.lower()
        self.assertIn("device mismatch", low)
        # Honest non-claim language required; forbid success-as-proven wording only
        self.assertTrue(
            "not met" in low or "not claimed" in low or "false" in low,
            "progress must deny S5/clinical success",
        )
        for bad in ("s5 passed", "s5 achieved", "s5 satisfied", "clinical_claim: true"):
            self.assertNotIn(bad, low)


if __name__ == "__main__":
    unittest.main()
