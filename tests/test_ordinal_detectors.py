#!/usr/bin/env python3
"""
Unit tests for native ordinal detectors (Option 1 OPC, Option 2 SDD).

Drives the real shipped functions in code/ordinal_detectors/.
Detectors are tested SEPARATELY — no fusion assertions.
No dependence on abs-z / continuous mean-variance alarm paths.
"""
from __future__ import annotations

import ast
import inspect
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
CODE = ROOT / "code"
sys.path.insert(0, str(CODE))

from ordinal_detectors.opc_detector import (  # noqa: E402
    opc_detect,
    opc_first_alarm_index,
)
from ordinal_detectors.sdd_detector import (  # noqa: E402
    empirical_distribution,
    kl_divergence_smoothed,
    sdd_detect,
    sdd_first_alarm_index,
    total_variation,
)


# ---------------------------------------------------------------------------
# Option 1 — Ordinal Persistence Collapse
# ---------------------------------------------------------------------------


class TestOPC:
    def test_collapse_then_persist_alarms(self):
        """Diverse stream then sustained single-symbol lock → OPC alarm."""
        K = 6
        L = 8
        theta_D = 0.35  # need |supp| <= 2 for K=6
        theta_R = 5
        # 40 mixed symbols (high diversity)
        rng = np.random.default_rng(0)
        head = rng.integers(0, K, size=40)
        # long collapse to one symbol
        tail = np.zeros(60, dtype=int)
        sigma = np.concatenate([head, tail])

        out = opc_detect(sigma, L=L, theta_D=theta_D, theta_R=theta_R, K=K)
        assert out["alarm"].sum() > 0
        # first alarm only after collapse has been sustained
        t_star, _ = opc_first_alarm_index(
            sigma, L=L, theta_D=theta_D, theta_R=theta_R, K=K, search_start=0
        )
        assert t_star is not None
        assert t_star >= 40  # not during the diverse head
        assert out["diversity"][t_star] <= theta_D
        assert out["persistence"][t_star] >= theta_R

    def test_high_diversity_no_alarm(self):
        """Uniform-ish random symbols over full alphabet → no OPC alarm."""
        K = 6
        rng = np.random.default_rng(1)
        sigma = rng.integers(0, K, size=200)
        out = opc_detect(sigma, L=8, theta_D=0.35, theta_R=5, K=K)
        # Extremely unlikely to stay at <=2 unique symbols for 5 consecutive
        # windows of length 8 under uniform draws; assert zero for this seed.
        assert out["alarm"].sum() == 0

    def test_diversity_is_support_over_K(self):
        sigma = np.array([0, 0, 0, 0, 1, 1, 1, 1], dtype=int)
        out = opc_detect(sigma, L=8, theta_D=0.5, theta_R=1, K=6)
        # At t=7, window is full: 2 unique / 6
        assert out["diversity"][7] == pytest.approx(2.0 / 6.0)

    def test_no_mean_variance_in_source(self):
        """Static: OPC code (AST) must not call continuous z / moment helpers."""
        src = Path(CODE / "ordinal_detectors" / "opc_detector.py").read_text()
        tree = ast.parse(src)
        # No imports from cctp_metrics_core / lead-time path
        for node in ast.walk(tree):
            if isinstance(node, (ast.Import, ast.ImportFrom)):
                mod = getattr(node, "module", None) or ""
                names = [a.name for a in getattr(node, "names", [])]
                blob = " ".join([mod] + names)
                assert "cctp_metrics" not in blob
                assert "detect_lead_time" not in blob
                assert "basal_stats" not in blob
        # No np.mean / np.std / np.var attribute calls in code
        for node in ast.walk(tree):
            if isinstance(node, ast.Attribute) and isinstance(node.value, ast.Name):
                if node.value.id == "np" and node.attr in {"mean", "std", "var"}:
                    pytest.fail(f"OPC must not call np.{node.attr}")

    def test_opc_independent_of_sdd_module(self):
        src = Path(CODE / "ordinal_detectors" / "opc_detector.py").read_text()
        assert "sdd" not in src.lower()
        assert "total_variation" not in src


# ---------------------------------------------------------------------------
# Option 2 — Symbolic Distribution Divergence
# ---------------------------------------------------------------------------


class TestSDD:
    def test_tv_zero_on_identical(self):
        p = np.array([0.5, 0.5, 0.0])
        assert total_variation(p, p) == pytest.approx(0.0)

    def test_tv_one_on_disjoint(self):
        p = np.array([1.0, 0.0, 0.0])
        q = np.array([0.0, 1.0, 0.0])
        assert total_variation(p, q) == pytest.approx(1.0)

    def test_reorganization_alarms(self):
        """Basal uses symbol 0; later regime uses symbol 1 → high TV alarm."""
        K = 4
        L_c = 20
        basal_len = 80
        # basal: almost all 0
        basal = np.zeros(basal_len, dtype=int)
        # post-basal: almost all 1
        post = np.ones(100, dtype=int)
        sigma = np.concatenate([basal, post])
        basal_range = (0, basal_len)

        out = sdd_detect(
            sigma,
            basal_range,
            L_c=L_c,
            theta_TV=0.35,
            theta_S=1,
            K=K,
            mask_basal=True,
        )
        assert out["alarm"].sum() > 0
        t_star, _ = sdd_first_alarm_index(
            sigma,
            basal_range,
            L_c=L_c,
            theta_TV=0.35,
            K=K,
            search_start=basal_len,
        )
        assert t_star is not None
        assert t_star >= basal_len + L_c - 1
        assert out["tv"][t_star] >= 0.35

    def test_same_distribution_no_alarm(self):
        """Stationary i.i.d. stream: basal and current match → no/rare TV alarm."""
        K = 4
        rng = np.random.default_rng(2)
        sigma = rng.integers(0, K, size=400)
        basal_range = (0, 100)
        out = sdd_detect(
            sigma,
            basal_range,
            L_c=50,
            theta_TV=0.45,  # moderately strict
            theta_S=3,
            K=K,
            mask_basal=True,
        )
        # With TV>=0.45 sustained 3 times on i.i.d. matching basal, expect zero
        assert out["alarm"].sum() == 0

    def test_p_basal_is_empirical(self):
        sigma = np.array([0, 0, 0, 1, 1, 2], dtype=int)
        p = empirical_distribution(sigma, 0, 6, K=3)
        assert p.sum() == pytest.approx(1.0)
        assert p[0] == pytest.approx(3 / 6)
        assert p[1] == pytest.approx(2 / 6)
        assert p[2] == pytest.approx(1 / 6)

    def test_kl_secondary_exists_but_not_in_alarm_path(self):
        """KL helper exists for diagnostics; sdd_detect source uses TV only."""
        p = np.array([0.5, 0.5])
        q = np.array([0.5, 0.5])
        assert kl_divergence_smoothed(p, q) == pytest.approx(0.0, abs=1e-8)
        src = Path(CODE / "ordinal_detectors" / "sdd_detector.py").read_text()
        # alarm path function body of sdd_detect must call total_variation
        assert "total_variation" in src
        tree = ast.parse(src)
        sdd_fn = None
        for node in tree.body:
            if isinstance(node, ast.FunctionDef) and node.name == "sdd_detect":
                sdd_fn = node
                break
        assert sdd_fn is not None
        called = {
            n.func.id
            for n in ast.walk(sdd_fn)
            if isinstance(n, ast.Call) and isinstance(n.func, ast.Name)
        }
        assert "total_variation" in called
        assert "kl_divergence_smoothed" not in called

    def test_no_mean_variance_alarm_on_continuous(self):
        """Static: SDD code (AST) must not call continuous z / moment helpers."""
        src = Path(CODE / "ordinal_detectors" / "sdd_detector.py").read_text()
        tree = ast.parse(src)
        for node in ast.walk(tree):
            if isinstance(node, (ast.Import, ast.ImportFrom)):
                mod = getattr(node, "module", None) or ""
                names = [a.name for a in getattr(node, "names", [])]
                blob = " ".join([mod] + names)
                assert "cctp_metrics" not in blob
                assert "detect_lead_time" not in blob
                assert "basal_stats" not in blob
        for node in ast.walk(tree):
            if isinstance(node, ast.Attribute) and isinstance(node.value, ast.Name):
                if node.value.id == "np" and node.attr in {"mean", "std", "var"}:
                    pytest.fail(f"SDD must not call np.{node.attr}")

    def test_sdd_independent_of_opc_module(self):
        src = Path(CODE / "ordinal_detectors" / "sdd_detector.py").read_text()
        assert "opc_detector" not in src
        assert "opc_detect" not in src
        assert "from .opc" not in src
        assert "import opc" not in src


# ---------------------------------------------------------------------------
# Package / formalization structure
# ---------------------------------------------------------------------------


class TestSeparationAndDocs:
    def test_init_exports_both_separately_no_fusion(self):
        import ordinal_detectors as od

        assert hasattr(od, "opc_detect")
        assert hasattr(od, "sdd_detect")
        assert not hasattr(od, "fused_detect")
        assert not hasattr(od, "combined_alarm")
        assert not hasattr(od, "joint_alarm")
        src = (CODE / "ordinal_detectors" / "__init__.py").read_text().lower()
        assert "independent" in src or "separately" in src or "no fused" in src
        # No callable that combines both
        assert "def fuse" not in src
        assert "def combined" not in src

    def test_formalization_doc_exists_and_ordered(self):
        doc = ROOT / "docs" / "ORDINAL_ALARM_DETECTORS.md"
        assert doc.is_file()
        text = doc.read_text()
        i1 = text.find("# Option 1")
        i2 = text.find("# Option 2")
        assert i1 >= 0 and i2 >= 0 and i1 < i2
        assert "Ordinal Persistence Collapse" in text
        assert "Symbolic Distribution Divergence" in text
        assert "Total Variation" in text
        assert r"\theta_D" in text or "theta_D" in text
        assert "0.25" in text and "0.45" in text  # diversity / TV ranges spirit
        # persistence range spirit from Contexto
        assert "4" in text and "6" in text
        # fusion explicitly rejected as deliverable (may name the forbidden form)
        assert "no fusion" in text.lower() or "No rule of the form" in text
        assert "Explicit non-goals" in text or "completely independent" in text.lower()
        # contrast with abs-z
        assert "abs-z" in text.lower() or "z-score" in text.lower()
        # basal moments of continuous series rejected
        assert "mean" in text.lower() and "variance" in text.lower()

    def test_opc_source_inspectable_signature(self):
        sig = inspect.signature(opc_detect)
        assert "theta_D" in sig.parameters
        assert "theta_R" in sig.parameters
        assert "L" in sig.parameters

    def test_sdd_source_inspectable_signature(self):
        sig = inspect.signature(sdd_detect)
        assert "theta_TV" in sig.parameters
        assert "L_c" in sig.parameters
        assert "basal" in sig.parameters
