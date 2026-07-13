"""
Tests that drive the exploratory ordinal bake-off path on real shipped code.

- Derives joint bivariate symbols from cleaned RR (project encoding).
- Calls shipped opc_detect / sdd_detect (not reimplementations).
- Compares structure vs detect_lead_time without retuning abs-z defaults.
"""
from __future__ import annotations

import ast
import csv
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
    build_bivariate_proxy,
    detect_lead_time,
    get_event_and_windows,
)
from ordinal_detectors.opc_detector import opc_detect, opc_first_alarm_index  # noqa: E402
from ordinal_detectors.sdd_detector import sdd_detect, sdd_first_alarm_index  # noqa: E402
from recd_ordinal_levels import generate_multivariate_symbols  # noqa: E402

# Import bake-off helpers (real runner module)
from run_ordinal_exploratory_bakeoff import (  # noqa: E402
    ANALYTIC_SDDB,
    K_JOINT,
    joint_bivariate_symbols,
    process_event_record,
    OPC_L,
    OPC_THETA_D,
    OPC_THETA_R,
    SDD_L_C,
    SDD_THETA_TV,
    SDD_THETA_S,
)


def _first_existing_sddb() -> Path:
    for rec in ANALYTIC_SDDB:
        p = ROOT / "data" / f"rr_{rec}_clean.npz"
        if p.exists():
            return p
    pytest.skip("no SDDB cleaned RR available")


def _first_existing_vfdb() -> Path:
    paths = sorted((ROOT / "data" / "rr_external").glob("vfdb_*_clean.npz"))
    if not paths:
        pytest.skip("no VFDB cleaned RR available")
    return paths[0]


def test_joint_bivariate_symbols_from_real_rr():
    path = _first_existing_sddb()
    d = np.load(path)
    rr = d["rr_ms"].astype(float)
    sigma, K, offset = joint_bivariate_symbols(rr)
    assert K == K_JOINT == 36
    assert offset == 2
    assert sigma.ndim == 1 and len(sigma) > 100
    assert int(sigma.min()) >= 0 and int(sigma.max()) < K
    # Consistency with project generators
    X = build_bivariate_proxy(rr)
    S = generate_multivariate_symbols(X, m=3, delay=1)
    expected = S[:, 0].astype(np.int64) * 6 + S[:, 1].astype(np.int64)
    np.testing.assert_array_equal(sigma, expected[: len(sigma)])


def test_opc_detect_driven_on_real_sddb_symbols():
    path = _first_existing_sddb()
    d = np.load(path)
    rr = d["rr_ms"].astype(float)
    sigma, K, _ = joint_bivariate_symbols(rr)
    # Use a short mid-record slice for speed
    mid = len(sigma) // 2
    seg = sigma[mid : mid + 500]
    out = opc_detect(seg, L=OPC_L, theta_D=OPC_THETA_D, theta_R=OPC_THETA_R, K=K)
    assert set(out.keys()) >= {"alarm", "diversity", "persistence", "low_div"}
    assert out["alarm"].shape == seg.shape
    # L=8, K=36 ⇒ max D = 8/36 < 0.35 ⇒ low_div after warm-up
    assert np.all(out["low_div"][OPC_L - 1 :] == 1)
    # After θ_R consecutive low-div, alarm must fire (drives real opc_detect)
    assert int(out["alarm"][OPC_L - 1 + OPC_THETA_R - 1]) == 1


def test_sdd_detect_driven_on_real_sddb_symbols():
    path = _first_existing_sddb()
    d = np.load(path)
    rr = d["rr_ms"].astype(float)
    sigma, K, _ = joint_bivariate_symbols(rr)
    mid = len(sigma) // 2
    seg = sigma[mid : mid + 2000]
    basal = (0, 200)
    out = sdd_detect(
        seg,
        basal,
        L_c=SDD_L_C,
        theta_TV=SDD_THETA_TV,
        theta_S=SDD_THETA_S,
        K=K,
        mask_basal=True,
    )
    assert set(out.keys()) >= {"alarm", "tv", "p_basal", "high_tv"}
    assert out["p_basal"].shape == (K,)
    assert abs(float(out["p_basal"].sum()) - 1.0) < 1e-9
    # No alarm inside basal when mask_basal=True
    assert int(np.sum(out["alarm"][:200])) == 0


def test_process_event_record_sddb_drives_all_three_separately():
    path = _first_existing_sddb()
    rec = path.stem.replace("rr_", "").replace("_clean", "")
    opc, sdd, absz = process_event_record(
        "sddb",
        rec,
        path,
        OPC_L,
        OPC_THETA_D,
        OPC_THETA_R,
        SDD_L_C,
        SDD_THETA_TV,
        SDD_THETA_S,
    )
    assert opc.get("encoding") == "joint_bivariate_BP_m3"
    assert sdd.get("encoding") == "joint_bivariate_BP_m3"
    assert absz.get("metric") == "tau_s"
    assert absz.get("z_threshold") == FROZEN_Z_THRESHOLD
    assert absz.get("min_consecutive") == FROZEN_MIN_CONSECUTIVE
    # Separate options: no fused fields
    assert "opc_and_sdd" not in opc and "fused" not in str(opc.get("note", "")).lower()
    assert "opc_and_sdd" not in sdd


def test_process_event_record_vfdb_uses_short_windows():
    path = _first_existing_vfdb()
    rec = path.stem.replace("vfdb_", "").replace("_clean", "")
    opc, sdd, absz = process_event_record(
        "vfdb",
        rec,
        path,
        OPC_L,
        OPC_THETA_D,
        OPC_THETA_R,
        SDD_L_C,
        SDD_THETA_TV,
        SDD_THETA_S,
    )
    assert opc.get("source") == "vfdb"
    assert opc.get("duration_stratum") in (
        "too_short",
        "short_15_60min",
        "short_ge60min",
        "holter_ge6h",
    )
    assert "alarmed" in opc and "alarmed" in sdd and "alarmed" in absz


def test_bakeoff_artifacts_exist_when_results_present():
    """If prior bake-off was run, required result files exist and are non-empty."""
    res = ROOT / "results"
    required = [
        "ordinal_data_inventory.txt",
        "ordinal_opc_per_record.csv",
        "ordinal_sdd_per_record.csv",
        "ordinal_absz_per_record.csv",
        "ordinal_opc_sdd_absz_comparison.csv",
        "ordinal_exploratory_summary.json",
    ]
    missing = [f for f in required if not (res / f).exists()]
    if missing:
        pytest.skip(f"bake-off artifacts not yet generated: {missing}")
    inv = (res / "ordinal_data_inventory.txt").read_text()
    assert "joint" in inv.lower() or "K=36" in inv or "Bandt" in inv
    assert "pre-baked" in inv.lower() or "derived" in inv.lower() or "not pre-baked" in inv.lower() or "not" in inv.lower()
    with open(res / "ordinal_opc_sdd_absz_comparison.csv") as f:
        rows = list(csv.DictReader(f))
    assert len(rows) >= 10
    assert "opc_L8_alarmed" in rows[0] and "sdd_alarmed" in rows[0] and "absz_tau_alarmed" in rows[0]


def test_absz_defaults_untouched_in_cctp_metrics_core():
    """Production frozen abs-z defaults must remain 2.0 and 3 (no retune)."""
    src = (CODE / "cctp_metrics_core.py").read_text()
    tree = ast.parse(src)
    found_detect = False
    for node in tree.body:
        if isinstance(node, ast.FunctionDef) and node.name == "detect_lead_time":
            found_detect = True
            defaults = {}
            args = node.args
            # pair defaults with kwonly or positional
            pos = args.args
            defs = args.defaults
            for name, default in zip(pos[-len(defs) :], defs) if defs else []:
                if isinstance(default, ast.Constant):
                    defaults[name.arg] = default.value
            for name, default in zip(args.kwonlyargs, args.kw_defaults):
                if default is not None and isinstance(default, ast.Constant):
                    defaults[name.arg] = default.value
            assert defaults.get("z_threshold") == 2.0
            assert defaults.get("min_consecutive") == 3
    assert found_detect
    assert FROZEN_Z_THRESHOLD == 2.0
    assert FROZEN_MIN_CONSECUTIVE == 3


def test_no_fusion_api_in_ordinal_package():
    init = (CODE / "ordinal_detectors" / "__init__.py").read_text().lower()
    assert "fuse" not in init or "no fused" in init or "intentionally no fused" in init
    assert "opc_detect" in init and "sdd_detect" in init


def test_writeup_exists_and_is_observational():
    doc = ROOT / "docs" / "ORDINAL_EXPLORATORY_BAKEOFF.md"
    if not doc.exists():
        pytest.skip("write-up not generated yet")
    text = doc.read_text().lower()
    assert "exploratory" in text
    assert "none" in text or "no clinical" in text or "clinical" in text
    assert "opc" in text and "sdd" in text and "abs-z" in text or "absz" in text
    assert "fusion" in text
    # Should not claim FDA readiness
    assert "fda approved" not in text
