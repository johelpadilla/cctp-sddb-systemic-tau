"""Shared import bootstrap for CCTP paper code (standalone + monorepo)."""
from __future__ import annotations

import os
import sys
from pathlib import Path

CODE_DIR = Path(__file__).resolve().parent
BASE_DIR = CODE_DIR.parent


def ensure_code_on_path() -> None:
    p = str(CODE_DIR)
    if p not in sys.path:
        sys.path.insert(0, p)


def import_systemictau_core():
    """Import systemictau.core (pip package first, then common local sources)."""
    env = os.environ.get("SYSTEMICTAU_SRC")
    candidates = []
    if env:
        candidates.append(env)
    # Common monorepo / editable layouts relative to this paper package
    candidates.extend(
        [
            CODE_DIR / "vendor" / "systemictau" / "src",
            BASE_DIR.parent.parent / "systemictau" / "src",  # grok-safe/systemictau
            BASE_DIR.parent.parent / "Gemini" / "systemictau" / "src",
            Path.home() / "grok-safe" / "Gemini" / "systemictau" / "src",
            Path.home() / "grok-safe" / "systemictau" / "src",
        ]
    )
    try:
        from systemictau.core import compute_taus, systemic_tau  # type: ignore

        return compute_taus, systemic_tau, True
    except Exception:
        pass
    for c in candidates:
        c = Path(c)
        if c.is_dir():
            s = str(c)
            if s not in sys.path:
                sys.path.insert(0, s)
            try:
                from systemictau.core import compute_taus, systemic_tau  # type: ignore

                return compute_taus, systemic_tau, True
            except Exception:
                continue
    return None, None, False


def import_recd_ordinal_levels():
    """Import nested RECD level functions (vendored next to this file)."""
    ensure_code_on_path()
    # Legacy sibling path used during pilot development
    legacy = BASE_DIR.parent / "Conversacion_Naturaleza_Tiempo"
    if legacy.is_dir():
        s = str(legacy)
        if s not in sys.path:
            sys.path.append(s)
    from recd_ordinal_levels import (  # type: ignore
        generate_multivariate_symbols,
        compute_phi1,
        compute_phi2,
        compute_phi3,
        high_level3_rate,
        compute_recd_from_conjunctions,
        compute_weighted_contributions,
        compute_lambda,
    )

    return {
        "generate_multivariate_symbols": generate_multivariate_symbols,
        "compute_phi1": compute_phi1,
        "compute_phi2": compute_phi2,
        "compute_phi3": compute_phi3,
        "high_level3_rate": high_level3_rate,
        "compute_recd_from_conjunctions": compute_recd_from_conjunctions,
        "compute_weighted_contributions": compute_weighted_contributions,
        "compute_lambda": compute_lambda,
    }
