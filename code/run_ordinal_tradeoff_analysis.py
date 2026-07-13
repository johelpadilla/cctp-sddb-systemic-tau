#!/usr/bin/env python3
"""
Combined sensitivity–specificity (FAR) trade-off across three detectors.

Joins *existing* bake-off and NSRDB FAR artifacts — no reprocessing of .wfdb
or full RR pipelines.

Primary arms (kept strictly separate; no fusion):
  - OPC L=50 companion  (opc_L50_companion / opc_L50)
  - SDD (TV)            (sdd)
  - Frozen abs-z τ_s    (absz_tau_s)

OPC L=8 is *not* a primary arm (K/L-invalid saturation on joint K=36).

Observational / exploratory only — no clinical, S5, FDA, or superiority claims.

Outputs under results/:
  ordinal_sensitivity_specificity_tradeoff.csv
  ordinal_tradeoff_by_cohort.csv
  ordinal_tradeoff_summary.json
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

BASE = Path(__file__).resolve().parent.parent
RES = BASE / "results"
DOCS = BASE / "docs"

# Key map: display name → (sensitivity JSON key, FAR pooled key)
PRIMARY_ARMS: List[Tuple[str, str, str]] = [
    ("opc_L50", "opc_L50_companion", "opc_L50"),
    ("sdd", "sdd", "sdd"),
    ("absz_tau_s", "absz_tau_s", "absz_tau_s"),
]

DISPLAY_NAMES = {
    "opc_L50": "OPC L=50",
    "sdd": "SDD (TV)",
    "absz_tau_s": "abs-z τ_s (frozen)",
}

# Affirmative forbidden framing (denials like "no S5 claim" are allowed)
FORBIDDEN_CLAIM_PHRASES = (
    "clinical superiority",
    "clinically superior",
    "is clinically superior",
    "fda ready",
    "fda-ready",
    "s5 achieved",
    "achieves s5",
    "ready for deployment",
    "device-matched superiority",
    "proven superior",
    "clinically validated for deployment",
)


def safe_relpath(path: Path, base: Path = BASE) -> str:
    """Return path relative to base when possible; else absolute string."""
    try:
        return str(path.resolve().relative_to(base.resolve()))
    except ValueError:
        return str(path.resolve())


def load_json(path: Path) -> dict:
    with open(path, encoding="utf-8") as f:
        return json.load(f)


def load_csv_rows(path: Path) -> List[dict]:
    with open(path, encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f))


def extract_sensitivity_block(exploratory: dict, sens_key: str) -> dict:
    """Pull SDDB / VFDB / all_events sensitivity from exploratory summary."""
    if sens_key not in exploratory:
        raise KeyError(f"Missing sensitivity arm '{sens_key}' in exploratory summary")
    block = exploratory[sens_key]
    out: Dict[str, Any] = {}
    for cohort in ("sddb", "vfdb", "all_events"):
        if cohort not in block:
            raise KeyError(f"Missing cohort '{cohort}' under '{sens_key}'")
        c = block[cohort]
        out[cohort] = {
            "n": int(c["n"]),
            "n_detected": int(c["n_detected"]),
            "sensitivity": float(c["sensitivity"]),
            "median_lead_h": float(c["median_lead_h"])
            if c.get("median_lead_h") is not None
            else None,
            "mean_lead_h": float(c["mean_lead_h"])
            if c.get("mean_lead_h") is not None
            else None,
        }
    return out


def extract_far_block(far_summary: dict, far_key: str) -> dict:
    """Pull pooled FAR fields for one detector from NSRDB FAR summary."""
    pooled = far_summary.get("pooled_far") or {}
    if far_key not in pooled:
        raise KeyError(f"Missing FAR arm '{far_key}' in pooled_far")
    p = pooled[far_key]
    stats = (far_summary.get("per_record_far_stats") or {}).get(far_key, {})
    return {
        "n_controls": int(p["n_controls"]),
        "total_episodes": float(p["total_episodes"]),
        "total_search_hours": float(p["total_search_hours"]),
        "far_per_24h": float(p["far_per_24h"]),
        "fraction_alarmed": float(p["fraction_alarmed"]),
        "mean_far_per_record": float(stats["mean"]) if stats else None,
        "median_far_per_record": float(stats["median"]) if stats else None,
        "min_far_per_record": float(stats["min"]) if stats else None,
        "max_far_per_record": float(stats["max"]) if stats else None,
        "std_far_per_record": float(stats["std"]) if stats else None,
    }


def sum_episodes_from_per_record(
    per_record_rows: List[dict], column: str
) -> float:
    """Sum episode counts from per-record FAR CSV (sanity cross-check)."""
    total = 0.0
    for row in per_record_rows:
        if row.get("skipped", "").lower() in ("true", "1"):
            continue
        total += float(row[column])
    return total


def balance_metrics(
    sens_all: float,
    far_per_24h: float,
    sens_sddb: float,
    sens_vfdb: float,
) -> dict:
    """
    Exploratory balance scores (not clinical decision metrics).

    - sens_per_far_unit: sens_all / far_per_24h  (hits per unit FAR; higher = more efficient)
    - far_per_sens_unit: far_per_24h / sens_all  (FAR cost of unit sensitivity)
    - geometric_balance: sqrt(sens_all * specificity_proxy) where specificity_proxy
      is a crude transform of FAR into (0,1] via 1/(1+far/24) — interpretive only
    """
    if far_per_24h <= 0:
        sens_per_far = float("inf") if sens_all > 0 else 0.0
        far_per_sens = 0.0
    else:
        sens_per_far = sens_all / far_per_24h
        far_per_sens = far_per_24h / sens_all if sens_all > 0 else float("inf")

    # Soft map FAR → [0,1] "specificity-like" score for ranking only:
    # FAR=0 → 1; FAR large → ~0. Not equal to classical specificity.
    specificity_proxy = 1.0 / (1.0 + far_per_24h / 24.0)
    geometric_balance = (sens_all * specificity_proxy) ** 0.5 if sens_all >= 0 else 0.0
    # Youden-like exploratory: sens + specificity_proxy - 1
    youden_like = sens_all + specificity_proxy - 1.0

    return {
        "sens_per_far_unit": sens_per_far,
        "far_per_sens_unit": far_per_sens,
        "specificity_proxy_1_over_1_plus_far24": specificity_proxy,
        "geometric_balance_sens_x_spec_proxy": geometric_balance,
        "youden_like_exploratory": youden_like,
        "sens_sddb_minus_vfdb_gap": sens_sddb - sens_vfdb,
    }


def build_tradeoff_rows(
    exploratory: dict, far_summary: dict
) -> List[Dict[str, Any]]:
    """Cross sensitivity and FAR into one row per primary detector."""
    rows: List[Dict[str, Any]] = []
    for arm_id, sens_key, far_key in PRIMARY_ARMS:
        sens = extract_sensitivity_block(exploratory, sens_key)
        far = extract_far_block(far_summary, far_key)
        bal = balance_metrics(
            sens_all=sens["all_events"]["sensitivity"],
            far_per_24h=far["far_per_24h"],
            sens_sddb=sens["sddb"]["sensitivity"],
            sens_vfdb=sens["vfdb"]["sensitivity"],
        )
        row: Dict[str, Any] = {
            "detector": arm_id,
            "detector_display": DISPLAY_NAMES[arm_id],
            "primary_arm": True,
            "fusion": False,
            # Sensitivity
            "sens_sddb": sens["sddb"]["sensitivity"],
            "n_sddb": sens["sddb"]["n"],
            "n_detected_sddb": sens["sddb"]["n_detected"],
            "median_lead_h_sddb": sens["sddb"]["median_lead_h"],
            "sens_vfdb": sens["vfdb"]["sensitivity"],
            "n_vfdb": sens["vfdb"]["n"],
            "n_detected_vfdb": sens["vfdb"]["n_detected"],
            "median_lead_h_vfdb": sens["vfdb"]["median_lead_h"],
            "sens_all_events": sens["all_events"]["sensitivity"],
            "n_all_events": sens["all_events"]["n"],
            "n_detected_all_events": sens["all_events"]["n_detected"],
            "median_lead_h_all_events": sens["all_events"]["median_lead_h"],
            # FAR / control specificity flavor
            "far_per_24h": far["far_per_24h"],
            "total_control_episodes": far["total_episodes"],
            "total_search_hours": far["total_search_hours"],
            "n_controls": far["n_controls"],
            "fraction_controls_alarmed": far["fraction_alarmed"],
            "median_far_per_record": far["median_far_per_record"],
            "mean_far_per_record": far["mean_far_per_record"],
            "min_far_per_record": far["min_far_per_record"],
            "max_far_per_record": far["max_far_per_record"],
            # Balance (exploratory)
            **bal,
        }
        rows.append(row)
    return rows


def build_cohort_summary_rows(tradeoff_rows: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """
    Long-format cohort table: one row per (detector, cohort).
    Cohorts: SDDB (events), VFDB (events), NSRDB (controls/FAR).
    """
    out: List[Dict[str, Any]] = []
    for r in tradeoff_rows:
        det = r["detector"]
        disp = r["detector_display"]
        out.append(
            {
                "detector": det,
                "detector_display": disp,
                "cohort": "SDDB",
                "role": "events_sensitivity",
                "n_records": r["n_sddb"],
                "n_detected_or_episodes": r["n_detected_sddb"],
                "sensitivity": r["sens_sddb"],
                "far_per_24h": None,
                "fraction_alarmed": None,
                "median_lead_h": r["median_lead_h_sddb"],
                "note": "Holter pre-VF event sensitivity (exploratory)",
            }
        )
        out.append(
            {
                "detector": det,
                "detector_display": disp,
                "cohort": "VFDB",
                "role": "events_sensitivity",
                "n_records": r["n_vfdb"],
                "n_detected_or_episodes": r["n_detected_vfdb"],
                "sensitivity": r["sens_vfdb"],
                "far_per_24h": None,
                "fraction_alarmed": None,
                "median_lead_h": r["median_lead_h_vfdb"],
                "note": "Telemetry pre-VF event sensitivity (exploratory)",
            }
        )
        out.append(
            {
                "detector": det,
                "detector_display": disp,
                "cohort": "all_events",
                "role": "events_sensitivity",
                "n_records": r["n_all_events"],
                "n_detected_or_episodes": r["n_detected_all_events"],
                "sensitivity": r["sens_all_events"],
                "far_per_24h": None,
                "fraction_alarmed": None,
                "median_lead_h": r["median_lead_h_all_events"],
                "note": "SDDB+VFDB pooled event sensitivity",
            }
        )
        out.append(
            {
                "detector": det,
                "detector_display": disp,
                "cohort": "NSRDB",
                "role": "controls_far",
                "n_records": r["n_controls"],
                "n_detected_or_episodes": int(r["total_control_episodes"]),
                "sensitivity": None,
                "far_per_24h": r["far_per_24h"],
                "fraction_alarmed": r["fraction_controls_alarmed"],
                "median_lead_h": None,
                "note": "Healthy Holter FAR; not device-matched; not classical specificity",
            }
        )
    return out


def rank_detectors(
    tradeoff_rows: List[Dict[str, Any]], key: str, higher_is_better: bool
) -> List[Dict[str, Any]]:
    """Observational ranking by a single numeric key."""
    ordered = sorted(
        tradeoff_rows,
        key=lambda r: (r[key] is None, -(r[key] or 0) if higher_is_better else (r[key] or 0)),
    )
    ranking = []
    for i, r in enumerate(ordered, start=1):
        ranking.append(
            {
                "rank": i,
                "detector": r["detector"],
                "detector_display": r["detector_display"],
                "value": r[key],
                "metric": key,
                "higher_is_better": higher_is_better,
            }
        )
    return ranking


def qualitative_notes(tradeoff_rows: List[Dict[str, Any]]) -> Dict[str, Any]:
    """Honest per-detector strengths/weaknesses and multi-criterion rankings."""
    by_id = {r["detector"]: r for r in tradeoff_rows}

    strengths_weaknesses = {
        "opc_L50": {
            "strengths": [
                "Lowest pooled FAR on NSRDB among the three primary arms under fixed params.",
                "Majority of NSRDB records have zero OPC episodes (median per-record FAR = 0).",
                "Collapse+persistence predicate is conceptually aligned with RECD ordinal structure.",
                "Cleaner basal behavior on SDDB in the bake-off vs saturated L=8.",
            ],
            "weaknesses": [
                "Lowest exploratory event sensitivity (SDDB 0.55, VFDB 0.36, all 0.42).",
                "Misses several Holters that SDD and abs-z detect.",
                "Still above informal S5 FAR target of ≤2/24h on this non-device-matched set (no S5 claim).",
                "Sensitive to L vs K alphabet interaction (L=8 is invalid for K=36).",
            ],
            "operating_point": "specificity-leaning / selective collapse detector",
        },
        "sdd": {
            "strengths": [
                "Highest exploratory event sensitivity among the three (SDDB 1.00, VFDB 0.95, all 0.97).",
                "Detects distributional reorganization (TV) rather than amplitude excursions.",
                "Matches or slightly exceeds abs-z hit rate on these public event sets.",
            ],
            "weaknesses": [
                "Highest pooled FAR on NSRDB (~46.3/24h) — near refractory-period ceiling on long Holters.",
                "Fraction of controls alarmed = 1.0 under these fixed params.",
                "Fixed early basal on ambulatory Holter dynamics is a hard specificity stress test.",
                "High sensitivity alone is not a usable operating point if control FAR saturates.",
            ],
            "operating_point": "sensitivity-leaning / distribution-shift detector",
        },
        "absz_tau_s": {
            "strengths": [
                "Frozen pilot baseline; parameters unchanged from production discovery freeze.",
                "High exploratory sensitivity (SDDB 1.00, VFDB 0.86, all 0.91).",
                "Intermediate FAR between OPC L=50 and SDD on NSRDB (~33.7/24h).",
                "Phase-2 FAR order of magnitude reproduced (sanity anchor).",
            ],
            "weaknesses": [
                "FAR remains high on healthy Holter relative to informal S5 (no S5 claim).",
                "All NSRDB controls alarmed at least once under z≥2 sustained ≥3.",
                "Continuous amplitude z-score is a different predicate than ordinal collapse/divergence.",
                "Three VFDB misses are in too_short pre-event strata.",
            ],
            "operating_point": "frozen continuous baseline; intermediate FAR vs high sens",
        },
    }

    rankings = {
        "best_specificity_lowest_far": rank_detectors(
            tradeoff_rows, "far_per_24h", higher_is_better=False
        ),
        "best_sensitivity_all_events": rank_detectors(
            tradeoff_rows, "sens_all_events", higher_is_better=True
        ),
        "best_sensitivity_sddb": rank_detectors(
            tradeoff_rows, "sens_sddb", higher_is_better=True
        ),
        "best_sensitivity_vfdb": rank_detectors(
            tradeoff_rows, "sens_vfdb", higher_is_better=True
        ),
        "best_sens_per_far_unit": rank_detectors(
            tradeoff_rows, "sens_per_far_unit", higher_is_better=True
        ),
        "best_geometric_balance": rank_detectors(
            tradeoff_rows, "geometric_balance_sens_x_spec_proxy", higher_is_better=True
        ),
        "best_youden_like_exploratory": rank_detectors(
            tradeoff_rows, "youden_like_exploratory", higher_is_better=True
        ),
        "fewest_control_episodes": rank_detectors(
            tradeoff_rows, "total_control_episodes", higher_is_better=False
        ),
    }

    # Objective-conditioned observational balance statements
    balance_by_objective = {
        "if_priority_is_lowest_false_alarm_burden": {
            "preferred_observational": "opc_L50",
            "rationale": (
                f"OPC L=50 pooled FAR {by_id['opc_L50']['far_per_24h']:.3f}/24h vs "
                f"abs-z {by_id['absz_tau_s']['far_per_24h']:.3f} and SDD {by_id['sdd']['far_per_24h']:.3f}; "
                "accepts substantially lower event hit rate."
            ),
        },
        "if_priority_is_maximum_event_hit_rate": {
            "preferred_observational": "sdd",
            "rationale": (
                f"SDD all-events sensitivity {by_id['sdd']['sens_all_events']:.3f} "
                f"(32/33) edges abs-z {by_id['absz_tau_s']['sens_all_events']:.3f} (30/33); "
                "FAR cost is highest of the three under these fixed params."
            ),
        },
        "if_priority_is_sens_per_unit_far": {
            "preferred_observational": "opc_L50",
            "rationale": (
                "Among fixed-param arms, OPC L=50 yields the highest sens_all / FAR ratio "
                "because FAR is ~9× lower than abs-z while retaining partial event detection."
            ),
        },
        "if_priority_is_frozen_baseline_continuity": {
            "preferred_observational": "absz_tau_s",
            "rationale": (
                "abs-z remains the only production-frozen continuous baseline; ordinal arms "
                "are methodological alternatives, not retunes of production thresholds."
            ),
        },
        "if_priority_is_balanced_hit_rate_with_intermediate_far": {
            "preferred_observational": "absz_tau_s",
            "rationale": (
                "abs-z sits between OPC (selective, low FAR) and SDD (near-ceiling FAR) with "
                "high event sensitivity; geometric/Youden-like proxies may still favor OPC "
                "because FAR dominates the denominator — interpret proxies cautiously."
            ),
        },
    }

    executive_summary = (
        "Exploratory trade-off on public PhysioNet cohorts (SDDB n=11, VFDB n=22 events; "
        "NSRDB n=18 controls, ~180 search-hours, refractory 0.5 h). Three detectors kept "
        "strictly separate with fixed parameters — no fusion, no production retune, no "
        "clinical/S5/FDA claim. "
        f"OPC L=50 is the most specific arm (FAR≈{by_id['opc_L50']['far_per_24h']:.2f}/24h, "
        f"{int(by_id['opc_L50']['total_control_episodes'])} episodes) but detects only "
        f"{by_id['opc_L50']['n_detected_all_events']}/{by_id['opc_L50']['n_all_events']} events "
        f"(sens≈{by_id['opc_L50']['sens_all_events']:.2f}). "
        f"SDD is the most sensitive (sens≈{by_id['sdd']['sens_all_events']:.2f}) with the "
        f"highest control FAR (≈{by_id['sdd']['far_per_24h']:.1f}/24h, "
        f"{int(by_id['sdd']['total_control_episodes'])} episodes). "
        f"Frozen abs-z is intermediate on FAR (≈{by_id['absz_tau_s']['far_per_24h']:.1f}/24h) "
        f"with high sensitivity (≈{by_id['absz_tau_s']['sens_all_events']:.2f}). "
        "Which arm looks 'balanced' depends on the objective: specificity → OPC L=50; "
        "hit rate → SDD; frozen pilot baseline → abs-z; sens-per-FAR efficiency → OPC L=50. "
        "NSRDB is not device-matched to VFDB/CU telemetry."
    )

    recommendations = [
        {
            "id": "A",
            "title": "Keep arms separate; do not fuse as primary claim",
            "detail": (
                "Continue evaluating OPC L=50, SDD, and abs-z independently. "
                "Any light fusion (e.g. sequential filter) would be a new experiment with "
                "its own multiplicity cost — not a free upgrade."
            ),
        },
        {
            "id": "B",
            "title": "If pursuing OPC, explore parameter/alphabet alignment — not abs-z retune",
            "detail": (
                "L vs K interaction is material (L=8 invalid for K=36). Modest grid on "
                "L, θ_D, θ_R on a held-out or bootstrap scheme could move the operating "
                "point without claiming clinical optimization. Do not retune production abs-z."
            ),
        },
        {
            "id": "C",
            "title": "If pursuing SDD, redesign basal/control policy before FAR claims",
            "detail": (
                "Current fixed early basal yields near-ceiling FAR on ambulatory Holters. "
                "Sliding basal, longer L_c, higher θ_TV, or θ_S>1 are natural specificity "
                "levers — each requires a pre-registered comparison, not post-hoc cherry-picking."
            ),
        },
        {
            "id": "D",
            "title": "Institutional device-matched controls (Tier A) for true FAR",
            "detail": (
                "Public NSRDB remains an interim upper-bound flavor only. Phase-2/S5-style "
                "specificity still needs partner Holter/telemetry matched to deployment setting."
            ),
        },
        {
            "id": "E",
            "title": "Optional light cascade (exploratory only, if ever tested)",
            "detail": (
                "A non-primary cascade such as 'SDD candidate → OPC confirm' or vice versa "
                "could be scored on the existing per-record tables without reprocessing RR, "
                "but only as a clearly labeled secondary analysis with no fusion claim in the "
                "primary manuscript arm."
            ),
        },
        {
            "id": "F",
            "title": "Manuscript framing",
            "detail": (
                "Present the three operating points as a trade-off surface, not a winner. "
                "Emphasize predicate differences (collapse vs distributional TV vs continuous z)."
            ),
        },
    ]

    return {
        "strengths_weaknesses": strengths_weaknesses,
        "rankings_observational": rankings,
        "balance_by_objective": balance_by_objective,
        "executive_summary": executive_summary,
        "recommendations": recommendations,
        "caveats": [
            "Exploratory public-data bake-off only; small event n (11+22).",
            "FAR ≠ classical 1−specificity; episode rate with 0.5 h refractory on Holter.",
            "NSRDB not device-matched to VFDB/CU; not institutional ICU telemetry.",
            "Timebases differ (strided τ_s vs symbol endpoints); FAR definition is shared.",
            "No threshold optimization / validation-set retune in this join analysis.",
            "No clinical utility, FDA readiness, or deployability is asserted; S5 (FAR ≤ 2/24h) is not asserted.",
            "OPC L=8 excluded as primary arm (parameter–alphabet saturation).",
            "No OPC∧SDD fusion scored.",
        ],
    }


def write_csv(path: Path, rows: List[Dict[str, Any]], fieldnames: Optional[List[str]] = None) -> None:
    if not rows:
        raise ValueError(f"No rows to write to {path}")
    fields = fieldnames or list(rows[0].keys())
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        for row in rows:
            out = {}
            for k in fields:
                v = row.get(k)
                if v is None:
                    out[k] = ""
                elif isinstance(v, float):
                    out[k] = f"{v:.10g}"
                else:
                    out[k] = v
            w.writerow(out)


def write_markdown_doc(
    path: Path,
    tradeoff_rows: List[Dict[str, Any]],
    cohort_rows: List[Dict[str, Any]],
    narrative: dict,
    meta: dict,
) -> None:
    lines: List[str] = []
    lines.append("# Ordinal vs abs-z: sensitivity–specificity trade-off (exploratory)")
    lines.append("")
    lines.append(f"*Generated: {meta['generated_at']}*")
    lines.append("")
    lines.append(
        "**Status:** Exploratory join of existing bake-off + NSRDB FAR results. "
        "**No clinical claim, no S5 claim, no superiority claim, no fusion.**"
    )
    lines.append("")
    lines.append("## Executive summary")
    lines.append("")
    lines.append(narrative["executive_summary"])
    lines.append("")
    lines.append("## Inputs (no reprocessing)")
    lines.append("")
    lines.append("| Artifact | Role |")
    lines.append("|----------|------|")
    for p in meta["inputs"]:
        lines.append(f"| `{p}` | existing |")
    lines.append("")
    lines.append("## Primary arms (strictly separate)")
    lines.append("")
    lines.append("| Detector | Sensitivity source key | FAR source key |")
    lines.append("|----------|------------------------|----------------|")
    lines.append("| OPC L=50 | `opc_L50_companion` | `opc_L50` |")
    lines.append("| SDD (TV) | `sdd` | `sdd` |")
    lines.append("| abs-z τ_s frozen | `absz_tau_s` | `absz_tau_s` |")
    lines.append("")
    lines.append("OPC L=8 is **not** a primary arm (K/L saturation with joint K=36).")
    lines.append("")
    lines.append("## Trade-off table (sensitivity vs FAR)")
    lines.append("")
    lines.append(
        "| Detector | Sens SDDB | Sens VFDB | Sens all | FAR /24h | "
        "Control episodes | Frac. controls alarmed | Sens/FAR |"
    )
    lines.append(
        "|----------|-----------|-----------|----------|----------|"
        "-------------------|------------------------|----------|"
    )
    for r in tradeoff_rows:
        lines.append(
            f"| {r['detector_display']} | "
            f"{r['sens_sddb']:.3f} ({r['n_detected_sddb']}/{r['n_sddb']}) | "
            f"{r['sens_vfdb']:.3f} ({r['n_detected_vfdb']}/{r['n_vfdb']}) | "
            f"{r['sens_all_events']:.3f} ({r['n_detected_all_events']}/{r['n_all_events']}) | "
            f"**{r['far_per_24h']:.3f}** | "
            f"{int(r['total_control_episodes'])} | "
            f"{r['fraction_controls_alarmed']:.3f} | "
            f"{r['sens_per_far_unit']:.4f} |"
        )
    lines.append("")
    lines.append("## Cohort summary")
    lines.append("")
    lines.append(
        "| Detector | Cohort | Role | n | Detected / episodes | Sensitivity | FAR /24h |"
    )
    lines.append(
        "|----------|--------|------|---|---------------------|-------------|----------|"
    )
    for r in cohort_rows:
        sens = f"{r['sensitivity']:.3f}" if r["sensitivity"] is not None else "—"
        far = f"{r['far_per_24h']:.3f}" if r["far_per_24h"] is not None else "—"
        lines.append(
            f"| {r['detector_display']} | {r['cohort']} | {r['role']} | "
            f"{r['n_records']} | {r['n_detected_or_episodes']} | {sens} | {far} |"
        )
    lines.append("")
    lines.append("## Observational rankings (multi-criterion)")
    lines.append("")
    for name, ranking in narrative["rankings_observational"].items():
        direction = "higher better" if ranking[0]["higher_is_better"] else "lower better"
        lines.append(f"### {name} ({direction})")
        lines.append("")
        for item in ranking:
            val = item["value"]
            if isinstance(val, float):
                val_s = f"{val:.6g}"
            else:
                val_s = str(val)
            lines.append(f"{item['rank']}. **{item['detector_display']}**: {val_s}")
        lines.append("")
    lines.append("## Strengths and weaknesses (qualitative)")
    lines.append("")
    for det, sw in narrative["strengths_weaknesses"].items():
        lines.append(f"### {DISPLAY_NAMES[det]}")
        lines.append("")
        lines.append(f"*Operating point:* {sw['operating_point']}")
        lines.append("")
        lines.append("**Strengths**")
        for s in sw["strengths"]:
            lines.append(f"- {s}")
        lines.append("")
        lines.append("**Weaknesses**")
        for w in sw["weaknesses"]:
            lines.append(f"- {w}")
        lines.append("")
    lines.append("## Balance by objective")
    lines.append("")
    for obj, block in narrative["balance_by_objective"].items():
        lines.append(f"- **{obj}** → observational preference: "
                     f"`{block['preferred_observational']}` — {block['rationale']}")
    lines.append("")
    lines.append("## Recommendations (non-binding)")
    lines.append("")
    for rec in narrative["recommendations"]:
        lines.append(f"### {rec['id']}. {rec['title']}")
        lines.append("")
        lines.append(rec["detail"])
        lines.append("")
    lines.append("## Caveats / non-claims")
    lines.append("")
    for c in narrative["caveats"]:
        lines.append(f"- {c}")
    lines.append("")
    lines.append("## Artifacts")
    lines.append("")
    for p in meta["outputs"]:
        lines.append(f"- `{p}`")
    lines.append(f"- Runner: `code/run_ordinal_tradeoff_analysis.py`")
    lines.append("")

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def run_analysis(
    results_dir: Path = RES,
    write_doc: bool = True,
) -> Dict[str, Any]:
    """Load existing artifacts, join trade-off, write outputs."""
    exploratory_path = results_dir / "ordinal_exploratory_summary.json"
    far_path = results_dir / "ordinal_nsrdb_far_summary.json"
    per_record_path = results_dir / "ordinal_nsrdb_far_per_record.csv"
    comparison_path = results_dir / "ordinal_opc_sdd_absz_comparison.csv"

    for p in (exploratory_path, far_path, per_record_path, comparison_path):
        if not p.is_file():
            raise FileNotFoundError(f"Required input missing: {p}")

    exploratory = load_json(exploratory_path)
    far_summary = load_json(far_path)
    per_record = load_csv_rows(per_record_path)

    tradeoff_rows = build_tradeoff_rows(exploratory, far_summary)
    cohort_rows = build_cohort_summary_rows(tradeoff_rows)
    narrative = qualitative_notes(tradeoff_rows)

    # Episode sum cross-check from per-record CSV
    episode_columns = {
        "opc_L50": "opc_n_episodes",
        "sdd": "sdd_n_episodes",
        "absz_tau_s": "absz_n_episodes",
    }
    episode_crosscheck = {}
    for arm_id, col in episode_columns.items():
        summed = sum_episodes_from_per_record(per_record, col)
        pooled = next(r for r in tradeoff_rows if r["detector"] == arm_id)[
            "total_control_episodes"
        ]
        episode_crosscheck[arm_id] = {
            "sum_per_record": summed,
            "pooled_total_episodes": pooled,
            "match": abs(summed - pooled) < 1e-6,
        }

    # Sanity: fusion flags false
    fusion_flags = {
        "exploratory_fusion": exploratory.get("params", {}).get("fusion", None),
        "far_fusion": far_summary.get("fusion", None),
        "tradeoff_any_fusion": any(r.get("fusion") for r in tradeoff_rows),
    }

    # Forbidden claim scan on narrative text
    blob = json.dumps(narrative).lower()
    forbidden_hits = [p for p in FORBIDDEN_CLAIM_PHRASES if p in blob]

    out_tradeoff_csv = results_dir / "ordinal_sensitivity_specificity_tradeoff.csv"
    out_cohort_csv = results_dir / "ordinal_tradeoff_by_cohort.csv"
    out_summary_json = results_dir / "ordinal_tradeoff_summary.json"
    out_doc = DOCS / "ORDINAL_SENSITIVITY_SPECIFICITY_TRADEOFF.md"

    write_csv(out_tradeoff_csv, tradeoff_rows)
    write_csv(out_cohort_csv, cohort_rows)

    summary: Dict[str, Any] = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "phase": "ordinal_sensitivity_specificity_tradeoff",
        "exploratory_only": True,
        "clinical_claim": False,
        "superiority_claim": False,
        "s5_claim": False,
        "fusion": False,
        "primary_arms": [a[0] for a in PRIMARY_ARMS],
        "excluded_from_primary": {
            "opc_L8": "K/L saturation with joint K=36; not a usable specificity probe",
        },
        "inputs": {
            "exploratory_summary": safe_relpath(exploratory_path),
            "far_summary": safe_relpath(far_path),
            "far_per_record": safe_relpath(per_record_path),
            "comparison_csv": safe_relpath(comparison_path),
            "bakeoff_doc": "docs/ORDINAL_EXPLORATORY_BAKEOFF.md",
            "far_doc": "docs/ORDINAL_NSRDB_FAR_COMPARISON.md",
        },
        "methodology": {
            "join_only": True,
            "reprocessed_wfdb": False,
            "sensitivity_source": "ordinal_exploratory_summary.json (SDDB+VFDB bake-off)",
            "far_source": "ordinal_nsrdb_far_summary.json (NSRDB n=18, refractory 0.5h)",
            "far_formula": "total_episodes / total_search_hours * 24",
            "device_mismatch": True,
            "device_mismatch_note": (
                "NSRDB is rhythm-healthy Holter ECG — NOT device-matched to VFDB/CU telemetry."
            ),
            "balance_metrics_note": (
                "sens_per_far_unit, specificity_proxy, geometric_balance, youden_like are "
                "exploratory ranking aids only — not clinical decision metrics."
            ),
        },
        "tradeoff_table": tradeoff_rows,
        "cohort_summary": cohort_rows,
        "episode_crosscheck_per_record_vs_pooled": episode_crosscheck,
        "fusion_flags": fusion_flags,
        "forbidden_claim_scan": {
            "phrases_checked": list(FORBIDDEN_CLAIM_PHRASES),
            "hits": forbidden_hits,
        },
        "narrative": narrative,
        "outputs": {
            "tradeoff_csv": safe_relpath(out_tradeoff_csv),
            "cohort_csv": safe_relpath(out_cohort_csv),
            "summary_json": safe_relpath(out_summary_json),
            "doc": safe_relpath(out_doc) if write_doc else None,
        },
    }

    with open(out_summary_json, "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)
        f.write("\n")

    if write_doc:
        write_markdown_doc(
            out_doc,
            tradeoff_rows,
            cohort_rows,
            narrative,
            meta={
                "generated_at": summary["generated_at"],
                "inputs": list(summary["inputs"].values()),
                "outputs": [
                    summary["outputs"]["tradeoff_csv"],
                    summary["outputs"]["cohort_csv"],
                    summary["outputs"]["summary_json"],
                    summary["outputs"]["doc"],
                ],
            },
        )

    return summary


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description="Join existing ordinal bake-off + FAR results into a trade-off analysis"
    )
    parser.add_argument(
        "--results-dir",
        type=Path,
        default=RES,
        help="Directory with ordinal_* result artifacts",
    )
    parser.add_argument(
        "--no-doc",
        action="store_true",
        help="Skip writing docs/ORDINAL_SENSITIVITY_SPECIFICITY_TRADEOFF.md",
    )
    args = parser.parse_args(argv)

    summary = run_analysis(results_dir=args.results_dir, write_doc=not args.no_doc)

    print("=== Ordinal sensitivity–specificity trade-off (exploratory) ===")
    print(f"Generated: {summary['generated_at']}")
    print(f"Fusion: {summary['fusion']} | Clinical claim: {summary['clinical_claim']}")
    print()
    print(
        f"{'Detector':<22} {'SensSDDB':>8} {'SensVFDB':>8} {'SensAll':>8} "
        f"{'FAR/24h':>10} {'Episodes':>9} {'Sens/FAR':>10}"
    )
    for r in summary["tradeoff_table"]:
        print(
            f"{r['detector_display']:<22} "
            f"{r['sens_sddb']:8.3f} {r['sens_vfdb']:8.3f} {r['sens_all_events']:8.3f} "
            f"{r['far_per_24h']:10.3f} {int(r['total_control_episodes']):9d} "
            f"{r['sens_per_far_unit']:10.4f}"
        )
    print()
    print("Episode cross-check (per-record sum vs pooled):")
    for arm, chk in summary["episode_crosscheck_per_record_vs_pooled"].items():
        status = "OK" if chk["match"] else "MISMATCH"
        print(f"  {arm}: sum={chk['sum_per_record']} pooled={chk['pooled_total_episodes']} [{status}]")
    print()
    print("Outputs:")
    for k, v in summary["outputs"].items():
        if v:
            print(f"  {k}: {v}")
    print()
    print("Executive summary:")
    print(summary["narrative"]["executive_summary"])
    return 0


if __name__ == "__main__":
    sys.exit(main())
