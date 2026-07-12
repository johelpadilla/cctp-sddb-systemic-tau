#!/usr/bin/env python3
"""
Extend cctp_batch_summary.csv with stratification fields and write stratified tables.

Inputs:
  results/cctp_batch_summary.csv
  data/records_inventory.csv (or sddb_full_inventory.csv)
  data/rr_{rec}_clean.npz  (authoritative event_hr via get_event_and_windows)

Outputs:
  results/cctp_cohort_stratified.csv   (one row per analytic record + flags)
  results/cctp_stratified_summary.csv  (group means)
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cctp_metrics_core import (
    event_type_label,
    resolve_event_timing_from_npz,
    stratified_summary,
    substrate_label,
)

BASE = Path(__file__).resolve().parent.parent
RES = BASE / "results"
DATA = BASE / "data"


def load_csv(path: Path) -> list:
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def write_csv(path: Path, rows: list, fieldnames=None):
    if not rows:
        path.write_text("")
        return
    fieldnames = fieldnames or list(rows[0].keys())
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def _finite_float(v):
    if v in (None, ""):
        return None
    try:
        x = float(v)
        if x != x:  # NaN
            return None
        return x
    except (TypeError, ValueError):
        return None


def resolve_record_event_timing(rec: str, batch_row: dict) -> dict:
    """
    Resolve event_hr / duration / event_type without using duration as a fake event.

    Priority:
      1. Clean RR npz + get_event_and_windows (authoritative)
      2. Finite event_hr already on the batch row + duration_h
      3. unknown
    """
    npz_path = DATA / f"rr_{rec}_clean.npz"
    if npz_path.exists():
        timing = resolve_event_timing_from_npz(rec, str(npz_path))
        return {
            "event_hr": timing["event_hr"],
            "duration_h": timing["duration_h"],
            "event_type": timing["event_type"],
            "event_hr_source": "npz",
        }

    event_hr = _finite_float(batch_row.get("event_hr"))
    dur = _finite_float(batch_row.get("duration_h"))
    if event_hr is not None and dur is not None:
        return {
            "event_hr": event_hr,
            "duration_h": dur,
            "event_type": event_type_label(rec, event_hr, dur),
            "event_hr_source": "batch_csv",
        }
    return {
        "event_hr": event_hr,
        "duration_h": dur,
        "event_type": "unknown",
        "event_hr_source": "missing",
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--batch", default=str(RES / "cctp_batch_summary.csv"))
    ap.add_argument("--inventory", default=str(DATA / "sddb_full_inventory.csv"))
    ap.add_argument("--fallback-inventory", default=str(DATA / "records_inventory.csv"))
    ap.add_argument("--out-cohort", default=str(RES / "cctp_cohort_stratified.csv"))
    ap.add_argument("--out-strata", default=str(RES / "cctp_stratified_summary.csv"))
    args = ap.parse_args()

    batch = load_csv(Path(args.batch))
    inv_path = Path(args.inventory)
    if not inv_path.exists():
        inv_path = Path(args.fallback_inventory)
    inv_rows = load_csv(inv_path)
    inv_by = {}
    for r in inv_rows:
        rid = str(r.get("record_id") or r.get("record") or "").strip()
        inv_by[rid] = r

    enriched = []
    for m in batch:
        rec = str(m.get("record", "")).strip()
        inv = inv_by.get(rec, {})
        rhythm = inv.get("rhythm", "")
        pacing = inv.get("pacing") or m.get("known_pacing_type") or "none"
        if str(m.get("pacing_detected", "")).lower() in ("true", "1", "yes"):
            if pacing in ("none", "", None):
                pacing = "intermittent"
        substrate = substrate_label(rhythm, str(pacing))

        timing = resolve_record_event_timing(rec, m)
        event_hr = timing["event_hr"]
        dur = timing["duration_h"] if timing["duration_h"] is not None else _finite_float(m.get("duration_h"))
        etype = timing["event_type"]

        def fget(k):
            v = m.get(k)
            if v in (None, ""):
                return ""
            try:
                return float(v)
            except ValueError:
                return v

        row = {
            **{k: m.get(k) for k in m},
            "record": rec,
            "duration_h": dur if dur is not None else m.get("duration_h"),
            "event_hr": event_hr if event_hr is not None else "",
            "event_hr_source": timing["event_hr_source"],
            "rhythm": rhythm,
            "pacing": pacing,
            "substrate": substrate,
            "event_type": etype,
            "pre_event_hours": inv.get("pre_event_hours", ""),
            "inventory_include": inv.get("inventory_include") or inv.get("include?", ""),
            "in_manuscript_n10": inv.get("in_manuscript_n10", "Yes"),
            "quality_flag_interp": (
                "high" if m.get("interp_frac") not in (None, "") and float(m["interp_frac"]) > 0.05 else "ok"
            ),
            "quality_flag_pacing": "paced" if substrate == "paced" else "none",
            "delta_tau": fget("delta_tau"),
            "delta_excess3": fget("delta_excess3"),
            "delta_var": fget("delta_var"),
            "delta_ar1": fget("delta_ar1"),
            "interp_frac": fget("interp_frac"),
            "cv_rr": fget("cv_rr"),
        }
        enriched.append(row)

    pref = [
        "record",
        "duration_h",
        "event_hr",
        "event_hr_source",
        "substrate",
        "event_type",
        "rhythm",
        "pacing",
        "delta_tau",
        "delta_excess3",
        "delta_var",
        "delta_ar1",
        "p_tau",
        "p_excess3",
        "p_tau_surrogate",
        "interp_frac",
        "pacing_detected",
        "known_pacing_type",
        "cv_rr",
        "n_beats",
        "quality_flag_interp",
        "quality_flag_pacing",
        "pre_event_hours",
        "in_manuscript_n10",
        "has_weighted",
    ]
    keys = list(dict.fromkeys(pref + [k for r in enriched for k in r.keys()]))
    write_csv(Path(args.out_cohort), enriched, fieldnames=keys)

    strata = []
    for gk in ("substrate", "event_type", "quality_flag_pacing"):
        strata.extend(stratified_summary(enriched, gk))
    write_csv(Path(args.out_strata), strata)

    n_int = sum(1 for r in enriched if r["event_type"] == "intermediate")
    n_term = sum(1 for r in enriched if r["event_type"] == "terminal")
    print(f"Wrote {args.out_cohort} ({len(enriched)} analytic rows)")
    print(f"Wrote {args.out_strata} ({len(strata)} stratum rows)")
    print(f"event_type counts: intermediate={n_int} terminal={n_term}")
    for r in enriched:
        print(
            f"  {r['record']}: substrate={r['substrate']} event_type={r['event_type']} "
            f"event_hr={r['event_hr']} src={r['event_hr_source']} "
            f"Δτ={r['delta_tau']} Δex3={r['delta_excess3']}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
