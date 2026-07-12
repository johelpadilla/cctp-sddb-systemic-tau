#!/usr/bin/env python3
"""
Build complete SDDB (23-record) inventory with inclusion/exclusion, annotation
availability, analytic-cohort membership, and process status.

Outputs:
  data/sddb_full_inventory.csv
  results/sddb_inventory_summary.json
"""
from __future__ import annotations

import argparse
import csv
import json
from datetime import datetime, timezone
from pathlib import Path

BASE = Path(__file__).resolve().parent.parent
DATA = BASE / "data"
RES = BASE / "results"
SDDB = DATA / "sddb"
RES.mkdir(parents=True, exist_ok=True)

# Authoritative N=10 manuscript analytic cohort
ANALYTIC_N10 = {"30", "31", "32", "35", "36", "38", "45", "47", "50", "51"}


def ann_status(rec: str) -> dict:
    hea = (SDDB / f"{rec}.hea").exists()
    atr = (SDDB / f"{rec}.atr").exists()
    ari = (SDDB / f"{rec}.ari").exists()
    dat = (SDDB / f"{rec}.dat").exists()
    dat_size = (SDDB / f"{rec}.dat").stat().st_size if dat else 0
    clean = (DATA / f"rr_{rec}_clean.npz").exists()
    return {
        "has_hea": hea,
        "has_atr": atr,
        "has_ari": ari,
        "has_annotation": atr or ari,
        "has_dat": dat and dat_size > 0,
        "dat_bytes": dat_size,
        "has_clean_rr": clean,
        "can_extract_rr": (atr or ari) and hea,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--inventory", default=str(DATA / "records_inventory.csv"))
    ap.add_argument("--out-csv", default=str(DATA / "sddb_full_inventory.csv"))
    ap.add_argument("--out-json", default=str(RES / "sddb_inventory_summary.json"))
    args = ap.parse_args()

    inv_path = Path(args.inventory)
    if not inv_path.exists():
        raise FileNotFoundError(inv_path)

    with open(inv_path, newline="") as f:
        base_rows = list(csv.DictReader(f))

    out_rows = []
    for r in base_rows:
        rec = str(r["record_id"]).strip()
        st = ann_status(rec)
        include = str(r.get("include?", "")).strip().lower() in ("yes", "y", "true", "1")
        in_n10 = rec in ANALYTIC_N10
        expansion_candidate = include and st["can_extract_rr"] and not st["has_clean_rr"]
        if in_n10 and st["has_clean_rr"]:
            process_status = "processed_analytic"
            exclusion_reason = ""
        elif st["has_clean_rr"] and not in_n10:
            process_status = "rr_extracted_expansion"
            exclusion_reason = "clean RR available; not in manuscript N=10 analytic set"
        elif not include:
            process_status = "excluded"
            exclusion_reason = r.get("reason", "")
        elif not st["can_extract_rr"]:
            process_status = "eligible_no_local_annotations"
            exclusion_reason = "include=Yes but no local .atr/.ari+.hea for RR extraction"
        else:
            process_status = "eligible_unprocessed"
            exclusion_reason = "annotations present; RR not yet extracted (expansion candidate)"

        row = {
            "record_id": rec,
            "duration_hours": r.get("duration_hours", ""),
            "has_vf_marker": r.get("has_vf_marker", ""),
            "pre_event_hours": r.get("pre_event_hours", ""),
            "has_atr_clinical": r.get("has_atr", ""),
            "pacing": r.get("pacing", "none"),
            "rhythm": r.get("rhythm", ""),
            "vf_onset": r.get("vf_onset", ""),
            "quality_notes": r.get("quality_notes", ""),
            "inventory_include": "Yes" if include else "No",
            "inventory_reason": r.get("reason", ""),
            "in_manuscript_n10": "Yes" if in_n10 else "No",
            "process_status": process_status,
            "exclusion_or_hold_reason": exclusion_reason or r.get("reason", ""),
            "expansion_candidate": "Yes" if expansion_candidate else "No",
            **{k: (int(v) if isinstance(v, bool) else v) for k, v in st.items()},
            "has_hea": int(st["has_hea"]),
            "has_atr": int(st["has_atr"]),
            "has_ari": int(st["has_ari"]),
            "has_annotation": int(st["has_annotation"]),
            "has_dat": int(st["has_dat"]),
            "has_clean_rr": int(st["has_clean_rr"]),
            "can_extract_rr": int(st["can_extract_rr"]),
        }
        # fix bools already cast
        for bk in ("has_hea", "has_atr", "has_ari", "has_annotation", "has_dat", "has_clean_rr", "can_extract_rr"):
            row[bk] = int(st[bk.replace("has_annotation", "has_annotation") if bk != "has_annotation" else "has_annotation"])
        row["has_hea"] = int(st["has_hea"])
        row["has_atr"] = int(st["has_atr"])
        row["has_ari"] = int(st["has_ari"])
        row["has_annotation"] = int(st["has_annotation"])
        row["has_dat"] = int(st["has_dat"])
        row["has_clean_rr"] = int(st["has_clean_rr"])
        row["can_extract_rr"] = int(st["can_extract_rr"])
        out_rows.append(row)

    # ensure all 23 SDDB ids 30-52 represented
    known = {r["record_id"] for r in out_rows}
    for i in range(30, 53):
        rid = str(i)
        if rid not in known:
            st = ann_status(rid)
            out_rows.append(
                {
                    "record_id": rid,
                    "duration_hours": "",
                    "has_vf_marker": "",
                    "pre_event_hours": "",
                    "has_atr_clinical": "",
                    "pacing": "unknown",
                    "rhythm": "unknown",
                    "vf_onset": "",
                    "quality_notes": "not in base inventory",
                    "inventory_include": "No",
                    "inventory_reason": "missing from records_inventory.csv",
                    "in_manuscript_n10": "Yes" if rid in ANALYTIC_N10 else "No",
                    "process_status": "unknown_not_in_base_inventory",
                    "exclusion_or_hold_reason": "not listed in records_inventory.csv",
                    "expansion_candidate": "No",
                    "has_hea": int(st["has_hea"]),
                    "has_atr": int(st["has_atr"]),
                    "has_ari": int(st["has_ari"]),
                    "has_annotation": int(st["has_annotation"]),
                    "has_dat": int(st["has_dat"]),
                    "dat_bytes": st["dat_bytes"],
                    "has_clean_rr": int(st["has_clean_rr"]),
                    "can_extract_rr": int(st["can_extract_rr"]),
                }
            )

    out_rows.sort(key=lambda x: int(x["record_id"]))
    fieldnames = list(out_rows[0].keys())
    out_csv = Path(args.out_csv)
    with open(out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(out_rows)

    n_total = len(out_rows)
    n_n10 = sum(1 for r in out_rows if r["in_manuscript_n10"] == "Yes")
    n_proc = sum(1 for r in out_rows if r["process_status"] == "processed_analytic")
    n_excl = sum(1 for r in out_rows if r["process_status"] == "excluded")
    n_exp = sum(1 for r in out_rows if r["expansion_candidate"] == "Yes")
    n_noann = sum(1 for r in out_rows if r["process_status"] == "eligible_no_local_annotations")

    summary = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "n_sddb_records": n_total,
        "n_manuscript_analytic": n_n10,
        "n_processed_with_clean_rr": n_proc,
        "n_excluded": n_excl,
        "n_expansion_candidates_with_local_ann": n_exp,
        "n_eligible_missing_local_annotations": n_noann,
        "analytic_record_ids": sorted(ANALYTIC_N10, key=int),
        "csv": str(out_csv.relative_to(BASE)),
        "notes": (
            "Full SDDB has 23 Holter records (30–52). Analytic cohort N=10 has cleaned RR npz. "
            "Expansion requires local annotations; many records only have inventory metadata without .atr/.ari on disk."
        ),
    }
    out_json = Path(args.out_json)
    with open(out_json, "w") as f:
        json.dump(summary, f, indent=2)

    print(f"Wrote {out_csv} ({n_total} rows)")
    print(f"Wrote {out_json}")
    print(
        f"Summary: total={n_total} analytic_n10={n_n10} processed_rr={n_proc} "
        f"excluded={n_excl} expansion_local={n_exp} eligible_no_ann={n_noann}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
