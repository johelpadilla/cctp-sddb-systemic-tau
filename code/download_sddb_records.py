#!/usr/bin/env python3
"""
download_sddb_records.py

Surgical downloader for SDDB .hea + .atr only (no .dat needed for RR extraction).
Reuses wfdb or falls back to stdlib urllib for reproducibility.

Usage:
  python3 code/download_sddb_records.py --records 30,31,35
  python3 code/download_sddb_records.py --list selected_records.txt
  python3 code/download_sddb_records.py --all-candidates   # from inventory with include?=Yes

After download, you can run extraction for each.
"""

import argparse
import os
import sys
import urllib.request
import urllib.error
import ssl
from pathlib import Path

# Surgical SSL patch for macOS envs with CERTIFICATE_VERIFY_FAILED on PhysioNet
try:
    _ssl_ctx = ssl._create_unverified_context()
    urllib.request.install_opener(urllib.request.build_opener(urllib.request.HTTPSHandler(context=_ssl_ctx)))
except Exception:
    pass

BASE = Path(__file__).resolve().parent.parent
DATA_SDD = BASE / "data" / "sddb"
DATA_SDD.mkdir(parents=True, exist_ok=True)

SDDB_BASE = "https://physionet.org/files/sddb/1.0.0"

def download_file(rec: str, ext: str, force: bool = False) -> bool:
    local = DATA_SDD / f"{rec}.{ext}"
    if local.exists() and not force:
        print(f"  [skip] {rec}.{ext} already present")
        return True
    url = f"{SDDB_BASE}/{rec}.{ext}"
    try:
        print(f"  downloading {rec}.{ext} ...")
        urllib.request.urlretrieve(url, local)
        print(f"  saved {local}")
        return True
    except urllib.error.HTTPError as e:
        print(f"  ERROR {rec}.{ext}: HTTP {e.code}")
        return False
    except Exception as e:
        print(f"  ERROR {rec}.{ext}: {e}")
        return False

def load_records_from_list(path: Path):
    recs = []
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        recs.append(line)
    return recs

def load_inventory_included():
    import csv
    inv = BASE / "data" / "records_inventory.csv"
    if not inv.exists():
        print("No inventory found. Run inventory generation first.")
        return []
    recs = []
    with inv.open() as f:
        for row in csv.DictReader(f):
            if row.get("include?") == "Yes":
                recs.append(row["record_id"])
    return recs

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--records", default="", help="Comma separated record ids e.g. 30,31,35")
    ap.add_argument("--list", default="", help="Path to selected_records.txt")
    ap.add_argument("--all-candidates", action="store_true", help="Use include?=Yes from inventory")
    ap.add_argument("--force", action="store_true", help="Redownload even if present")
    args = ap.parse_args()

    recs = []
    if args.records:
        recs = [r.strip() for r in args.records.split(",") if r.strip()]
    elif args.list:
        recs = load_records_from_list(Path(args.list))
    elif args.all_candidates:
        recs = load_inventory_included()
    else:
        # default to current selected
        sel = BASE / "selected_records.txt"
        if sel.exists():
            recs = load_records_from_list(sel)
        else:
            recs = ["30", "35"]

    print(f"Target records: {recs}")
    ok = 0
    for r in recs:
        print(f"\nRecord {r}:")
        hea_ok = download_file(r, "hea", args.force)
        atr_ok = download_file(r, "atr", args.force)
        if hea_ok and atr_ok:
            ok += 1
    print(f"\nDone. {ok}/{len(recs)} records have .hea+.atr present or downloaded.")

if __name__ == "__main__":
    main()
