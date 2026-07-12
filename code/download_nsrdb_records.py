#!/usr/bin/env python3
"""
download_nsrdb_records.py

Surgical downloader for PhysioNet NSRDB .hea + .atr only
(beat annotations suffice for RR extraction; no .dat required).

Usage:
  python3 code/download_nsrdb_records.py --remaining
  python3 code/download_nsrdb_records.py --records 16773,16786
  python3 code/download_nsrdb_records.py --all
"""
from __future__ import annotations

import argparse
import ssl
import urllib.error
import urllib.request
from pathlib import Path

try:
    _ssl_ctx = ssl._create_unverified_context()
    urllib.request.install_opener(
        urllib.request.build_opener(urllib.request.HTTPSHandler(context=_ssl_ctx))
    )
except Exception:
    pass

BASE = Path(__file__).resolve().parent.parent
DATA_NSR = BASE / "data" / "nsrdb"
DATA_NSR.mkdir(parents=True, exist_ok=True)

NSRDB_BASE = "https://physionet.org/files/nsrdb/1.0.0"

# Full NSRDB record list (PhysioNet nsrdb 1.0.0)
FULL_NSRDB = [
    "16265",
    "16272",
    "16273",
    "16420",
    "16483",
    "16539",
    "16773",
    "16786",
    "16795",
    "17052",
    "17453",
    "18177",
    "18184",
    "19088",
    "19090",
    "19093",
    "19140",
    "19830",
]


def download_file(rec: str, ext: str, force: bool = False) -> bool:
    local = DATA_NSR / f"{rec}.{ext}"
    if local.exists() and not force:
        print(f"  [skip] {rec}.{ext} already present")
        return True
    url = f"{NSRDB_BASE}/{rec}.{ext}"
    try:
        print(f"  downloading {rec}.{ext} ...")
        urllib.request.urlretrieve(url, local)
        print(f"  saved {local}")
        return True
    except urllib.error.HTTPError as e:
        print(f"  ERROR {rec}.{ext}: HTTP {e.code}")
        if local.exists() and local.stat().st_size == 0:
            local.unlink()
        return False
    except Exception as e:
        print(f"  ERROR {rec}.{ext}: {e}")
        if local.exists() and local.stat().st_size == 0:
            local.unlink()
        return False


def present_records() -> list[str]:
    return sorted({p.stem for p in DATA_NSR.glob("*.hea") if (DATA_NSR / f"{p.stem}.atr").exists()})


def update_records_lists(present: list[str]) -> None:
    (DATA_NSR / "RECORDS").write_text("\n".join(present) + ("\n" if present else ""))
    (DATA_NSR / "RECORDS.local").write_text("\n".join(present) + ("\n" if present else ""))
    (DATA_NSR / "RECORDS.full").write_text("\n".join(FULL_NSRDB) + "\n")


def main() -> int:
    ap = argparse.ArgumentParser(description="Download PhysioNet NSRDB hea+atr for controls")
    ap.add_argument("--records", default="", help="Comma-separated record ids")
    ap.add_argument("--remaining", action="store_true", help="Download only missing FULL_NSRDB")
    ap.add_argument("--all", action="store_true", help="Ensure all FULL_NSRDB hea+atr present")
    ap.add_argument("--force", action="store_true", help="Redownload even if present")
    args = ap.parse_args()

    if args.records:
        recs = [r.strip() for r in args.records.split(",") if r.strip()]
    elif args.remaining or args.all:
        have = set(present_records())
        if args.remaining:
            recs = [r for r in FULL_NSRDB if r not in have]
        else:
            recs = list(FULL_NSRDB)
    else:
        have = set(present_records())
        recs = [r for r in FULL_NSRDB if r not in have]

    print(f"Target NSRDB records: {recs} ({len(recs)})")
    ok = 0
    fail = 0
    for r in recs:
        print(f"\nRecord {r}:")
        hea_ok = download_file(r, "hea", args.force)
        atr_ok = download_file(r, "atr", args.force)
        if hea_ok and atr_ok:
            ok += 1
        else:
            fail += 1

    present = present_records()
    update_records_lists(present)
    print(f"\nDone. newly_ok={ok} fail={fail}; local present hea+atr={len(present)}/{len(FULL_NSRDB)}")
    print(f"Present: {present}")
    return 0 if fail == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
