#!/usr/bin/env python3
"""
extract_rr.py

Reusable surgical RR extractor + cleaner for SDDB records.
Input: data/sddb/{rec}.hea + .atr (no .dat needed)
Output: data/rr_{rec}_clean.npz with:
  rr_ms, t_sec, fs_ann, vfon_sec, vfon_beat_idx, n_beats, total_hours,
  plus quality: n_invalid, interp_frac, pacing_detected, known_pacing_type

Usage:
  python3 code/extract_rr.py --record 31
  python3 code/extract_rr.py --record all   # from selected_records.txt or inventory include=Yes

Cleans:
- RR intervals outside [250, 2000] ms marked invalid
- Linear interpolation across invalid gaps (small gaps)
- vfon parsed from .hea "#vfon: HH:MM:SS"
- Auto-detects pacing (comments + known list + RR regularity heuristic)
- Extracts quality metadata for batch reports
"""

import argparse
import os
import re
import numpy as np
from pathlib import Path
import wfdb

BASE = Path(__file__).resolve().parent.parent
DATA_DIR = BASE / "data"
SDDB_DIR = DATA_DIR / "sddb"
OUT_DIR = DATA_DIR
OUT_DIR.mkdir(parents=True, exist_ok=True)

def parse_vfon_from_hea(hea_path: Path) -> float:
    """Parse '#vfon: HH:MM:SS' or similar from .hea comments -> seconds from start."""
    if not hea_path.exists():
        return None
    txt = hea_path.read_text()
    # Primary: #vfon: HH:MM:SS  (handles spaces, as seen in SDDB headers)
    m = re.search(r"#\s*vfon\s*:\s*(\d+):(\d+):(\d+)", txt, re.IGNORECASE)
    if not m:
        m = re.search(r"vfon[:\s]+(\d+):(\d+):(\d+)", txt, re.IGNORECASE)
    if m:
        h, m_, s = int(m.group(1)), int(m.group(2)), int(m.group(3))
        return h * 3600 + m_ * 60 + s
    return None

def parse_duration_from_hea(hea_path: Path) -> float:
    """Try to parse approximate duration in hours from header first line (samples/fs)."""
    if not hea_path.exists():
        return None
    txt = hea_path.read_text()
    # e.g. "30 2 250 22099250 12:00:00"  --> 4th token often total samples
    m = re.search(r"^\s*\S+\s+\d+\s+(\d+)\s+(\d+)", txt, re.MULTILINE)
    if m:
        try:
            fs = float(m.group(1))
            nsamp = float(m.group(2))
            return nsamp / fs / 3600.0
        except Exception:
            pass
    return None

KNOWN_PACING = {
    "40": "continuous",
    "32": "intermittent",
    "43": "intermittent",
    "49": "intermittent",
    "51": "intermittent",
}  # from PhysioNet SDDB clinical table + research


def extract_and_clean(record: str, min_rr=250.0, max_rr=2000.0) -> dict:
    rec_path = str(SDDB_DIR / record)
    hea_path = SDDB_DIR / f"{record}.hea"
    atr_path = SDDB_DIR / f"{record}.atr"
    ari_path = SDDB_DIR / f"{record}.ari"

    # Accept .atr (preferred) or .ari (unaudited fallback); .hea always required
    if not hea_path.exists() or not (atr_path.exists() or ari_path.exists()):
        raise FileNotFoundError(f"Missing .hea or usable annotation (.atr/.ari) for {record} in {SDDB_DIR}")

    # Read annotations (beat times) - try atr (audited) then ari (unaudited) fallback
    ann = None
    annotator = "atr"
    try:
        ann = wfdb.rdann(rec_path, "atr")
    except Exception:
        try:
            ann = wfdb.rdann(rec_path, "ari")
            annotator = "ari"
            print(f"  {record}: .atr not found, using .ari (unaudited)")
        except Exception as e2:
            raise FileNotFoundError(f"No usable annotations (.atr or .ari) for {record}: {e2}")
    samples = np.asarray(ann.sample, dtype=float)
    fs = 250.0  # all SDDB are 250 Hz
    print(f"  {record}: using annotator={annotator}, n_samples={len(samples)}")

    if len(samples) < 10:
        raise ValueError(f"Too few annotations for {record}")

    # RR intervals in ms
    rr = np.diff(samples) / fs * 1000.0
    t_sec = samples[1:] / fs   # time of each RR (at second beat of interval)

    # vfon (robust)
    vfon_sec = parse_vfon_from_hea(hea_path)
    if vfon_sec is None:
        # fallback: use last annotation time
        vfon_sec = float(samples[-1] / fs)
        print(f"  WARNING: no vfon comment found for {record}, using end of anns")

    # duration heuristic from hea if available
    dur_h = parse_duration_from_hea(hea_path) or (float(t_sec[-1]) / 3600.0 if len(t_sec) > 0 else 0.0)

    # Clean: mark invalid RR
    valid = (rr >= min_rr) & (rr <= max_rr)
    n_invalid = int((~valid).sum())
    n_total = len(rr)
    interp_frac = float(n_invalid) / n_total if n_total > 0 else 0.0
    if n_invalid > 0:
        print(f"  {record}: {n_invalid}/{n_total} RR outside [{min_rr},{max_rr}] ms → interpolating ({100*interp_frac:.1f}%)")

    # Linear interp for invalid
    rr_clean = rr.copy().astype(float)
    if n_invalid > 0 and n_invalid < len(rr) * 0.3:  # safety: don't interp too many
        idx = np.arange(len(rr))
        rr_clean[~valid] = np.interp(idx[~valid], idx[valid], rr_clean[valid])
    else:
        # conservative: leave as-is or NaN (but keep for continuity in tau/recd)
        rr_clean[~valid] = np.nan
        # later code uses get_valid anyway, but for compatibility fill forward-ish
        # simple ffill
        mask = ~np.isfinite(rr_clean)
        rr_clean[mask] = np.interp(np.flatnonzero(mask),
                                   np.flatnonzero(~mask),
                                   rr_clean[~mask]) if (~mask).any() else rr_clean

    # Build t at RR positions (use original t for beats after interp)
    t_clean = t_sec.copy()

    # Compute vfon_beat_idx (closest to vfon_sec)
    vfon_beat_idx = int(np.argmin(np.abs(t_clean - vfon_sec)))

    total_h = float(t_clean[-1] / 3600.0) if len(t_clean) > 0 else 0.0

    # === Pacing detection (research-informed) ===
    txt_lower = hea_path.read_text().lower() if hea_path.exists() else ""
    pacing_comment = "pace" in txt_lower or "paced" in txt_lower
    known_pace = KNOWN_PACING.get(record, "none")
    # regularity heuristic on cleaned RR (paced often extremely regular)
    rr_finite = rr_clean[np.isfinite(rr_clean)]
    cv = (np.std(rr_finite) / np.mean(rr_finite)) if len(rr_finite) > 10 and np.mean(rr_finite) > 0 else 1.0
    regularity_paced = cv < 0.06  # very low CV suggests pacing or strong regularity
    pacing_detected = bool(pacing_comment or known_pace != "none" or regularity_paced)

    out = {
        "rr_ms": rr_clean.astype(np.float32),
        "t_sec": t_clean.astype(np.float32),
        "fs_ann": np.float32(fs),
        "vfon_sec": np.float32(vfon_sec),
        "vfon_beat_idx": np.int32(vfon_beat_idx),
        "n_beats": np.int32(len(rr_clean)),
        "total_hours": np.float32(total_h),
        # quality + metadata (new for batch reports)
        "n_invalid": np.int32(n_invalid),
        "interp_frac": np.float32(interp_frac),
        "pacing_detected": np.bool_(pacing_detected),
        "known_pacing_type": known_pace,
        "cv_rr": np.float32(cv),
    }
    return out

def save_npz(record: str, data: dict):
    out_path = DATA_DIR / f"rr_{record}_clean.npz"
    np.savez_compressed(out_path, **data)
    pace_str = data.get("known_pacing_type", "none")
    pdet = bool(data.get("pacing_detected", False))
    inv = int(data.get("n_invalid", 0))
    ifrac = 100 * float(data.get("interp_frac", 0))
    print(f"  Saved {out_path}  (n={data['n_beats']}, total_h={data['total_hours']:.2f})")
    print(f"    quality: invalid={inv} ({ifrac:.1f}%), pacing_detected={pdet} (type={pace_str}), cv={data.get('cv_rr',0):.4f}")
    return out_path

def load_selected_records():
    sel = BASE / "selected_records.txt"
    if sel.exists():
        return [line.strip() for line in sel.read_text().splitlines()
                if line.strip() and not line.strip().startswith("#")]
    inv = BASE / "data" / "records_inventory.csv"
    if inv.exists():
        import csv
        with inv.open() as f:
            return [r["record_id"] for r in csv.DictReader(f) if r.get("include?") == "Yes"]
    return ["30", "35"]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--record", default="", help="Record id or 'all'")
    ap.add_argument("--min-rr", type=float, default=250.0)
    ap.add_argument("--max-rr", type=float, default=2000.0)
    args = ap.parse_args()

    if args.record.lower() in ("all", ""):
        recs = load_selected_records()
    else:
        recs = [r.strip() for r in args.record.split(",") if r.strip()]

    print(f"Extracting RR for: {recs}")
    for rec in recs:
        print(f"\n=== {rec} ===")
        try:
            data = extract_and_clean(rec, min_rr=args.min_rr, max_rr=args.max_rr)
            save_npz(rec, data)
        except Exception as e:
            print(f"  FAILED for {rec}: {e}")

if __name__ == "__main__":
    main()
