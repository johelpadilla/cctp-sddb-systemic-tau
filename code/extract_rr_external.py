#!/usr/bin/env python3
"""
RR extraction for external Phase-1 databases (VFDB, CUDB, NSRDB).

Clean rules match discovery SDDB extract_rr.py: RR in [250, 2000] ms with
linear interpolation of short invalid gaps. No threshold retuning.

Sources
-------
- VFDB: rhythm-change annotations only → QRS via wfdb XQRS on ECG; event onset
  from first (VT|VFL|VF|ASYS) rhythm label.
- CUDB: beat annotations + rhythm markers; event onset from VF/VT-related aux.
- NSRDB: beat annotations only (negative-control Holters; no event).

Outputs cleaned npz under data/rr_external/{db}_{rec}_clean.npz
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import wfdb
from wfdb.processing import XQRS

BASE = Path(__file__).resolve().parent.parent
DATA = BASE / "data"
OUT_DIR = DATA / "rr_external"
OUT_DIR.mkdir(parents=True, exist_ok=True)

MIN_RR_MS = 250.0
MAX_RR_MS = 2000.0

# Malignant / terminal rhythm labels seen in VFDB/CU aux notes
MALIGNANT_PAT = re.compile(
    r"\b(VF|VFIB|VFL|V_FLUTTER|VT|VTACT|ASY|ASYS|VENTFIB)\b",
    re.IGNORECASE,
)


def _clean_aux(s: str) -> str:
    if s is None:
        return ""
    return str(s).replace("\x00", "").strip()


def _clean_rr(rr: np.ndarray, t_sec: np.ndarray) -> Tuple[np.ndarray, np.ndarray, int, float]:
    """Apply frozen [250, 2000] ms clean + linear interp (same spirit as extract_rr)."""
    rr = np.asarray(rr, dtype=float)
    t_sec = np.asarray(t_sec, dtype=float)
    valid = (rr >= MIN_RR_MS) & (rr <= MAX_RR_MS) & np.isfinite(rr)
    n_invalid = int((~valid).sum())
    n_total = len(rr)
    interp_frac = float(n_invalid) / n_total if n_total > 0 else 0.0
    rr_clean = rr.copy()
    if n_invalid > 0 and n_invalid < len(rr) * 0.5 and valid.any():
        idx = np.arange(len(rr))
        rr_clean[~valid] = np.interp(idx[~valid], idx[valid], rr_clean[valid])
    elif n_invalid > 0 and valid.any():
        mask = ~np.isfinite(rr_clean) | ~valid
        rr_clean[mask] = np.interp(
            np.flatnonzero(mask), np.flatnonzero(~mask), rr_clean[~mask]
        )
    return rr_clean.astype(np.float32), t_sec.astype(np.float32), n_invalid, interp_frac


def parse_malignant_onsets(
    samples: np.ndarray,
    symbols: List[str],
    aux_notes: Optional[List[str]],
    fs: float,
) -> List[Dict]:
    """Return list of malignant rhythm onsets {sample, sec, label}."""
    events = []
    aux_notes = aux_notes or [""] * len(samples)
    for i, (samp, sym) in enumerate(zip(samples, symbols)):
        aux = _clean_aux(aux_notes[i] if i < len(aux_notes) else "")
        label = None
        # MIT rhythm change annotations use '+' with aux like '(VF'
        if sym == "+" and aux:
            # strip leading '('
            body = aux.lstrip("(").upper()
            if MALIGNANT_PAT.search(body) or body.startswith(
                ("VF", "VT", "VFL", "ASY")
            ):
                label = body.split()[0] if body else aux
        # CU-style: sometimes '[' = VF start, '!' = ventricular flutter, etc.
        elif sym in ("[", "!", "*"):
            label = {"[": "VF", "!": "VFL", "*": "VT"}.get(sym, sym)
        elif aux and MALIGNANT_PAT.search(aux):
            label = _clean_aux(aux)
        if label:
            events.append(
                {
                    "sample": int(samp),
                    "sec": float(samp) / float(fs),
                    "label": label,
                    "symbol": sym,
                }
            )
    return events


def rr_from_beat_ann(
    record_path: str,
    annotator: str = "atr",
    fs_fallback: float = 250.0,
) -> Tuple[np.ndarray, np.ndarray, float, object]:
    """RR (ms) and t_sec from beat-level annotations; returns (rr, t, fs, ann)."""
    ann = wfdb.rdann(record_path, annotator)
    fs = float(ann.fs) if ann.fs else fs_fallback
    samples = np.asarray(ann.sample, dtype=float)
    # Keep only beat-like symbols if many non-beats present
    beat_syms = set("NLRBAaJSVeEjF/")  # common MIT beat codes + paced
    if ann.symbol is not None:
        mask = np.array(
            [s in beat_syms or s in ("N", "V", "S", "F", "Q", "/") for s in ann.symbol]
        )
        # If almost nothing matches (rhythm-only DB), fall through empty
        if mask.sum() >= 10:
            samples = samples[mask]
    if len(samples) < 10:
        raise ValueError(f"Too few beat annotations in {record_path} ({annotator})")
    rr = np.diff(samples) / fs * 1000.0
    t_sec = samples[1:] / fs
    return rr, t_sec, fs, ann


def rr_from_xqrs(record_path: str, channel: int = 0) -> Tuple[np.ndarray, np.ndarray, float]:
    """Detect QRS with XQRS and return RR (ms), t_sec, fs."""
    rec = wfdb.rdrecord(record_path, channels=[channel])
    fs = float(rec.fs)
    sig = np.asarray(rec.p_signal[:, 0], dtype=float)
    # replace nan
    if np.isnan(sig).any():
        nans = np.isnan(sig)
        sig[nans] = np.interp(np.flatnonzero(nans), np.flatnonzero(~nans), sig[~nans])
    xqrs = XQRS(sig=sig, fs=fs)
    xqrs.detect(verbose=False)
    peaks = np.asarray(xqrs.qrs_inds, dtype=float)
    if len(peaks) < 10:
        raise ValueError(f"XQRS found too few peaks for {record_path}")
    rr = np.diff(peaks) / fs * 1000.0
    t_sec = peaks[1:] / fs
    return rr, t_sec, fs


def extract_vfdb(rec: str, db_dir: Path = DATA / "vfdb") -> dict:
    """VFDB: XQRS RR + first malignant rhythm onset as event."""
    path = str(db_dir / rec)
    hea = db_dir / f"{rec}.hea"
    if not hea.exists():
        raise FileNotFoundError(hea)
    # rhythm annotations
    ann = wfdb.rdann(path, "atr")
    fs_ann = float(ann.fs) if ann.fs else 250.0
    events = parse_malignant_onsets(
        ann.sample, list(ann.symbol or []), list(ann.aux_note or []), fs_ann
    )
    if not events:
        # try any non-N rhythm
        events = []
        for i, aux in enumerate(ann.aux_note or []):
            a = _clean_aux(aux).lstrip("(").upper()
            if a and not a.startswith("N") and a not in ("NOISE", "NSR"):
                events.append(
                    {
                        "sample": int(ann.sample[i]),
                        "sec": float(ann.sample[i]) / fs_ann,
                        "label": a,
                        "symbol": ann.symbol[i] if ann.symbol else "+",
                    }
                )
    if not events:
        raise ValueError(f"No malignant/event onset found for VFDB {rec}")
    # first event after at least ~30s of recording preferred; else first
    first = events[0]
    for e in events:
        if e["sec"] >= 30.0:
            first = e
            break
    event_sec = float(first["sec"])
    event_label = str(first["label"])

    rr, t_sec, fs = rr_from_xqrs(path, channel=0)
    # keep only pre-event RR (+ tiny post for span bookkeeping); metrics use pre-event
    # store full RR series; event anchor separate
    rr_c, t_c, n_inv, ifrac = _clean_rr(rr, t_sec)
    total_h = float(t_c[-1] / 3600.0) if len(t_c) else 0.0
    pre_h = event_sec / 3600.0
    return {
        "rr_ms": rr_c,
        "t_sec": t_c,
        "fs_ann": np.float32(fs),
        "vfon_sec": np.float32(event_sec),
        "vfon_beat_idx": np.int32(int(np.argmin(np.abs(t_c - event_sec)))),
        "n_beats": np.int32(len(rr_c)),
        "total_hours": np.float32(total_h),
        "n_invalid": np.int32(n_inv),
        "interp_frac": np.float32(ifrac),
        "source_db": "vfdb",
        "record_id": rec,
        "event_label": event_label,
        "pre_event_hours": np.float32(pre_h),
        "rr_method": "xqrs",
    }


def extract_cudb(rec: str, db_dir: Path = DATA / "cudb") -> dict:
    """CUDB: prefer beat ann RR; fall back to XQRS; event from rhythm/VF markers."""
    path = str(db_dir / rec)
    hea = db_dir / f"{rec}.hea"
    if not hea.exists():
        raise FileNotFoundError(hea)
    ann = wfdb.rdann(path, "atr")
    fs_ann = float(ann.fs) if ann.fs else 250.0
    events = parse_malignant_onsets(
        ann.sample, list(ann.symbol or []), list(ann.aux_note or []), fs_ann
    )
    # CU often marks VF with '[' 
    if not events:
        for i, sym in enumerate(ann.symbol or []):
            if sym in ("[", "!", "*"):
                events.append(
                    {
                        "sample": int(ann.sample[i]),
                        "sec": float(ann.sample[i]) / fs_ann,
                        "label": {"[": "VF", "!": "VFL", "*": "VT"}[sym],
                        "symbol": sym,
                    }
                )
    if not events:
        raise ValueError(f"No VF/VT onset found for CUDB {rec}")
    first = events[0]
    for e in events:
        if e["sec"] >= 15.0:
            first = e
            break
    event_sec = float(first["sec"])
    event_label = str(first["label"])

    rr_method = "beat_ann"
    try:
        rr, t_sec, fs, _ = rr_from_beat_ann(path, "atr", fs_fallback=fs_ann)
        if len(rr) < 50:
            raise ValueError("sparse beats")
    except Exception:
        rr, t_sec, fs = rr_from_xqrs(path, channel=0)
        rr_method = "xqrs"

    rr_c, t_c, n_inv, ifrac = _clean_rr(rr, t_sec)
    total_h = float(t_c[-1] / 3600.0) if len(t_c) else 0.0
    return {
        "rr_ms": rr_c,
        "t_sec": t_c,
        "fs_ann": np.float32(fs),
        "vfon_sec": np.float32(event_sec),
        "vfon_beat_idx": np.int32(int(np.argmin(np.abs(t_c - event_sec)))),
        "n_beats": np.int32(len(rr_c)),
        "total_hours": np.float32(total_h),
        "n_invalid": np.int32(n_inv),
        "interp_frac": np.float32(ifrac),
        "source_db": "cudb",
        "record_id": rec,
        "event_label": event_label,
        "pre_event_hours": np.float32(event_sec / 3600.0),
        "rr_method": rr_method,
    }


def extract_nsrdb_control(rec: str, db_dir: Path = DATA / "nsrdb") -> dict:
    """NSRDB negative control: beat-ann RR, no event (vfon_sec = end of series)."""
    path = str(db_dir / rec)
    hea = db_dir / f"{rec}.hea"
    if not hea.exists():
        raise FileNotFoundError(hea)
    # header for fs
    try:
        header = wfdb.rdheader(path)
        fs_fb = float(header.fs)
    except Exception:
        fs_fb = 128.0
    rr, t_sec, fs, _ = rr_from_beat_ann(path, "atr", fs_fallback=fs_fb)
    rr_c, t_c, n_inv, ifrac = _clean_rr(rr, t_sec)
    total_h = float(t_c[-1] / 3600.0) if len(t_c) else 0.0
    # no VF — place synthetic end-anchor at end for schema compatibility
    return {
        "rr_ms": rr_c,
        "t_sec": t_c,
        "fs_ann": np.float32(fs),
        "vfon_sec": np.float32(t_c[-1] if len(t_c) else 0.0),
        "vfon_beat_idx": np.int32(len(rr_c) - 1 if len(rr_c) else 0),
        "n_beats": np.int32(len(rr_c)),
        "total_hours": np.float32(total_h),
        "n_invalid": np.int32(n_inv),
        "interp_frac": np.float32(ifrac),
        "source_db": "nsrdb",
        "record_id": rec,
        "event_label": "NONE_CONTROL",
        "pre_event_hours": np.float32(total_h),
        "rr_method": "beat_ann",
        "is_control": np.bool_(True),
    }


def save_npz(data: dict) -> Path:
    db = data.get("source_db", "ext")
    rec = data.get("record_id", "unk")
    out = OUT_DIR / f"{db}_{rec}_clean.npz"
    # numpy can't save plain str easily in all versions — store as 0-d object
    payload = {}
    for k, v in data.items():
        if isinstance(v, str):
            payload[k] = np.array(v)
        else:
            payload[k] = v
    np.savez_compressed(out, **payload)
    return out


def list_records(db_dir: Path) -> List[str]:
    recs_file = db_dir / "RECORDS"
    if recs_file.exists():
        return [ln.strip() for ln in recs_file.read_text().splitlines() if ln.strip()]
    return sorted({p.stem for p in db_dir.glob("*.hea")})


def main() -> int:
    ap = argparse.ArgumentParser(description="Extract cleaned RR for Phase-1 external DBs")
    ap.add_argument("--db", choices=["vfdb", "cudb", "nsrdb", "all"], default="all")
    ap.add_argument("--records", default="", help="comma list or empty=all local")
    args = ap.parse_args()

    jobs = []
    if args.db in ("vfdb", "all"):
        jobs.append(("vfdb", DATA / "vfdb", extract_vfdb))
    if args.db in ("cudb", "all"):
        jobs.append(("cudb", DATA / "cudb", extract_cudb))
    if args.db in ("nsrdb", "all"):
        jobs.append(("nsrdb", DATA / "nsrdb", extract_nsrdb_control))

    want = {r.strip() for r in args.records.split(",") if r.strip()} or None
    n_ok, n_fail = 0, 0
    for db_name, db_dir, fn in jobs:
        if not db_dir.exists():
            print(f"skip {db_name}: missing {db_dir}")
            continue
        recs = list_records(db_dir)
        if want:
            recs = [r for r in recs if r in want]
        print(f"\n=== {db_name}: {len(recs)} records ===")
        for rec in recs:
            # need at least hea; for xqrs need dat; for nsrdb atr is enough
            if not (db_dir / f"{rec}.hea").exists():
                print(f"  {rec}: missing hea")
                n_fail += 1
                continue
            try:
                data = fn(rec, db_dir)
                out = save_npz(data)
                print(
                    f"  {rec}: ok n={data['n_beats']} total_h={float(data['total_hours']):.3f} "
                    f"pre_h={float(data['pre_event_hours']):.3f} event={data.get('event_label')} "
                    f"method={data.get('rr_method')} inv={int(data['n_invalid'])} "
                    f"({100*float(data['interp_frac']):.1f}%) -> {out.name}"
                )
                n_ok += 1
            except Exception as e:
                print(f"  {rec}: FAIL {type(e).__name__}: {e}")
                n_fail += 1
    print(f"\nDone: ok={n_ok} fail={n_fail}")
    return 0 if n_ok > 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
