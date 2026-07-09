# CCTP / SDDB: Systemic Tau and ordinal RECD before spontaneous VF

**Analysis code and companion data for the manuscript**

> *Context-Dependent Relational Reorganization of Heart Rate Dynamics Precedes Spontaneous Ventricular Fibrillation: Systemic Tau and Ordinal RECD Evidence from the Sudden Cardiac Death Holter Database*

**Author:** Johel Padilla-Villanueva, DrPH  
University of Puerto Rico, Medical Sciences Campus  
ORCID: [0000-0002-5797-6931](https://orcid.org/0000-0002-5797-6931)

**Companion libraries**

- [`systemictau`](https://github.com/johelpadilla/systemictau) — Systemic Tau core  
- Nested ordinal RECD levels are vendored in `code/recd_ordinal_levels.py` for exact paper reproducibility

**PhysioNet data**

- Sudden Cardiac Death Holter Database (SDDB)  
  DOI: [10.13026/C2W306](https://doi.org/10.13026/C2W306)  
  https://physionet.org/content/sddb/1.0.0/

---

## Repository contents

| Path | Description |
|------|-------------|
| `code/` | Full analysis pipeline used in the paper |
| `data/rr_*_clean.npz` | Cleaned RR series for the *N* = 10 analytic cohort |
| `data/sddb/` | PhysioNet headers + beat annotations (`.hea`, `.atr`/`.ari`; no large `.dat`) |
| `data/records_inventory.csv` | Record inventory / inclusion notes |
| `selected_records.txt` | Final *N* = 10 record list |
| `results/` | Batch summary CSV + per-record JSON outputs |
| `figures/` | Per-record and batch figures (including manuscript panels) |
| `manuscript/` | Manuscript source (Markdown + bibliography + embedded figures) |

### Analysis scripts (`code/`)

| Script | Role |
|--------|------|
| `download_sddb_records.py` | Download SDDB annotations from PhysioNet |
| `extract_rr.py` | RR extraction, cleaning, quality metadata |
| `analyze_cctp_pilot.py` | Systemic Tau + classical EWS + diagnostic figures |
| `run_cctp_surrogates.py` | Phase-shuffle surrogates for Δτ_s |
| `run_recd_on_rr.py` | Nested ordinal levels Φ₁–Φ₃ + excess3 |
| `run_recd_weighted_on_rr.py` | λ-weighted RECD contributions from \|τ_s\| |
| `run_cctp_batch.py` | Full batch orchestrator + comparative plots |
| `recd_ordinal_levels.py` | Nested RECD level math (vendored) |
| `_bootstrap.py` | Robust imports (`systemictau` + RECD) |

---

## Quick start

```bash
git clone https://github.com/johelpadilla/cctp-sddb-systemic-tau.git
cd cctp-sddb-systemic-tau

python3 -m venv .venv
source .venv/bin/activate   # Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

### Reproduce the paper batch (*N* = 10)

Cleaned RR files for records `30,31,32,35,36,38,45,47,50,51` are included. From the repo root:

```bash
python3 code/run_cctp_batch.py \
  --records 30,31,32,35,36,38,45,47,50,51 \
  --theta3 0.08 --high-thresh 0.65 --lambda-relative --force
```

Main consolidated table:

```text
results/cctp_batch_summary.csv
```

Batch figures:

```text
figures/batch/
```

### Single-record examples

```bash
python3 code/analyze_cctp_pilot.py --record 38
python3 code/run_cctp_surrogates.py --record 38
python3 code/run_recd_on_rr.py --record 38 --theta3 0.08 --high-thresh 0.65
python3 code/run_recd_weighted_on_rr.py --record 38 \
  --theta3 0.08 --high-thresh 0.65 --lambda-relative
```

### From raw PhysioNet annotations

```bash
python3 code/download_sddb_records.py --records 30,31,32,35,36,38,45,47,50,51
python3 code/extract_rr.py --record 38
```

If `systemictau` is installed from a local clone instead of PyPI:

```bash
export SYSTEMICTAU_SRC=/path/to/systemictau/src
```

---

## Analysis parameters (paper)

| Parameter | Value |
|-----------|--------|
| Analytic cohort | 30, 31, 32, 35, 36, 38, 45, 47, 50, 51 |
| Proxy | *X* = [z(RR), z(\|ΔRR\|)] |
| τ_s window / stride | *W*_τ = 101 beats, stride 5 |
| Ordinal RECD | *m* = 3, delay 1, *w*_φ = 101, stride 5 |
| θ₃ (re-calibrated) | 0.08 |
| High-level threshold | 0.65 |
| λ | record-relative \|τ_s\| / max\|τ_s\| |
| Surrogates | phase-shuffle, *n* = 8 per record |
| Primary Level-3 metric | continuous excess3 (not binary high-level rate) |

---

## Citation

If you use this code or derived results, please cite the manuscript and the Systemic Tau software release:

```bibtex
@software{Padilla2026CCTPcode,
  author  = {Padilla-Villanueva, Johel},
  title   = {{CCTP}/{SDDB}: Systemic Tau and ordinal {RECD} before spontaneous {VF}},
  year    = {2026},
  url     = {https://github.com/johelpadilla/cctp-sddb-systemic-tau},
  note    = {Analysis code accompanying the SDDB Holter manuscript}
}
```

Also cite PhysioNet SDDB ([Goldberger et al., 2000](https://doi.org/10.1161/01.CIR.101.23.e215); [Greenwald, 1986](https://dspace.mit.edu/handle/1721.1/28139)).

---

## License

- **Code** in this repository: MIT (see `LICENSE`)
- **PhysioNet SDDB** data: subject to PhysioNet / original database terms  
- **`systemictau`**: MIT (separate repository)

## Contact

Johel Padilla-Villanueva — [ORCID](https://orcid.org/0000-0002-5797-6931) · [github.com/johelpadilla](https://github.com/johelpadilla)
