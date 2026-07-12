# CCTP / SDDB: Systemic Tau and ordinal RECD before spontaneous VF

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21270699.svg)](https://doi.org/10.5281/zenodo.21270699)

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
| `code/` | Full analysis pipeline used in the paper + Phase-1 external validation |
| `data/rr_*_clean.npz` | Cleaned RR series for the *N* = 10 analytic cohort (+ `rr_44` expansion) |
| `data/rr_external/` | Cleaned RR for VFDB / CUDB / NSRDB Phase-1 extracts |
| `data/sddb/` | PhysioNet headers + beat annotations (`.hea`, `.atr`/`.ari`; no large `.dat`) |
| `data/vfdb/`, `data/cudb/`, `data/nsrdb/` | External DB headers + annotations (`.dat` omitted; re-download from PhysioNet) |
| `data/records_inventory.csv` | Record inventory / inclusion notes |
| `selected_records.txt` | Final *N* = 10 record list |
| `results/` | Batch, stratified, lead-time, head-to-head, and Phase-1 external tables |
| `figures/` | Per-record, batch, and publication figures |
| `docs/` | External validation plan/report, Jul-12 interpretation, Copilot draft |
| `tests/` | Unit tests for detector / FAR / short-window helpers |
| `manuscript/` | Manuscript source (Markdown + PDF + bibliography + figures) |

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
| `cctp_metrics_core.py` | Pure lead-time / detector / stratification helpers |
| `run_sddb_inventory.py` | Full 23-record inventory + process status |
| `run_cohort_stratified.py` | Stratified cohort table (substrate / event type) |
| `run_leadtime_detector.py` | Lead-time + detector performance (τ_s, excess3, var, AR1) |
| `run_ews_head2head.py` | Head-to-head relational vs classic EWS table |
| `run_publication_figures.py` | Publication-ready stratification + lead-time figures |
| `extract_rr_external.py` | RR extract for VFDB (XQRS) / CUDB / NSRDB controls |
| `run_external_validation_phase1.py` | Frozen external sensitivity + control FAR (Phase 1) |

### Extension docs

| Path | Description |
|------|-------------|
| `docs/EXTERNAL_VALIDATION_PLAN.md` | Independent Holter targets, frozen metrics, success criteria |
| `docs/EXTERNAL_VALIDATION_PHASE1_REPORT.md` | Phase-1 results (VFDB n=11, NSRDB FAR; no clinical claim) |
| `docs/JUL12_RESULTS_INTERPRETATION.md` | Discovery stratified / lead-time / H2H reading |
| `docs/CLINICAL_COPILOT_DRAFT.md` | Copilot inputs/displays + research-use limits |
| `tests/test_leadtime_detector.py` | Unit tests for pure detector/lead-time functions |
| `tests/test_far_and_short_windows.py` | FAR / short-window helper tests |

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

### Full SDDB inventory, lead-time, EWS comparison (extension)

```bash
python3 code/run_sddb_inventory.py
python3 code/run_cohort_stratified.py
python3 code/run_leadtime_detector.py \
  --records 30,31,32,35,36,38,45,47,50,51 --theta3 0.08
python3 code/run_ews_head2head.py
python3 code/run_publication_figures.py
python3 tests/test_leadtime_detector.py
```

Key tables: `data/sddb_full_inventory.csv`, `results/cctp_cohort_stratified.csv`,
`results/leadtime_per_record.csv`, `results/ews_head2head.csv`,
`figures/publication/`.

### External Validation Phase 1 (frozen params; no retune)

Cleaned external RR under `data/rr_external/` is included. Raw VFDB/CU `.dat`
waveforms are **not** in the repo (re-download from PhysioNet if re-running XQRS).

```bash
# Re-run analysis on included cleaned RR (no PhysioNet download needed)
python3 code/run_external_validation_phase1.py
python3 -m pytest tests/test_leadtime_detector.py tests/test_far_and_short_windows.py -q
```

Outputs: `results/external_phase1_*.{csv,json}`, report
`docs/EXTERNAL_VALIDATION_PHASE1_REPORT.md`. **No clinical/deployability claim**
(NSRDB FAR ≫ clinical tolerance under the frozen abs-*z* rule).

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
  author    = {Padilla-Villanueva, Johel},
  title     = {{CCTP}/{SDDB}: Systemic Tau and ordinal {RECD} before spontaneous {VF}},
  year      = {2026},
  version   = {v1.0.1},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.21270699},
  url       = {https://github.com/johelpadilla/cctp-sddb-systemic-tau}
}
```

**DOIs:** version [10.5281/zenodo.21270699](https://doi.org/10.5281/zenodo.21270699) · concept [10.5281/zenodo.21270698](https://doi.org/10.5281/zenodo.21270698)

Also cite PhysioNet SDDB ([Goldberger et al., 2000](https://doi.org/10.1161/01.CIR.101.23.e215); [Greenwald, 1986](https://dspace.mit.edu/handle/1721.1/28139)).

---

## License

- **Code** in this repository: MIT (see `LICENSE`)
- **PhysioNet SDDB** data: subject to PhysioNet / original database terms  
- **`systemictau`**: MIT (separate repository)

## Contact

Johel Padilla-Villanueva — [ORCID](https://orcid.org/0000-0002-5797-6931) · [github.com/johelpadilla](https://github.com/johelpadilla)
