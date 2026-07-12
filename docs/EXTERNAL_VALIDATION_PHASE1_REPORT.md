# External Validation Phase 1 — Progress Report

**Date:** 2026-07-12  
**Plan of record:** `docs/EXTERNAL_VALIDATION_PLAN.md`  
**Status:** Phase 1 **started and reported** (not a completed multi-center validation)  
**Clinical / FDA / deployability claim:** **NONE**

Frozen discovery parameters (Jul-12 2026; **not retuned** on validation data):

| Parameter | Value |
|-----------|--------|
| θ₃ | 0.08 |
| high-threshold (RECD reporting) | 0.65 |
| W_TAU | 101 |
| stride | 5 |
| W_EWS (comparators) | 501 |
| Detector | abs-z ≥ 2.0 sustained ≥ 3 consecutive windows |
| RR clean | [250, 2000] ms |

---

## 1. Data footprint

| Source | Local status | Role |
|--------|--------------|------|
| **PhysioNet VFDB** | 22/22 records (hea/dat/atr) under `data/vfdb/` | Independent short pre-VF/VT episodes |
| **PhysioNet CUDB** | 35/35 records under `data/cudb/` | Independent; almost all **too short** for ≥15 min pre-event |
| **PhysioNet NSRDB** | 6 Holters (hea+atr) under `data/nsrdb/` | Negative controls (no VF) for FAR |
| **SDDB rec 44** | Existing `data/rr_44_clean.npz` | Optional **internal extension** (not independent) |
| Cleaned RR | `data/rr_external/*_clean.npz` (62 files) | VFDB XQRS; CU/NSRDB beat ann |

Inventory with inclusion/exclusion: `results/external_phase1_inventory.csv` (63 rows including SDDB extension when enabled).

### Inclusion summary (events)

| DB | Included | Excluded (main reason) |
|----|----------|------------------------|
| VFDB | **11** | 11: pre-event &lt; 15 min and/or high interp |
| CUDB | **0** | All: pre-event &lt; 15 min (episodes ~8 min total) |
| SDDB 44 | 1 (extension) | — |

**Independent processable event n = 11** (all VFDB, duration stratum `short_15_60min`). Target n≥20 not reached from public short DBs under the ≥15 min pre-event floor — reported honestly, not invented.

---

## 2. Pipeline applied (frozen)

Entry points:

```bash
python3 code/extract_rr_external.py --db all
python3 code/run_external_validation_phase1.py
# optional internal extension:
python3 code/run_external_validation_phase1.py --include-sddb-extension
```

Core math reused from `code/cctp_metrics_core.py` + `recd_ordinal_levels.py` + `systemictau`. Short-DB windows: pre-event thirds (basal / buffer / approach) via `short_db_windows`. Lead-time via `detect_lead_time` (same abs-z rule as discovery). FAR via `count_alarm_episodes` + `false_alarm_rate` on NSRDB.

Artifacts:

| File | Content |
|------|---------|
| `results/external_phase1_inventory.csv` | Inclusion/exclusion + strata |
| `results/external_phase1_per_record.csv` | Per-record lead-time + deltas |
| `results/external_phase1_sensitivity.json` | Sensitivity / lead-time by metric |
| `results/external_phase1_controls.csv` | Control episode counts |
| `results/external_phase1_far.json` | FAR / 24 h |
| `results/external_phase1_summary.json` | Master summary (clinical_claim=false) |

---

## 3. Independent sensitivity & lead-time (VFDB n=11)

| Metric | n | Detected | Sensitivity | Median lead (h) | Median lead (min) |
|--------|---|----------|-------------|-----------------|-------------------|
| **τ_s** | 11 | 11 | **1.00** | 0.173 | ~10.4 |
| **excess3** | 11 | 9 | **0.82** | 0.127 | ~7.6 |
| AR(1) | 11 | 6 | 0.55 | 0.236 | ~14.2 |
| variance | 11 | 6 | 0.55 | 0.236 | ~14.1 |

Sign concordance τ_s vs excess3 (Δ approach−basal): **7/11 = 0.64**.

### vs pre-specified success bars (S1–S6) — **not claimed as pass/fail package**

| Criterion | Preliminary read |
|-----------|------------------|
| **S3 Lead time** (short DB ≥15 min median) | **Not met** — median leads ~8–10 min |
| **S4 Detector sensitivity ≥0.60** | **Met on this n=11 set** for τ_s and excess3 |
| **S2 Concordance ≥0.70** | **Not met** (0.64) on this short set |
| **S5 FAR ≤2 / 24 h** | **Failed hard** on NSRDB controls (below) |
| **S1 / S6** | Not fully tested (no phase-shuffle suite; exploratory H2H only) |

Sensitivity 1.0 / 0.82 on a **tiny short-recording set** must **not** be read as clinical performance.

---

## 4. False-alarm rate (NSRDB controls)

| Metric | n controls | Search hours | Episodes | **FAR / 24 h** | Fraction alarmed |
|--------|------------|--------------|----------|----------------|------------------|
| τ_s | 6 | ~60 | 86 | **~34.4** | 1.0 |
| excess3 | 6 | ~60 | 72 | **~28.8** | 1.0 |

**Interpretation (honest):** The frozen abs-z rule is **not specific** on healthy long Holters. Every control produced ≥1 episode. S5 (≤2 alarms/24 h) is **not** approached. This was the main discovery gap (FAR=nan on SDDB N=10) and is now quantified — **detrimental for deployability**, informative for the paper.

Caveats:

- NSRDB is rhythm-healthy Holter, **not** device-matched to VFDB telemetry.
- Basal = early ~2 h; search capped at 12 h/record; refractory 0.5 h between episodes.
- High FAR may partly reflect non-stationarity / diurnal HRV vs a short basal reference — not “VF precursors.”

FAR is **numeric and honest**, not nan: controls were available. Specificity for clinical use remains **unproven / currently poor** under this detector rule.

---

## 5. Strengths, limitations, paper value

### Strengths
1. Independent public sources downloaded and inventoried with exclusion reasons.
2. Frozen params echoed in every Phase-1 JSON artifact — no retuning.
3. First quantitative **control-arm FAR** for the Jul-12 detector rule.
4. External short-DB sensitivity table for τ_s / excess3 (n=11).
5. Reproducible CLI + unit tests for FAR/short-window helpers.

### Limitations
1. n=11 independent events (target ≥20 not met); all short_15_60min.
2. CUDB contributes **zero** processable events under ≥15 min pre-event floor.
3. FAR far above clinical tolerance; detector not deployable.
4. VFDB RR from XQRS (rhythm-only annotations) — QRS error possible.
5. Lead times are minutes, not multi-hour SDDB discovery leads — strata not comparable without care.
6. No institutional Holters / IRB Phase 2 yet.

### Value for the research program
Phase 1 **does not green-light** a clinical Copilot. It **does** show: (a) relational metrics can still fire before short public VT/VF onsets under frozen rules; (b) the same rule floods healthy Holters with alarms — so the next scientific bottleneck is **specificity / multi-metric fusion / calibrated basals**, not more discovery sensitivity on SDDB alone.

---

## 6. Recommended next step

1. **Phase 2 (highest impact):** Institutional or device-matched non-event Holters + longer pre-VF recordings; keep params frozen; redesign detector only under a **pre-registered** change protocol if S5 remains impossible.  
2. Optional: exploratory multi-metric fusion **labeled exploratory** (not primary frozen claim).  
3. Do **not** ship Copilot alarms with current FAR.  
4. Optional manuscript note: add a short “External pilot (Phase 1)” paragraph citing these tables without overclaim.

### Next-session prompt

```
Continue CCTP External Validation after Phase 1
(docs/EXTERNAL_VALIDATION_PHASE1_REPORT.md).

Facts: VFDB n=11 independent short events (τs sens=1.0, excess3=0.82,
median lead ~8–10 min); NSRDB FAR ~29–34 / 24h (S5 fail); CU all excluded
(pre<15min); params frozen θ3=0.08 W_TAU=101 abs-z≥2×3.

Highest-impact next: Phase 2 specificity — institutional/matched controls
and/or pre-registered detector fusion WITHOUT retuning on SDDB discovery.
Do not claim clinical deployability.
Root: Investigaciones/Cardiac_CCTP_Pilot
```

---

## 7. Code & tests

- `code/extract_rr_external.py` — VFDB/CU/NSRDB RR extraction  
- `code/run_external_validation_phase1.py` — Phase 1 analysis entry  
- `code/cctp_metrics_core.py` — `short_db_windows`, `count_alarm_episodes`, `false_alarm_rate`, frozen constants  
- `tests/test_far_and_short_windows.py` + `tests/test_leadtime_detector.py`
