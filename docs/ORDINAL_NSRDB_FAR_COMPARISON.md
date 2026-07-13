# Ordinal vs abs-z FAR on NSRDB controls

*Generated: 2026-07-13T17:43:18.508962+00:00*

**Status:** Exploratory / interim public control FAR comparison. No clinical claim, no S5 claim, no superiority claim.

## Purpose

Measure False Alarm Rate (alarms per 24 h with refractory episode counting) for two native ordinal detectors and the frozen abs-z baseline on the same 18 NSRDB negative-control Holters, with Phase-2–aligned windowing.

Options remain **completely separate** (no OPC∧SDD fusion).

## Methodology (reproducible)

### Entry point

```bash
python code/run_ordinal_far_comparison.py --write-doc
```

### Data

- Controls: `data/rr_external/nsrdb_*_clean.npz` (18 files found, 18 processed).
- Cap: first `12.0` hours per record (Phase 2 default).
- Basal hours parameter: `2.0` (window = `(0.25, min(basal_hours, 0.25×total_h))`, same as Phase 1/2).
- Search: after basal end → end of capped recording.

### Symbol encoding (ordinal arms)

- Bivariate proxy from cleaned RR (project `build_bivariate_proxy`).
- Bandt–Pompe m=3, delay=1 → per-channel alphabet 6.
- Joint code σ = s0×6 + s1 ∈ {0,…,35}, **K=36** (same as bake-off).

### Detectors (fixed params; not retuned for FAR)

| Detector | Parameters |
|----------|------------|
| **OPC L=50** | L=50, θ_D=0.35, θ_R=5, K=36 |
| **SDD** | L_c=50, θ_TV=0.35, θ_S=1, fixed early basal, mask_basal=True |
| **abs-z τ_s** | z≥2.0, min_consecutive=3, W_TAU=101, stride=5 (**frozen**) |

### Episode counting and FAR

- **Refractory period:** 0.5 h after each episode.
- **abs-z:** `count_alarm_episodes` — abs-z ≥ 2 sustained ≥ 3 samples, then refractory.
- **OPC / SDD:** `count_binary_alarm_episodes` on the detector’s binary alarm stream (persistence already enforced by θ_R / θ_S); True outside refractory → episode.
- **FAR formula:** `total_episodes / total_search_hours × 24` via `false_alarm_rate`.
- **Pooled FAR** sums episodes and hours across records; **per-record FAR** = n_episodes_i / search_hours_i × 24.

### Comparability notes

- abs-z uses continuous strided τ_s timebase; ordinal uses symbol-endpoint times. Shared FAR *definition* enables comparison; sample indices are not identical.
- NSRDB is healthy Holter, **not** device-matched to VFDB/CU telemetry (same caveat as Phase 2).

## Results

### Pooled FAR

| Detector | Total episodes | Search hours | Pooled FAR (/24h) | Fraction records alarmed |
|----------|----------------|--------------|-------------------|--------------------------|
| OPC L=50 | 28 | 180.00 | **3.733** | 0.389 |
| SDD | 347 | 180.00 | **46.267** | 1.000 |
| abs-z τ_s | 253 | 180.00 | **33.734** | 1.000 |

### Per-record FAR statistics

| Detector | Mean | Median | Min | Max | Std |
|----------|------|--------|-----|-----|-----|
| OPC L=50 | 3.733 | 0.000 | 0.000 | 16.800 | 5.726 |
| SDD | 46.267 | 45.601 | 45.600 | 48.001 | 1.075 |
| abs-z τ_s | 33.734 | 33.600 | 19.200 | 45.601 | 7.042 |

### Per-record table

| Record | OPC ep | OPC FAR | SDD ep | SDD FAR | abs-z ep | abs-z FAR | Search h |
|--------|--------|---------|--------|---------|----------|-----------|----------|
| 16265 | 0 | 0.00 | 19 | 45.60 | 15 | 36.00 | 10.00 |
| 16272 | 0 | 0.00 | 19 | 45.60 | 16 | 38.40 | 10.00 |
| 16273 | 0 | 0.00 | 20 | 48.00 | 13 | 31.20 | 10.00 |
| 16420 | 6 | 14.40 | 19 | 45.60 | 10 | 24.00 | 10.00 |
| 16483 | 6 | 14.40 | 19 | 45.60 | 18 | 43.20 | 10.00 |
| 16539 | 0 | 0.00 | 20 | 48.00 | 14 | 33.60 | 10.00 |
| 16773 | 0 | 0.00 | 19 | 45.60 | 16 | 38.40 | 10.00 |
| 16786 | 0 | 0.00 | 20 | 48.00 | 12 | 28.80 | 10.00 |
| 16795 | 0 | 0.00 | 19 | 45.60 | 8 | 19.20 | 10.00 |
| 17052 | 4 | 9.60 | 19 | 45.60 | 12 | 28.80 | 10.00 |
| 17453 | 0 | 0.00 | 19 | 45.60 | 10 | 24.00 | 10.00 |
| 18177 | 0 | 0.00 | 20 | 48.00 | 19 | 45.60 | 10.00 |
| 18184 | 2 | 4.80 | 19 | 45.60 | 12 | 28.80 | 10.00 |
| 19088 | 2 | 4.80 | 19 | 45.60 | 14 | 33.60 | 10.00 |
| 19090 | 0 | 0.00 | 19 | 45.60 | 16 | 38.40 | 10.00 |
| 19093 | 0 | 0.00 | 19 | 45.60 | 14 | 33.60 | 10.00 |
| 19140 | 1 | 2.40 | 20 | 48.00 | 16 | 38.40 | 10.00 |
| 19830 | 7 | 16.80 | 19 | 45.60 | 18 | 43.20 | 10.00 |

## Qualitative observations

- Pooled FAR (alarms/24h): OPC L=50 = 3.733, SDD = 46.267, abs-z τ_s = 33.734 (n=18).
- OPC L=50 shows lower pooled FAR than frozen abs-z on these NSRDB Holters (observational; not a clinical superiority claim).
- SDD shows higher pooled FAR than frozen abs-z (more distribution-shift triggers on healthy Holter dynamics).
- Among ordinal rules, OPC L=50 appears more specific (lower FAR) than SDD under these fixed params — consistent with bake-off basal cleanliness.
- Timebases differ: abs-z on strided τ_s vs symbol-endpoint times for OPC/SDD; FAR remains comparable by shared definition (episodes/search_h×24), not identical indices.
- S5 (FAR ≤ 2/24h) is NOT claimed. NSRDB is not device-matched ICU/telemetry control.

## Specificity ranking (observational)

Rank by **lower pooled FAR** under these fixed params on NSRDB (not a clinical ranking):

1. **OPC L=50**: 3.733 / 24h
2. **abs-z τ_s**: 33.734 / 24h
3. **SDD**: 46.267 / 24h

### Recommendation (specificity-only, interim)

If the sole question is which **ordinal** rule looks more specific on healthy Holters under these defaults: prefer the lower-FAR ordinal arm above. Neither ordinal rule is promoted to production; abs-z remains the frozen baseline. Sensitivity on SDDB/VFDB was explored separately in the ordinal bake-off (`docs/ORDINAL_EXPLORATORY_BAKEOFF.md`). FAR and sensitivity trade-offs must be weighed jointly before any detector choice.

## Artifacts

- `results/ordinal_nsrdb_far_per_record.csv`
- `results/ordinal_nsrdb_far_by_detector.csv`
- `results/ordinal_nsrdb_far_summary.json`
- Runner: `code/run_ordinal_far_comparison.py`
- Episode helper: `count_binary_alarm_episodes` in `code/cctp_metrics_core.py`

## Non-claims

- No S5 (FAR ≤ 2/24h).
- No clinical utility / FDA / deployability.
- No fusion of OPC and SDD.
- No retune of abs-z or ordinal thresholds based on this run.
