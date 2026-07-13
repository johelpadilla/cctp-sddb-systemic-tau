# Ordinal vs abs-z: sensitivity–specificity trade-off (exploratory)

*Generated: 2026-07-13T17:54:41.456235+00:00*

**Status:** Exploratory join of existing bake-off + NSRDB FAR results. **No clinical claim, no S5 claim, no superiority claim, no fusion.**

## Executive summary

Exploratory trade-off on public PhysioNet cohorts (SDDB n=11, VFDB n=22 events; NSRDB n=18 controls, ~180 search-hours, refractory 0.5 h). Three detectors kept strictly separate with fixed parameters — no fusion, no production retune, no clinical/S5/FDA claim. OPC L=50 is the most specific arm (FAR≈3.73/24h, 28 episodes) but detects only 14/33 events (sens≈0.42). SDD is the most sensitive (sens≈0.97) with the highest control FAR (≈46.3/24h, 347 episodes). Frozen abs-z is intermediate on FAR (≈33.7/24h) with high sensitivity (≈0.91). Which arm looks 'balanced' depends on the objective: specificity → OPC L=50; hit rate → SDD; frozen pilot baseline → abs-z; sens-per-FAR efficiency → OPC L=50. NSRDB is not device-matched to VFDB/CU telemetry.

## Inputs (no reprocessing)

| Artifact | Role |
|----------|------|
| `results/ordinal_exploratory_summary.json` | existing |
| `results/ordinal_nsrdb_far_summary.json` | existing |
| `results/ordinal_nsrdb_far_per_record.csv` | existing |
| `results/ordinal_opc_sdd_absz_comparison.csv` | existing |
| `docs/ORDINAL_EXPLORATORY_BAKEOFF.md` | existing |
| `docs/ORDINAL_NSRDB_FAR_COMPARISON.md` | existing |

## Primary arms (strictly separate)

| Detector | Sensitivity source key | FAR source key |
|----------|------------------------|----------------|
| OPC L=50 | `opc_L50_companion` | `opc_L50` |
| SDD (TV) | `sdd` | `sdd` |
| abs-z τ_s frozen | `absz_tau_s` | `absz_tau_s` |

OPC L=8 is **not** a primary arm (K/L saturation with joint K=36).

## Trade-off table (sensitivity vs FAR)

| Detector | Sens SDDB | Sens VFDB | Sens all | FAR /24h | Control episodes | Frac. controls alarmed | Sens/FAR |
|----------|-----------|-----------|----------|----------|-------------------|------------------------|----------|
| OPC L=50 | 0.545 (6/11) | 0.364 (8/22) | 0.424 (14/33) | **3.733** | 28 | 0.389 | 0.1136 |
| SDD (TV) | 1.000 (11/11) | 0.955 (21/22) | 0.970 (32/33) | **46.267** | 347 | 1.000 | 0.0210 |
| abs-z τ_s (frozen) | 1.000 (11/11) | 0.864 (19/22) | 0.909 (30/33) | **33.734** | 253 | 1.000 | 0.0269 |

## Cohort summary

| Detector | Cohort | Role | n | Detected / episodes | Sensitivity | FAR /24h |
|----------|--------|------|---|---------------------|-------------|----------|
| OPC L=50 | SDDB | events_sensitivity | 11 | 6 | 0.545 | — |
| OPC L=50 | VFDB | events_sensitivity | 22 | 8 | 0.364 | — |
| OPC L=50 | all_events | events_sensitivity | 33 | 14 | 0.424 | — |
| OPC L=50 | NSRDB | controls_far | 18 | 28 | — | 3.733 |
| SDD (TV) | SDDB | events_sensitivity | 11 | 11 | 1.000 | — |
| SDD (TV) | VFDB | events_sensitivity | 22 | 21 | 0.955 | — |
| SDD (TV) | all_events | events_sensitivity | 33 | 32 | 0.970 | — |
| SDD (TV) | NSRDB | controls_far | 18 | 347 | — | 46.267 |
| abs-z τ_s (frozen) | SDDB | events_sensitivity | 11 | 11 | 1.000 | — |
| abs-z τ_s (frozen) | VFDB | events_sensitivity | 22 | 19 | 0.864 | — |
| abs-z τ_s (frozen) | all_events | events_sensitivity | 33 | 30 | 0.909 | — |
| abs-z τ_s (frozen) | NSRDB | controls_far | 18 | 253 | — | 33.734 |

## Observational rankings (multi-criterion)

### best_specificity_lowest_far (lower better)

1. **OPC L=50**: 3.73338
2. **abs-z τ_s (frozen)**: 33.7337
3. **SDD (TV)**: 46.2672

### best_sensitivity_all_events (higher better)

1. **SDD (TV)**: 0.969697
2. **abs-z τ_s (frozen)**: 0.909091
3. **OPC L=50**: 0.424242

### best_sensitivity_sddb (higher better)

1. **SDD (TV)**: 1
2. **abs-z τ_s (frozen)**: 1
3. **OPC L=50**: 0.545455

### best_sensitivity_vfdb (higher better)

1. **SDD (TV)**: 0.954545
2. **abs-z τ_s (frozen)**: 0.863636
3. **OPC L=50**: 0.363636

### best_sens_per_far_unit (higher better)

1. **OPC L=50**: 0.113635
2. **abs-z τ_s (frozen)**: 0.026949
3. **SDD (TV)**: 0.0209586

### best_geometric_balance (higher better)

1. **abs-z τ_s (frozen)**: 0.614744
2. **OPC L=50**: 0.605914
3. **SDD (TV)**: 0.575503

### best_youden_like_exploratory (higher better)

1. **abs-z τ_s (frozen)**: 0.324792
2. **SDD (TV)**: 0.31125
3. **OPC L=50**: 0.289626

### fewest_control_episodes (lower better)

1. **OPC L=50**: 28
2. **abs-z τ_s (frozen)**: 253
3. **SDD (TV)**: 347

## Strengths and weaknesses (qualitative)

### OPC L=50

*Operating point:* specificity-leaning / selective collapse detector

**Strengths**
- Lowest pooled FAR on NSRDB among the three primary arms under fixed params.
- Majority of NSRDB records have zero OPC episodes (median per-record FAR = 0).
- Collapse+persistence predicate is conceptually aligned with RECD ordinal structure.
- Cleaner basal behavior on SDDB in the bake-off vs saturated L=8.

**Weaknesses**
- Lowest exploratory event sensitivity (SDDB 0.55, VFDB 0.36, all 0.42).
- Misses several Holters that SDD and abs-z detect.
- Still above informal S5 FAR target of ≤2/24h on this non-device-matched set (no S5 claim).
- Sensitive to L vs K alphabet interaction (L=8 is invalid for K=36).

### SDD (TV)

*Operating point:* sensitivity-leaning / distribution-shift detector

**Strengths**
- Highest exploratory event sensitivity among the three (SDDB 1.00, VFDB 0.95, all 0.97).
- Detects distributional reorganization (TV) rather than amplitude excursions.
- Matches or slightly exceeds abs-z hit rate on these public event sets.

**Weaknesses**
- Highest pooled FAR on NSRDB (~46.3/24h) — near refractory-period ceiling on long Holters.
- Fraction of controls alarmed = 1.0 under these fixed params.
- Fixed early basal on ambulatory Holter dynamics is a hard specificity stress test.
- High sensitivity alone is not a usable operating point if control FAR saturates.

### abs-z τ_s (frozen)

*Operating point:* frozen continuous baseline; intermediate FAR vs high sens

**Strengths**
- Frozen pilot baseline; parameters unchanged from production discovery freeze.
- High exploratory sensitivity (SDDB 1.00, VFDB 0.86, all 0.91).
- Intermediate FAR between OPC L=50 and SDD on NSRDB (~33.7/24h).
- Phase-2 FAR order of magnitude reproduced (sanity anchor).

**Weaknesses**
- FAR remains high on healthy Holter relative to informal S5 (no S5 claim).
- All NSRDB controls alarmed at least once under z≥2 sustained ≥3.
- Continuous amplitude z-score is a different predicate than ordinal collapse/divergence.
- Three VFDB misses are in too_short pre-event strata.

## Balance by objective

- **if_priority_is_lowest_false_alarm_burden** → observational preference: `opc_L50` — OPC L=50 pooled FAR 3.733/24h vs abs-z 33.734 and SDD 46.267; accepts substantially lower event hit rate.
- **if_priority_is_maximum_event_hit_rate** → observational preference: `sdd` — SDD all-events sensitivity 0.970 (32/33) edges abs-z 0.909 (30/33); FAR cost is highest of the three under these fixed params.
- **if_priority_is_sens_per_unit_far** → observational preference: `opc_L50` — Among fixed-param arms, OPC L=50 yields the highest sens_all / FAR ratio because FAR is ~9× lower than abs-z while retaining partial event detection.
- **if_priority_is_frozen_baseline_continuity** → observational preference: `absz_tau_s` — abs-z remains the only production-frozen continuous baseline; ordinal arms are methodological alternatives, not retunes of production thresholds.
- **if_priority_is_balanced_hit_rate_with_intermediate_far** → observational preference: `absz_tau_s` — abs-z sits between OPC (selective, low FAR) and SDD (near-ceiling FAR) with high event sensitivity; geometric/Youden-like proxies may still favor OPC because FAR dominates the denominator — interpret proxies cautiously.

## Recommendations (non-binding)

### A. Keep arms separate; do not fuse as primary claim

Continue evaluating OPC L=50, SDD, and abs-z independently. Any light fusion (e.g. sequential filter) would be a new experiment with its own multiplicity cost — not a free upgrade.

### B. If pursuing OPC, explore parameter/alphabet alignment — not abs-z retune

L vs K interaction is material (L=8 invalid for K=36). Modest grid on L, θ_D, θ_R on a held-out or bootstrap scheme could move the operating point without claiming clinical optimization. Do not retune production abs-z.

### C. If pursuing SDD, redesign basal/control policy before FAR claims

Current fixed early basal yields near-ceiling FAR on ambulatory Holters. Sliding basal, longer L_c, higher θ_TV, or θ_S>1 are natural specificity levers — each requires a pre-registered comparison, not post-hoc cherry-picking.

### D. Institutional device-matched controls (Tier A) for true FAR

Public NSRDB remains an interim upper-bound flavor only. Phase-2/S5-style specificity still needs partner Holter/telemetry matched to deployment setting.

### E. Optional light cascade (exploratory only, if ever tested)

A non-primary cascade such as 'SDD candidate → OPC confirm' or vice versa could be scored on the existing per-record tables without reprocessing RR, but only as a clearly labeled secondary analysis with no fusion claim in the primary manuscript arm.

### F. Manuscript framing

Present the three operating points as a trade-off surface, not a winner. Emphasize predicate differences (collapse vs distributional TV vs continuous z).

## Caveats / non-claims

- Exploratory public-data bake-off only; small event n (11+22).
- FAR ≠ classical 1−specificity; episode rate with 0.5 h refractory on Holter.
- NSRDB not device-matched to VFDB/CU; not institutional ICU telemetry.
- Timebases differ (strided τ_s vs symbol endpoints); FAR definition is shared.
- No threshold optimization / validation-set retune in this join analysis.
- No clinical utility, FDA readiness, or deployability is asserted; S5 (FAR ≤ 2/24h) is not asserted.
- OPC L=8 excluded as primary arm (parameter–alphabet saturation).
- No OPC∧SDD fusion scored.

## Artifacts

- `results/ordinal_sensitivity_specificity_tradeoff.csv`
- `results/ordinal_tradeoff_by_cohort.csv`
- `results/ordinal_tradeoff_summary.json`
- `docs/ORDINAL_SENSITIVITY_SPECIFICITY_TRADEOFF.md`
- Runner: `code/run_ordinal_tradeoff_analysis.py`
