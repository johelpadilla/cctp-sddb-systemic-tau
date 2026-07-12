# Jul-12 2026 results interpretation (SDDB N=10)

**Date:** 2026-07-12  
**Source of truth:** `results/cctp_*stratified*`, `results/leadtime_*`, `results/ews_head2head*`  
**Decision:** **Option A** — fold stratified / lead-time / head-to-head into the manuscript (not external validation or Copilot prototype this session).

Frozen parameters: $\theta_3=0.08$, high-threshold $=0.65$, $W_{\tau}=101$, $W_{\mathrm{EWS}}=501$, stride $=5$, relative $\lambda$, detector rule abs-$z\geq 2$ for $\geq 3$ consecutive windows.

---

## 1. Corrected event geometry

After fixing empty `event_hr` (npz timing via `resolve_event_timing_from_npz` / inventory):

| Label | n | Records |
|-------|---|---------|
| intermediate | **6** | 30, 32, 38, 45, 47, 50 |
| terminal | **4** | 31, 35, 36, 51 |

All 10 rows have `event_hr_source=npz`. Substrates: sinus 5, AF 3, paced 2.

---

## 2. Stratified deltas (approach − basal)

From `cctp_stratified_summary.csv` (means):

| Stratum | n | mean Δτ_s | mean Δexcess3 | mean Δvar | mean ΔAR(1) |
|---------|---|-----------|---------------|-----------|-------------|
| AF | 3 | **+0.056** | +0.0028 | +11269 | +0.013 |
| paced | 2 | **−0.040** | **−0.021** | −446 | +0.036 |
| sinus | 5 | +0.013 | +0.008 | +12256 | −0.150 |
| intermediate | 6 | +0.019 | +0.009 | +15766 | −0.097 |
| terminal | 4 | +0.009 | −0.011 | −101 | −0.013 |

**Reading:** paced substrate consistently tightens/loosens in the **negative** relational direction (both τ_s and excess3). AF tends **positive** for τ_s. Sinus is mixed (median Δτ_s ≈ 0). Intermediate events show larger mean variance rises; terminal group-mean excess3 is slightly negative — do **not** over-claim a universal terminal signature with n=4.

---

## 3. Lead-time detector

Alarm: first sustained abs-z departure of the metric series from its basal mean/SD (`z_threshold=2`, `min_consecutive=3`). Same basal/approach windows as the manuscript pipeline.

| Metric | Sensitivity | Median lead (h) | Mean lead (h) | FAR |
|--------|-------------|-----------------|---------------|-----|
| excess3 | 1.0 | **6.86** | 5.91 | **nan** |
| τ_s | 1.0 | **5.88** | 5.63 | nan |
| AR(1) | 1.0 | 5.97 | 5.19 | nan |
| variance | 1.0 | **3.90** | 4.22 | nan |

Cumulative detection rate (fraction with lead ≥ horizon):

| Horizon | excess3 | τ_s | AR(1) | var |
|---------|---------|-----|-------|-----|
| 2 h | 1.0 | 1.0 | 1.0 | 0.9 |
| 3 h | 0.9 | 0.9 | 0.9 | 0.7 |
| 6 h | **0.7** | 0.4 | 0.5 | **0.2** |

**Reading:** On this discovery cohort, relational metrics (especially excess3) alarm earlier in the median and retain higher detection at 6 h than variance. **Sensitivity 1.0 is not clinically deployable** — every record is a true event case; FAR is undefined without control Holters. Shortest τ_s lead is record 47 (~2.47 h), a borderline discordant case.

---

## 4. Head-to-head concordance (sign of Δ)

| Pair | Concordance |
|------|-------------|
| τ_s vs excess3 | **8/10 = 0.8** |
| τ_s vs AR(1) | 7/10 = 0.7 |
| τ_s vs var | 5/10 = 0.5 |
| excess3 vs var | 5/10 = 0.5 |
| var vs AR(1) | 2/10 = 0.2 |

Effect-size signs: τ_s and excess3 each positive in 5/10; var and AR(1) positive in 6/10 (AR(1) “positive” = rise, which is uncommon under classical CSD expectation that AR rises; empirically AR often falls).

**Reading:** Relational pair (τ_s, excess3) is the only high-concordance pair. Classic var and AR(1) are poorly aligned with each other (0.2), underscoring that pre-VF Holter dynamics are not a single CSD story.

Discordant relational records remain **47** and **50** (small / borderline effect sizes).

---

## 5. Strengths, limitations, paper value

### Strengths
1. Correct event geometry unlocks honest intermediate vs terminal strata (was previously terminal-heavy by bug).
2. Multi-hour median lead times for excess3/τ_s with frozen abs-z rule — useful **hypothesis** for online monitors.
3. Direction concordance 0.8 between independent relational constructions (Kendall τ_s vs nested ordinal excess3).
4. Context-dependence by substrate is now table-ready (AF vs paced polarity).
5. Publication figures exist: stratified deltas, lead-time, concordance.

### Limitations (must stay explicit)
1. **N=10 discovery** on a single public corpus (SDDB ≤23).
2. **FAR = nan** — no control Holters; sensitivity 1.0 is optimistic by design.
3. abs-z vs basal may fire on non-specific non-stationarity; not a validated clinical alarm.
4. Borderline 47/50 still discordance; group terminal excess3 mean is fragile (n=4).
5. Thresholds frozen from light re-cal; not re-tuned here.

### Paper value (what Jul-12 adds beyond the N=10 batch table)
- Moves from “sign concordance exists” to **how early** a simple detector would have fired (hours, not just basal/approach contrast).
- Places classical EWS in a **concordance matrix**, not only narrative comparison.
- Stratifies by **substrate + event_type** after correct VF anchors — supports the “context-dependent reorganization” thesis with numbers.

---

## 6. Decision and next session

| Option | Verdict this session |
|--------|----------------------|
| **A. Manuscript update** | **EXECUTED** — highest impact; figures/tables ready; paper was the bottleneck. |
| B. External validation (VFDB/CU) | Deferred — plan already in `docs/EXTERNAL_VALIDATION_PLAN.md`; kick off next. |
| C. Clinical Copilot prototype | Deferred — draft only in `docs/CLINICAL_COPILOT_DRAFT.md`. |

### Next-session prompt (copy-paste)

```
Continue CCTP after Jul-12 manuscript update (stratified + lead-time + H2H in
manuscript/CCTP_SDBB_manuscript.md and docs/JUL12_RESULTS_INTERPRETATION.md).

Highest-impact next work: EXTERNAL VALIDATION Phase 1.
- Follow docs/EXTERNAL_VALIDATION_PLAN.md
- PhysioNet VFDB and/or CU ventricular tachyarrhythmia DB + non-event control Holters
- Reuse frozen params (θ3=0.08, high-thresh=0.65, W_TAU=101, abs-z detector)
- Primary deliverable: control-arm FAR estimate + external sensitivity table
- Do NOT retune thresholds; do NOT claim clinical deployability
Root: Investigaciones/Cardiac_CCTP_Pilot
```
