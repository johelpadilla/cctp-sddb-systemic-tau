# External Validation Plan — Systemic Tau + RECD before spontaneous VF

**Status:** research plan; **Phase 1 reported 2026-07-12** (`docs/EXTERNAL_VALIDATION_PHASE1_REPORT.md`); **Phase 2 planning + interim full-NSRDB FAR (n=18) 2026-07-12** (`docs/EXTERNAL_VALIDATION_PHASE2_PLAN.md`, `docs/EXTERNAL_VALIDATION_PHASE2_PROGRESS.md`) — institutional Tier A still open; **no clinical claim**; **S5 not claimed**  
**Primary metrics (frozen):** Δτ_s and continuous Δexcess3 on bivariate RR proxy  
**Paper parameters (do not re-tune on validation sets):** θ₃ = 0.08, high-thresh = 0.65, W_τ = 101, stride = 5, λ-relative for weighted RECD (secondary)

---

## 1. Objective

Test whether relational reorganization detected by Systemic Tau (τ_s) and nested ordinal RECD excess3 appears **before spontaneous ventricular fibrillation (VF)** in **independent** Holter / long-term ECG sources outside the PhysioNet Sudden Cardiac Death Holter Database (SDDB) discovery cohort (N = 10 analytic records).

Success is defined prospectively below; no threshold re-optimization on validation data is permitted for primary endpoints.

---

## 2. Target independent data sources

| Priority | Source | Rationale | Access notes |
|----------|--------|-----------|--------------|
| 1 | **PhysioNet MIT-BIH Malignant Ventricular Arrhythmia Database (VFDB)** | Annotated malignant VT/VF episodes; short but public | Open PhysioNet |
| 2 | **PhysioNet CU Ventricular Tachyarrhythmia Database** | VF/VT terminations; public | Open PhysioNet |
| 3 | **PhysioNet Sudden Cardiac Death Holter Database — held-out records** | Remaining SDDB records with usable annotations not in N=10 (e.g. 33, 43, 44 if extractable) | Same license as discovery; treat as **internal extension**, not fully independent |
| 4 | **Institutional Holter archives (ICD / telemetry pre-arrest)** with timed VF or sustained VT | True external clinical validation | IRB + de-identification required |
| 5 | **Public long-term ECG with documented SCA** (e.g. selected THEW / hospital repositories when available) | Longer pre-event windows | Data-use agreements |

**Minimum independent target for primary claim:** ≥20 additional pre-VF Holters from sources 1–2 and/or 4, with ≥3 h continuous RR before the indexed event when possible (shorter windows allowed with explicit sensitivity strata).

---

## 3. Inclusion / exclusion criteria

### Inclusion
- Continuous ECG or beat annotations sufficient to reconstruct RR (ms).
- Documented spontaneous VF **or** sustained VT degenerating to VF with a clear onset time.
- Pre-event recording ≥ 60 min preferred; ≥ 15 min absolute minimum for short-DB strata.
- Adult human subjects.

### Exclusion
- Continuously paced rhythm without usable intrinsic RR.
- No usable beat annotations or >15% interpolated / invalid RR after standard cleaning.
- Induced VF only (EP lab induction) for primary analysis (may form a separate exploratory arm).
- Event time uncertain beyond ±2 min.

### Substrate strata (pre-specified)
- Sinus / mostly sinus  
- Atrial fibrillation  
- Intermittent pacing (flagged, not primary stratum)

---

## 4. Frozen analysis pipeline

1. Extract RR as in `code/extract_rr.py` (250–2000 ms, linear interp of short gaps).  
2. Bivariate proxy: \(X = [z(\mathrm{RR}),\, z(|\Delta\mathrm{RR}|)]\).  
3. τ_s via `systemictau.compute_taus` (W=101, stride=5).  
4. excess3 via nested ordinal RECD (m=3, delay=1, W=101, θ₃=0.08).  
5. Classic EWS comparators: rolling variance and lag-1 AR (W_EWS=501, stride=5).  
6. Windows: basal = early stable segment; approach = last 3 h before event (or last available third if duration < 6 h).  
7. Lead-time detector: absolute z-score vs basal ≥ 2.0 sustained for ≥ 3 consecutive metric samples (`cctp_metrics_core.detect_lead_time`).

**Primary endpoints (frozen)**
- Sign and magnitude of Δτ_s and Δexcess3 (approach − basal).  
- Lead-time distribution for τ_s and excess3.  
- Sensitivity at the frozen z-rule; specificity / FAR only on negative-control Holters.

**Secondary**
- Concordance of τ_s vs excess3.  
- Head-to-head vs var / AR(1) (same windows).  
- Surrogate phase-shuffle tests (n ≥ 8) on a random subset.

---

## 5. Negative controls and false-alarm rate

| Control type | Definition | Use |
|--------------|------------|-----|
| Matched non-event Holters | Same device/duration, no sustained VT/VF | FAR / specificity |
| Within-record remote windows | Early segment far from any arrhythmia | Secondary FAR estimate |
| AF-only non-VF Holters | AF without malignant ventricular events | Substrate-matched specificity |

Report FAR per 24 h of monitoring at the frozen z-rule. If controls are unavailable in a public phase, state sensitivity-only results and defer FAR to institutional data.

---

## 6. Success criteria (pre-specified)

| Criterion | Pass bar |
|-----------|----------|
| **S1 Relational signal** | Median \|Δτ_s\| or \|Δexcess3\| significantly > phase-shuffle null on ≥70% of independent events (direction free / context-dependent allowed) |
| **S2 Concordance** | Sign concordance τ_s–excess3 ≥ 0.70 on independent set |
| **S3 Lead time** | Among detections, median lead time ≥ 15 min (short DBs) or ≥ 30 min (Holter ≥ 6 h) |
| **S4 Detector** | Sensitivity ≥ 0.60 at frozen z-rule on independent pre-VF set |
| **S5 Specificity** | FAR ≤ 2 alarms / 24 h on negative-control Holters (when available) |
| **S6 Superiority (exploratory)** | Relational metrics not required to beat variance on every case; report side-by-side AUC/lead-time honestly |

Failing S1–S2 → revise theory claims before clinical tooling. Failing only S4–S5 → keep mechanistic paper, tighten detector or require multi-metric fusion.

---

## 7. Sample size (order-of-magnitude)

For sensitivity 0.60 with 95% CI width ±0.15 (Wilson), ~40 independent events are desirable; **phase-1 external** target is n≥20 with transparent CIs. SDDB extension alone (remaining processable records) is **supporting**, not sufficient for S4–S5.

---

## 8. Governance and reporting

- Register analysis code commit hash before unblinding institutional labels when possible.  
- Report all excluded records with reasons (mirror `data/sddb_full_inventory.csv` schema).  
- No FDA/clinical-deployment claim until multi-source S4–S5 are met.  
- Cite PhysioNet data DOIs and institutional approvals.

---

## 9. Timeline (suggested)

1. **Phase 0 (done in-repo):** SDDB inventory + N=10 pipeline + lead-time / EWS head-to-head + frozen code.  
2. **Phase 1 (done & reported 2026-07-12):** Public VFDB/CU + NSRDB controls; independent n=11; FAR ~29–34/24 h (S5 fail). See `docs/EXTERNAL_VALIDATION_PHASE1_REPORT.md`.  
3. **Phase 2 (2026-07-12):** **Specificity / realistic FAR** — plan + IRB checklist; **public interim done** (full NSRDB n=18, primary FAR τ_s ~33.7 / excess3 ~32.3 per 24h, device mismatch remains; S5 not met). Next: institutional / device-matched controls. Plan: `docs/EXTERNAL_VALIDATION_PHASE2_PLAN.md`; progress: `docs/EXTERNAL_VALIDATION_PHASE2_PROGRESS.md`.  
4. **Phase 3:** Pre-specified manuscript addendum / validation paper (after Phase 2 data).

**Phase 2 scientific goal (not more short-DB sensitivity):** estimate FAR on higher-quality controls; do not retune θ₃ / high-threshold / W_TAU / abs-z primary rule; no clinical claim until multi-source S4–S5 are met under adequate controls.

---

## 10. Related repo artifacts

- `data/sddb_full_inventory.csv` — 23-record inventory with exclusion reasons  
- `results/cctp_batch_summary.csv` / `cctp_cohort_stratified.csv` — discovery cohort  
- `results/leadtime_*` / `ews_head2head_*` — detector and comparator baselines  
- `results/external_phase1_*` — Phase 1 FAR / sensitivity (v1.1.0)  
- `results/phase2_public_control_inventory.csv` — Phase 2 candidate public/institutional control inventory  
- `docs/EXTERNAL_VALIDATION_PHASE2_PLAN.md` — Phase 2 specificity plan  
- `docs/PHASE2_IRB_DATA_CHECKLIST.md` — IRB / data-request checklist  
- `code/cctp_metrics_core.py` — pure lead-time / detector / FAR functions (unit-tested)
