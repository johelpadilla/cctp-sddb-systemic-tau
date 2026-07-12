# External Validation Phase 2 — Specificity / Realistic FAR

**Status:** planning started 2026-07-12 (post v1.1.0)  
**Base release:** GitHub `v1.1.0` (commit `0a4cdbc`) · Zenodo DOI [10.5281/zenodo.21326738](https://doi.org/10.5281/zenodo.21326738)  
**Parent plan:** `docs/EXTERNAL_VALIDATION_PLAN.md`  
**Phase 1 report:** `docs/EXTERNAL_VALIDATION_PHASE1_REPORT.md`  
**Clinical / FDA / deployability claim:** **NONE**

---

## 0. Scientific goal (what Phase 2 is — and is not)

| Phase 2 **is** | Phase 2 is **not** |
|----------------|--------------------|
| Estimating **specificity / false-alarm rate (FAR)** on **higher-quality negative controls** | More short-DB sensitivity hunting (VFDB/CU) as the main deliverable |
| Preferring **device-matched** or **institutional** non-event Holters over larger public n | Claiming S5 (FAR ≤ 2 / 24 h) without adequate controls |
| Keeping **primary detector frozen** from discovery / Phase 1 | Retuning θ₃, high-threshold, W_TAU, or abs-z rule on validation or control data |
| Documenting IRB / de-ID / inclusion paths for real clinical archives | Shipping Copilot alarms or clinical tooling |

**Primary scientific bottleneck (from Phase 1):** the frozen abs-z rule is **sensitive but not specific**. Phase 2 exists to measure FAR under controls that better match the clinical acquisition setting — not to chase another public short-episode sensitivity table.

---

## 1. Why Phase 1 forces this focus (honest numbers)

All figures below are from `results/external_phase1_summary.json` / Phase 1 report. **No clinical claim.**

| Fact | Value | Implication for Phase 2 |
|------|-------|-------------------------|
| Independent processable events | **VFDB n = 11** (all `short_15_60min`) | Lead-time strata are minutes, not multi-hour Holters |
| τ_s sensitivity | **1.00** (11/11) | Sensitivity is **not** the binding constraint |
| excess3 sensitivity | **≈ 0.82** (9/11) | Same |
| Median lead τ_s / excess3 | **~10.4 / ~7.6 min** | Short DBs cap lead time; not comparable to SDDB multi-hour leads |
| CUDB processable events | **0** (all pre-event &lt; 15 min) | Public short VF corpora are nearly exhausted under floor |
| NSRDB controls | **n = 6**, ~60 search hours | FAR is **quantified**, not nan |
| FAR τ_s / excess3 | **~34.4 / ~28.8 per 24 h** | **S5 (≤ 2 / 24 h) failed hard** |
| Fraction of controls with ≥1 alarm | **1.0** | Rule floods healthy Holters |
| Control match quality | NSRDB = rhythm-healthy Holter; **not device-matched** to VFDB telemetry | FAR may over- or under-state clinical rate |
| S1–S6 package | **Not claimed** as pass | Paper can report pilot; cannot claim deployability |

**Bottom line:** Phase 1 closed the “FAR = nan” gap with an honest high FAR. The next scientific step is **better controls and (only if pre-registered) exploratory specificity machinery** — not more VFDB sensitivity rows.

---

## 2. Frozen primary parameters (sole primary rule)

These are the **only** primary analysis rule for Phase 2 endpoints. They must appear in every Phase-2 artifact. **Do not retune** on any validation or control set.

| Parameter | Frozen value | Source constant |
|-----------|--------------|-----------------|
| θ₃ | **0.08** | `FROZEN_THETA3` |
| high-threshold (RECD reporting) | **0.65** | `FROZEN_HIGH_THRESHOLD` |
| W_TAU | **101** | `W_TAU` |
| stride | **5** | pipeline default |
| W_EWS (comparators) | **501** | manuscript convention |
| Detector | abs-z ≥ **2.0** sustained ≥ **3** consecutive windows | `FROZEN_Z_THRESHOLD`, `FROZEN_MIN_CONSECUTIVE` |
| RR clean | **[250, 2000] ms** | extractors |
| Episode refractory (FAR counting) | **0.5 h** (reporting convention from Phase 1; not a “threshold retune”) | `count_alarm_episodes` default |

Any change to the primary rule requires a **pre-registered protocol amendment** and a new labeled analysis stream (not silent replacement of the frozen primary).

---

## 3. Priority target sources (quality over quantity)

Order is intentional: **matched / institutional first**, public interim only as scaffolding.

### Tier A — highest scientific value (primary Phase 2 path)

| Priority | Source | Role | Why preferred | Access |
|----------|--------|------|---------------|--------|
| **A1** | **Institutional Holter / telemetry archives** with timed spontaneous VT/VF **and** matched non-event recordings | Primary event + control arms | Device, population, and annotation practice match clinical use | IRB + de-identification; DUA |
| **A2** | **Device-matched non-event controls** (same Holter brand/firmware, same sampling, same RR pipeline as event cases) | FAR / specificity | Removes Phase-1 NSRDB–VFDB device mismatch | Usually institutional |
| **A3** | **Longer pre-VF series** (≥3 h continuous RR preferred; ≥6–24 h ideal) | Lead-time strata comparable to SDDB discovery | Phase-1 leads are truncated by short public episodes | Institutional or selected public long SCA series under DUA |

### Tier B — public interim controls (useful, limited claim)

| Priority | Source | Role | Quality notes | Status in-repo |
|----------|--------|------|---------------|----------------|
| **B1** | **PhysioNet NSRDB** remaining subjects (full list has **18**; Phase 1 used **6**) | Expand search hours under same healthy-Holter caveat | Still **not** device-matched to VFDB; better n, same mismatch class | 6 local under `data/nsrdb/`; 12 not downloaded |
| **B2** | **Within-record remote windows** on long institutional (or long public) series far from any malignant event | Secondary FAR | Requires careful event-free window definition | Not yet operationalized for Phase 2 |
| **B3** | **AF-without-VF Holters** (e.g. selected Long-Term AF Database segments) | Substrate-matched **specificity** exploratory arm | Not “healthy NSR”; may elevate or change FAR | Not downloaded; inventory only |
| **B4** | THEW / other long-term ECG repositories with documented SCA | Long pre-event + controls when DUA allows | Access friction; quality varies | Planned / DUA only |

### Tier C — explicitly deprioritized for Phase 2 primary goal

| Source | Why not primary |
|--------|-----------------|
| More VFDB short episodes | Sensitivity already high; pre-event floor already binding |
| CUDB expansion | Phase 1: **0** processable under ≥15 min pre-event |
| Retuning thresholds on NSRDB to “pass S5” | **Forbidden** for primary claim |
| MIT-BIH Arrhythmia DB mixed strips | Short, heterogeneous annotations; poor FAR denominator |

**Inventory file:** `results/phase2_public_control_inventory.csv` (candidate rows + inclusion flags; no institutional PHI).

---

## 4. Specificity strategies (primary frozen vs exploratory)

### 4.1 Primary reporting stream (mandatory)

- Same abs-z ≥ 2 × 3 rule as Phase 1 / discovery.
- FAR = alarm **episodes** / search hours × 24 (same as `false_alarm_rate`).
- Report FAR separately for τ_s and excess3.
- Label control quality: device match (yes/no), population (healthy / inpatient / AF / etc.), hours searched.
- **Success bar S5 remains** FAR ≤ 2 / 24 h on **adequate** negative controls — Phase 2 may still fail S5 honestly.

### 4.2 Reporting conventions (not retunes)

These may be varied for **sensitivity tables** only if pre-declared; primary row keeps Phase-1 defaults:

| Convention | Phase-1 default | Notes |
|------------|-----------------|-------|
| Basal window length | Early ~2 h (control) | Longer / multi-segment basal is **exploratory** |
| Search cap per record | 12 h | Can extend when full Holters available |
| Refractory between episodes | 0.5 h | Changing refractory changes episode count; report both if varied |
| Max hours used | 12 h Phase 1 | Institutional 24 h+ preferred |

### 4.3 Exploratory / pre-registered only (must not replace primary)

Any of the following **must** be labeled `exploratory` or `pre_registered_amendment` and **must not** overwrite the frozen primary claim:

1. **Better basal referencing** — longer basal, rolling basal, multi-day circadian-aware baseline, excluding sleep/wake transitions from basal stats.
2. **Multi-metric fusion** — e.g. require τ_s **and** excess3 concurrent alarms; AND/OR with classical var only as comparator.
3. **Substrate-stratified FAR** — sinus vs AF vs paced controls reported separately.
4. **Direction-aware z** — signed alarms only (weaker default for relational polarity).
5. **Longer min_consecutive or higher z** — **threshold changes are not primary**; if studied, treat as a new detector version with a new name and pre-registration, never as “the” frozen rule.

**Rule of thumb:** if it would make S5 easier by changing the detector, it is either (a) a new pre-registered detector or (b) exploratory — never silent primary.

---

## 5. Data requirements

### 5.1 Negative-control (non-event) definition

A control recording is eligible for **primary FAR** if **all** hold:

1. Continuous ECG or beat annotations sufficient to reconstruct RR (ms) after the same cleaning as events.
2. **No** sustained VT/VF (and preferably no sustained VT) during the searched window; document arrhythmia burden if available.
3. Usable duration ≥ **4 h** preferred for primary FAR; ≥ **1 h** absolute minimum only for pilot inventory (flag as `duration_marginal`).
4. Adult human subjects.
5. Interpolation / invalid RR after cleaning ≤ **15%** (same as Phase 1 event floor spirit).
6. Device / sampling metadata recorded when known (even if “unknown”).

**Matched control (preferred):** same device family, similar recording duration, same RR extraction path as the event arm when events exist in the same archive.

### 5.2 Event-arm (longer pre-VF) definition — secondary Phase 2 aim

- Spontaneous VF or sustained VT → VF with onset time uncertainty ≤ ±2 min.
- Pre-event continuous RR ≥ **3 h** preferred (stratum `holter_ge3h`); ≥ **15 min** remains absolute floor for short strata only.
- Same RR clean and frozen metrics as Phase 1.
- Induced EP-lab VF excluded from primary event arm.

### 5.3 IRB / ethics / de-identification

| Requirement | Detail |
|-------------|--------|
| IRB / ethics approval | Required before any institutional identifiable data leave the clinical environment |
| De-identification | Remove names, MRNs, exact dates if policy requires; keep relative time-to-event and duration |
| Data use agreement | When sharing outside the approved protocol or to collaborators |
| Minimal necessary data | Prefer beat annotations + RR timestamps over full raw ECG if sufficient |
| Provenance log | Source system, device model if known, extraction date, commit hash of analysis code |
| No PHI in git | Raw institutional waveforms **never** committed; only de-ID aggregates / synthetic examples |

Checklist: `docs/PHASE2_IRB_DATA_CHECKLIST.md`.

### 5.4 Minimum analysis protocol for a Phase 2 FAR report

1. Freeze code commit hash **before** unblinding institutional labels when possible.  
2. Build inventory CSV with include/exclude reasons (schema in `results/phase2_public_control_inventory.csv` header + institutional extension columns).  
3. Run **primary** frozen detector FAR.  
4. Optionally run **pre-registered exploratory** streams in separate JSON keys.  
5. Report hours searched, n controls, FAR/24 h, fraction alarmed, device-match flag.  
6. Explicitly set `clinical_claim: false` until S4–S5 are met under adequate controls.

---

## 6. Phase 2 deliverables (staged)

| Stage | Deliverable | This session |
|-------|-------------|--------------|
| **P2.0** | Planning docs + frozen-param echo + Phase-1 bottleneck statement | **Yes** |
| **P2.1** | Public interim control inventory (quality flags; optional NSRDB expansion plan) | **Yes (inventory)** |
| **P2.2** | IRB / data-request checklist for institutional partners | **Yes** |
| **P2.3** | Optional: download remaining NSRDB subjects and re-run **primary** FAR only (same mismatch caveat) | Next session if chosen |
| **P2.4** | Institutional control cohort schema filled (no PHI) + first de-ID RR batch | Requires real partner access |
| **P2.5** | Phase 2 results report + manuscript addendum | Later; **not** this session |

---

## 7. Success criteria mapping (unchanged bars; honest status)

| Criterion | Pass bar | Phase 1 status | Phase 2 intent |
|-----------|----------|----------------|----------------|
| S1 Relational signal | Surrogate support ≥70% independent events | Not fully tested externally | Secondary until n and duration improve |
| S2 Concordance | ≥0.70 τ_s–excess3 | 0.64 on short VFDB | Revisit on longer pre-VF |
| S3 Lead time | ≥15 min short / ≥30 min Holter | Short-DB medians ~8–10 min | Needs longer pre-event series |
| S4 Detector sensitivity | ≥0.60 | Met on n=11 short set only | Maintain under frozen rule; do not hunt |
| **S5 Specificity** | **FAR ≤ 2 / 24 h** | **Failed (~29–34)** | **Primary Phase 2 target** |
| S6 Superiority | Exploratory H2H | Exploratory only | Keep exploratory |

---

## 8. Explicit non-claims

- No clinical decision support, FDA, or Copilot **deployability**.  
- No claim that expanding NSRDB alone “solves” specificity.  
- No claim that fusion “will” pass S5 until measured under primary + exploratory labels.  
- Public short databases are **not** the path to multi-hour lead-time external claims.

---

## 9. Related artifacts

| Path | Role |
|------|------|
| `docs/EXTERNAL_VALIDATION_PLAN.md` | Master validation plan (Phases 0–3) |
| `docs/EXTERNAL_VALIDATION_PHASE1_REPORT.md` | Phase 1 numeric source of truth |
| `docs/PHASE2_IRB_DATA_CHECKLIST.md` | IRB / data request checklist |
| `results/phase2_public_control_inventory.csv` | Candidate public non-event inventory |
| `results/external_phase1_*.json/csv` | Phase 1 FAR / sensitivity |
| `code/cctp_metrics_core.py` | Frozen constants + FAR helpers |
| `code/run_external_validation_phase1.py` | Reusable Phase-1 entry (extend carefully for Phase 2) |
| `HANDOFF.md` | Session status + next prompt |

---

## 10. Recommended immediate next actions (ordered)

1. **Partner path (highest impact):** complete IRB checklist items with a clinical collaborator; define device-matched non-event pull criteria (Tier A).  
2. **Public interim (optional, low claim):** download remaining NSRDB records (12) and re-estimate primary FAR under frozen rule — still report device-mismatch caveat; do **not** retune.  
3. **Schema only until data exist:** keep institutional inventory empty of PHI; fill metadata rows as partners agree.  
4. **Do not** start Copilot alarm shipping or manuscript “validation success” language.

### Next-session prompt (copy-paste)

```
Continue CCTP External Validation Phase 2 (specificity / FAR).
Base: v1.1.0 (Zenodo 10.5281/zenodo.21326738). Plan:
docs/EXTERNAL_VALIDATION_PHASE2_PLAN.md

Done this cycle: Phase 2 plan + IRB checklist + public control inventory.
Frozen primary only: θ3=0.08, high=0.65, W_TAU=101, abs-z≥2×3. No retune.

Concrete next (pick one, quality > quantity):
(A) Download remaining PhysioNet NSRDB subjects and re-run PRIMARY FAR only
    (same NSRDB mismatch caveat; still not device-matched), OR
(B) Fill institutional control request using docs/PHASE2_IRB_DATA_CHECKLIST.md
    (preferred if a clinical partner is available).

Do not claim clinical deployability; do not retune thresholds.
Root: /Users/johelpadilla/grok-safe/Investigaciones/Cardiac_CCTP_Pilot
```
