# Phase 2 — Institutional / Device-Matched Data Request (Partner Brief)

**Document type:** Partner-facing research data request (draft skeleton)  
**Study:** CCTP / SDDB relational early-warning methods — External Validation Phase 2  
**Base release:** v1.1.0 · commit `0a4cdbc` · DOI [10.5281/zenodo.21326738](https://doi.org/10.5281/zenodo.21326738)  
**Status:** Draft ready for partner discussion · **IRB/DUA not completed** · **no institutional data on hand**  
**Companions:** `docs/PHASE2_IRB_DATA_CHECKLIST.md` · `docs/EXTERNAL_VALIDATION_PHASE2_PLAN.md` · `docs/EXTERNAL_VALIDATION_PHASE2_PROGRESS.md`

**Honesty boundary (read first)**

| Statement | Status |
|-----------|--------|
| Clinical decision support / real-time alarms | **Not requested and not claimed** |
| Deployability / FDA / bedside use | **None** |
| S5 (FAR ≤ 2 / 24 h) | **Not met** on public healthy Holters; **not claimed** |
| Primary detector thresholds | **Frozen** — will **not** be retuned on institutional data |
| Institutional PHI | **Zero** in public git; none requested until ethics + DUA |

This brief is a **request skeleton**, not a signed agreement and not evidence that any partner has exported data.

---

## 1. Why we are asking (scientific motivation)

Phase 1 external validation (public PhysioNet) showed:

- Independent short VT/VF episodes (VFDB): frozen detector is **sensitive** on usable pre-event windows.
- Negative controls (NSRDB healthy Holter, expanded Phase 2 interim **n = 18**): primary FAR remains **high** (~33.7 / 24 h for τ_s; ~32.3 / 24 h for excess3 under a 12 h search cap).
- **Device mismatch:** healthy ambulatory Holter ≠ clinical telemetry / device-matched non-event segments. Public expansion did **not** fix specificity estimation quality.

**Phase 2 scientific goal:** estimate **specificity / false-alarm rate (FAR)** under **higher-quality negative controls** — preferably **device-matched** and from a **clinical non-event** population — using the **same frozen primary rule**. Quality of controls matters more than adding more healthy Holters.

---

## 2. Candidate partner / archive classes (Tier A)

We are **not** inventing named hospital sites or research IDs here. Suitable partner classes include any of the following that can provide de-identified continuous ECG under local ethics/DUA:

| Class ID | Partner / archive class | Why useful for Tier A | Typical access path |
|----------|-------------------------|----------------------|---------------------|
| **P1** | Hospital or health-system **Holter lab archive** (outpatient / ambulatory) | Long duration (often ~24 h); strong for FAR denominators | IRB/ethics + DUA; data steward export |
| **P2** | Inpatient **telemetry / monitored-bed ECG archive** | Closer to malignant-arrhythmia clinical context; device often closer to VFDB-like telemetry | Same + clinical engineering / informatics |
| **P3** | **Device-vendor or core-lab** Holter/telemetry repository used by a hospital (brand/firmware documented) | Best path to **device-matched** controls when event cases exist on the same platform | DUA with institution and/or vendor terms |
| **P4** | Electrophysiology / arrhythmia service **clinical research archive** (retrospective secondary use) | May hold both non-event long records and spontaneous VT/VF with timed onset | IRB secondary-use protocol |
| **P5** | Multi-center **research ECG repository** already under governance (if local site can sponsor secondary analysis) | Scale + metadata standards | Existing governance + amendment |

**Preferred match order:** (P3 or same-device P1/P2) → any clinical non-event long ECG with known device → last resort heterogeneous clinical Holters (still better than healthy-only public NSRDB if device is documented).

**Out of scope for first pull as primary FAR controls:** random short ECG strips; EP-lab **induced** VF labeled as “controls”; pediatric-only series without a separate protocol; any export that cannot remove direct identifiers.

---

## 3. What we request — non-event controls (primary FAR arm)

### 3.1 Inclusion criteria (must / preferred)

Aligned with checklist §C (`docs/PHASE2_IRB_DATA_CHECKLIST.md`):

| # | Criterion | Level |
|---|-----------|--------|
| C1 | Adult subjects | **Required** |
| C2 | Continuous ECG **or** reliable beat annotations sufficient to reconstruct RR | **Required** |
| C3 | **No** sustained VT/VF during the exported analysis window | **Required** (primary FAR) |
| C4 | Duration **≥ 4 h preferred**; **≥ 24 h ideal** (diurnal coverage) | Preferred |
| C5 | Device model and sampling rate known (or explicitly “unknown”) | Preferred |
| C6 | **Same device family** as any concurrent event cases when an event arm exists (**device-matched**) | **Strongly preferred** |
| C7 | Substrate label if known (sinus / AF / paced / mixed) | Preferred |
| C8 | Noise / invalid RR after cleaning not expected to exceed ~15% | Preferred |
| C9 | Target sample: **n ≈ 10–20 good-quality** matched controls first | **Quality over quantity** |

### 3.2 Explicit exclusions (first institutional pull)

- Sustained VT/VF anywhere in the exported control window  
- Induced EP-lab VF  
- Continuously paced without usable intrinsic RR (unless a separate paced exploratory stratum is agreed in writing)  
- Pediatric-only without separate protocol  
- Strips too short to support a basal window (~2 h early) + search period

### 3.3 Initial sample-size target (quality-first)

| Target | Rationale |
|--------|-----------|
| **n ≈ 10–20** non-event recordings of good quality | Enough search hours to re-estimate FAR with a **clinical / device-matched** caveat; avoids a large heterogeneous dump |
| Prefer **full Holter (≥12–24 h)** over many 1 h fragments | Stable FAR denominator; diurnal structure |
| Expand only after first quality cohort is analyzed under frozen rule | Do not trade match quality for n |

We will **not** treat public NSRDB n=18 as a substitute for this cohort.

---

## 4. Optional concurrent pull — long pre-VF / event arm

If the same archive can supply events **without delaying** the control pull:

| # | Criterion | Level |
|---|-----------|--------|
| D1 | Spontaneous VF or sustained VT degenerating to VF | Required for event arm |
| D2 | Onset time uncertainty ≤ ±2 min | Required |
| D3 | Pre-event continuous recording **≥ 3 h preferred** (≥ 15 min only as short stratum, labeled) | Preferred |
| D4 | Induced EP-lab VF excluded from primary | Required |
| D5 | Continuously paced without usable intrinsic RR excluded | Required |

Event arm supports lead-time strata closer to discovery SDDB geometry; it is **not** required to start the control FAR path.

---

## 5. Proposed de-identification scheme

Aligned with checklist §B. Partner applies local policy; we request the **minimum** needed for RR-based analysis.

### 5.1 Remove / never transfer

- Name, MRN, address, phone, email, insurance IDs  
- Free-text fields that may re-identify (scrub or drop)  
- Absolute calendar dates if policy requires (use relative time instead)  
- Full raw waveform **unless** RR/annotations are insufficient (see §5.3)

### 5.2 Preferred research identifiers and time

| Field practice | Recommendation |
|----------------|----------------|
| Subject / recording ID | Opaque **research_id** only (no MRN fragment) |
| Time base | Relative seconds/hours from recording start; for events, **time-to-onset** relative hours |
| Dates | Shifted or omitted per IRB/DUA |
| Site codes | Optional non-identifying site code if multi-site |

### 5.3 Preferred data products (in order)

1. **Best:** beat annotations or cleaned **RR interval series** + event timestamps (if any) + metadata table  
2. **Acceptable:** de-identified annotation files sufficient to rebuild RR offline  
3. **Only if needed:** raw ECG under access control, **outside public git**, never pushed to GitHub  

### 5.4 What may enter the public research repository later

- Aggregate tables (FAR, hours, n, inclusion counts)  
- De-identified research IDs **if** DUA allows  
- **Never:** raw institutional waveforms, MRNs, names, or unshifted clinical dates  

### 5.5 Retention

Data retention and destruction follow the partner DUA and local IRB; analysis code commit hash should be recorded before unblinding labels when feasible (checklist A5–A6).

---

## 6. Minimum metadata export (one row per recording)

```
research_id, is_control, has_vf_or_sustained_vt, event_onset_rel_h,
total_duration_h, pre_event_h, device_model, sampling_hz,
substrate, annotation_source, export_date, notes
```

Do **not** include MRN, name, or absolute calendar dates if policy forbids.  
Schema matches `docs/PHASE2_IRB_DATA_CHECKLIST.md` §E.

---

## 7. Analysis constraints (frozen primary rule — non-negotiable for primary)

Share this section with the collaborator as written:

1. **Primary detector is frozen:**
   - θ₃ = **0.08**
   - high-threshold = **0.65**
   - W_TAU = **101**
   - Detector: **abs-z ≥ 2** sustained **≥ 3** consecutive windows  
2. **No retuning** of primary thresholds on the institutional set for the primary analysis.  
3. FAR definition: alarm **episodes** / search hours × 24 (same convention as Phase 1 / Phase 2 public interim).  
4. Exploratory fusion / basal variants only if **pre-registered** and labeled separately — never silent replacement of the primary rule.  
5. Results remain **research-only**. Meeting multi-source S4–S5 under adequate controls would still not authorize clinical deployment in this phase.  
6. **No clinical decision support** and **no real-time alarm deployment** are part of this request.

Public interim reference (for context only, not a clinical FAR): full NSRDB n=18, τ_s FAR ~33.7/24 h, excess3 ~32.3/24 h, **device_mismatch = true**, **s5_claim = false**.

---

## 8. Ethics / governance checklist pointer

Before any export:

| Gate | Document |
|------|----------|
| Protocol + no-CDS language | Checklist §A (`PHASE2_IRB_DATA_CHECKLIST.md`) |
| De-ID + DUA | Checklist §B |
| Inclusion / exclusion | This brief §§3–4 + checklist §C–D |
| Metadata | §6 + checklist §E |
| Analysis freeze | §7 + checklist §F |

**Status of this draft:** content prepared for partner discussion. Fields A1–A6 / B1–B7 remain **not completed** until a real partner and ethics path are identified. No fabricated IRB numbers or DUA signatures.

---

## 9. Suggested outreach email (copy-paste skeleton)

**Subject:** Research data request — de-identified long-term ECG for frozen-parameter early-warning methods study (no clinical decision support)

> Dear [Data steward / collaborator],  
>  
> We request a **de-identified** research extract of long-term ECG / Holter (or telemetry) recordings for a **methods** study of relational early-warning signals prior to malignant ventricular arrhythmia. This phase is **research only**: **no clinical decision support**, **no real-time alarms**, and **no claim of clinical deployability**.  
>  
> **Primary need — non-event controls (quality first):**  
> - Adult recordings **without** sustained VT/VF in the export window  
> - Preferably **≥ 4–24 h** continuous  
> - Beat annotations or RR series preferred over raw ECG when sufficient  
> - **Same device family** as any event cases when available (**device-matched**)  
> - Initial target **n ≈ 10–20** good-quality controls (not a large heterogeneous dump)  
>  
> **Optional concurrent need — events:** spontaneous VT/VF with clear onset and **≥ 3 h** pre-event when available (induced EP-lab VF excluded from primary).  
>  
> **Analysis constraints:** detection thresholds are **frozen** (θ₃=0.08, high-threshold=0.65, W_TAU=101, abs-z ≥ 2 sustained ≥ 3) and will **not** be retuned on these data for the primary analysis. PHI will not be placed in public repositories.  
>  
> We can share our full checklist (`PHASE2_IRB_DATA_CHECKLIST.md`) and this request brief, plus public interim results (PhysioNet healthy Holter FAR under device-mismatch caveat).  
>  
> Thank you for considering this secondary-use research extract.  
>  
> Sincerely,  
> [Investigator]

---

## 10. What success looks like for this request (not S5 claims)

| Milestone | Honest success definition |
|-----------|---------------------------|
| Partner engaged | Named class P1–P5 path identified; ethics route known |
| Export received | De-ID RR/annotations + metadata for n≈10–20 controls; device documented when possible |
| Analysis | Primary frozen FAR reported with device-match flag; `clinical_claim: false` |
| S5 | Only discuss if measured under **adequate** matched controls; **do not claim** from public NSRDB |

---

## 11. Next human steps (after this draft)

1. Identify **one concrete partner** in classes P1–P5 (hospital Holter lab, telemetry archive, or research EP service) — **not** inventing names in-repo until real contact exists.  
2. Complete checklist §A–B with that partner’s IRB/DUA language.  
3. Agree export format (RR/annotations first) and research_id scheme.  
4. Only then fill inventory Tier A rows with real (de-ID) research IDs — never fabricate.  
5. Run frozen primary FAR; keep public NSRDB as interim reference only.

---

## 12. Document control

| Item | Value |
|------|--------|
| Path | `docs/PHASE2_INSTITUTIONAL_DATA_REQUEST.md` |
| Derived from | `PHASE2_IRB_DATA_CHECKLIST.md`, `EXTERNAL_VALIDATION_PHASE2_PLAN.md` Tier A, Phase 2 interim progress |
| Institutional data present | **No** |
| Primary rule | Frozen as §7 |
| Clinical claim | **false** |
