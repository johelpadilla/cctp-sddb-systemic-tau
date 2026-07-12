# Phase 2 — IRB / Data-Request Checklist

**Purpose:** Prepare a **device-matched or institutional** non-event control cohort (and optional longer pre-VF series) for External Validation Phase 2.  
**Base release:** v1.1.0 · frozen primary detector only · **no clinical claim**.  
**Companion plan:** `docs/EXTERNAL_VALIDATION_PHASE2_PLAN.md`  
**Partner-facing request brief (sendable skeleton):** `docs/PHASE2_INSTITUTIONAL_DATA_REQUEST.md`

Use this checklist when talking to a clinical partner, IRB office, or data steward. Prefer **quality over quantity**.

**Progress note (2026-07-12):** Partner brief drafted; public NSRDB interim FAR complete (n=18, high FAR, device mismatch). Checklist §A–B items below remain **open until a real partner and ethics path are identified**. Do not invent institutional IDs or DUA signatures.

---

## A. Protocol & ethics

| # | Item | Status (Y/N/NA) | Notes |
|---|------|-----------------|-------|
| A1 | Written research protocol describing secondary use of Holter/telemetry ECG for relational EWS research | N | Partner brief + Phase 2 plan exist; formal protocol text still needed with partner |
| A2 | IRB / ethics committee review path identified (full / expedited / exempt determination) | N | Pending partner institution |
| A3 | Waiver of consent requested if retrospective only (justify minimal risk + impracticability) | N | Pending ethics path |
| A4 | Protocol states **no clinical decision support** and **no real-time alarm deployment** in this phase | Y (draft) | Explicit in partner brief §7 / §9 and this checklist §F |
| A5 | Analysis code commit hash will be recorded before unblinding labels when feasible | N | Process stated; no institutional labels yet |
| A6 | Data retention and destruction plan documented | N | Follow partner DUA when executed |

---

## B. De-identification & sharing

| # | Item | Status | Notes |
|---|------|--------|-------|
| B1 | Direct identifiers removed (name, MRN, address, phone, email) | Proposed | Scheme in partner brief §5; not yet applied to real export |
| B2 | Dates shifted or relative time-to-event only (per local policy) | Proposed | Relative hours / shifted dates |
| B3 | Free-text fields scrubbed of re-identifying content | Proposed | |
| B4 | Preferred export: beat annotations / RR series + event timestamps (not full raw ECG) if sufficient | Proposed | RR-first order in partner brief §5.3 |
| B5 | If raw ECG needed: stored outside git; access-controlled; **never** pushed to GitHub | Proposed | Zero institutional PHI in repo today |
| B6 | DUA / data-sharing agreement required? (yes → attach template) | N | Required when partner named; no DUA yet |
| B7 | Only aggregate tables and de-ID research IDs may enter the public repo | Y (policy) | Enforced in brief + inventory (Tier A templates only) |

---

## C. Control (non-event) inclusion — data steward request

Request Holters / long-term ECGs that meet:

| # | Criterion | Required |
|---|-----------|----------|
| C1 | Adult subjects | Yes |
| C2 | Continuous ECG or reliable beat annotations for RR reconstruction | Yes |
| C3 | **No** sustained VT/VF during the exported window | Yes (primary FAR) |
| C4 | Duration ≥ 4 h preferred; ≥ 24 h ideal for diurnal coverage | Preferred |
| C5 | Device model / sampling rate known (or explicitly “unknown”) | Preferred |
| C6 | Same device family as event cases when event arm exists (**device-matched**) | **Strongly preferred** |
| C7 | Substrate label if known (sinus / AF / paced / mixed) | Preferred |
| C8 | Interpolation or noise burden not expected to exceed ~15% invalid RR after cleaning | Preferred |
| C9 | Target **n**: start with quality n≈10–20 matched controls rather than hundreds of heterogeneous strips | Quality first |

**Explicit non-goals for the first institutional pull:** random short strips, EP-lab induced VF as “controls,” pediatric-only series without a separate protocol.

---

## D. Event arm (optional concurrent pull)

| # | Criterion | Required |
|---|-----------|----------|
| D1 | Spontaneous VF or sustained VT degenerating to VF | Yes for event arm |
| D2 | Onset time uncertainty ≤ ±2 min | Yes |
| D3 | Pre-event continuous recording ≥ 3 h preferred (≥ 15 min absolute floor only as short stratum) | Preferred |
| D4 | Induced EP-lab VF excluded from primary | Yes |
| D5 | Continuously paced without usable intrinsic RR excluded | Yes |

---

## E. Minimum metadata fields (de-ID research table)

Export one row per recording (CSV/JSON) with at least:

```
research_id, is_control, has_vf_or_sustained_vt, event_onset_rel_h,
total_duration_h, pre_event_h, device_model, sampling_hz,
substrate, annotation_source, export_date, notes
```

Do **not** include MRN, name, or absolute calendar dates if policy forbids.

---

## F. Analysis constraints (share with collaborator)

1. Primary detector is **frozen**: θ₃=0.08, high-threshold=0.65, W_TAU=101, abs-z ≥ 2 sustained ≥ 3 windows.  
2. **No retuning** of primary thresholds on the institutional set.  
3. FAR definition: alarm episodes / search hours × 24, same as Phase 1 public pilot.  
4. Exploratory fusion / basal variants only if pre-registered and labeled separately.  
5. Results remain research-only until multi-source S4–S5 are met.

---

## G. Suggested request email skeleton

> We request a de-identified research extract of long-term ECG / Holter recordings for a frozen-parameter early-warning methods study (no clinical decision support).  
> **Controls:** adult recordings without sustained VT/VF, preferably ≥4–24 h, same device family as any event cases, with beat annotations or RR series.  
> **Events (optional):** spontaneous VT/VF with clear onset and ≥3 h pre-event when available.  
> We will not retune detection thresholds on these data for the primary analysis and will not place PHI in public repositories.

---

## H. Public interim path (if no institutional partner yet)

| # | Action | Claim limit | Status |
|---|--------|-------------|--------|
| H1 | Expand PhysioNet NSRDB beyond the 6 Phase-1 Holters (full list = 18 subjects) | Still **not** device-matched; FAR remains a public healthy-Holter estimate | **Done** (n=18; τ_s ~33.7/24h; excess3 ~32.3/24h) |
| H2 | Inventory AF-without-VF public series as exploratory substrate controls | Separate exploratory arm only | Schema only (placeholder row) |
| H3 | Do **not** invent institutional rows | Inventory stays empty for Tier A until real access | **Enforced** (RESEARCH_ID_TBD only) |

See `results/phase2_public_control_inventory.csv`.  
**Active path now:** Tier A via `docs/PHASE2_INSTITUTIONAL_DATA_REQUEST.md` (do not further expand healthy Holter as main arm).
