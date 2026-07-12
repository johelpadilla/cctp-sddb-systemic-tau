# CCTP — Handoff (authoritative status)

**Date**: 2026-07-12 (v1.2.0 public close-out: manuscript final + Phase 2 interim)  
**Root**: `/Users/johelpadilla/grok-safe/Investigaciones/Cardiac_CCTP_Pilot/`  
**GitHub**: https://github.com/johelpadilla/cctp-sddb-systemic-tau  
**Release**: **v1.2.0** — manuscript final + Phase 2 public interim  
**Zenodo concept DOI**: https://doi.org/10.5281/zenodo.21270698  
**Prior version DOI (v1.1.0)**: https://doi.org/10.5281/zenodo.21326738  

**Status snapshot**
| Workstream | State |
|------------|--------|
| SDDB discovery N=10 (stratified + lead-time + H2H) | **Done** |
| Manuscript + PDF (Phase 1 §3.9 + Phase 2 §3.10) | **Done** (v1.2.0) |
| External Validation Phase 1 (VFDB + NSRDB FAR) | **Done & reported** |
| External Validation Phase 2 public interim (NSRDB n=18) | **Done & reported** |
| Tier A institutional data-request draft | **Prepared** (no partner data yet) |
| GitHub push + Zenodo archive | **v1.2.0** (auto-archive via GitHub↔Zenodo) |
| Clinical Copilot | Draft only — **do not deploy** |

**Frozen params** (never retune on validation without pre-registration):  
θ₃=0.08, high-threshold=0.65, W_τ=101, W_EWS=501, stride=5, relative λ, detector abs-z≥2 × 3 consecutive windows, RR clean [250, 2000] ms.

**Clinical / FDA / deployability claim**: **NONE**. Phase 1 and Phase 2 interim FAR fail S5 by a wide margin; **S5 is not claimed**.

---

## v1.2.0 public close-out

1. Polished manuscript Abstract / §3.10 / Discussion / Limitations / Conclusions for honest Phase 2 interim FAR without diluting the discovery relational message.
2. Regenerated preprint PDF (`manuscript/CCTP_SDBB_manuscript.pdf`).
3. Shipped Phase 2 artifacts, Tier A request brief, IRB checklist, planning docs, and structural tests.
4. Version metadata: `CITATION.cff`, `.zenodo.json`, `README.md` → 1.2.0 messaging.

### Primary FAR (frozen rule)

| Fact | Phase 1 | Phase 2 interim |
|------|---------|-----------------|
| NSRDB n | 6 | **18** |
| Search hours (capped 12 h/record) | ~60 | **~180** |
| τ_s FAR / 24 h | ~34.4 | **~33.7** |
| excess3 FAR / 24 h | ~28.8 | **~32.3** |
| Device-matched to VFDB? | **No** | **No** |
| S5 (FAR ≤2/24h) | Failed | **Still failed / not claimed** |
| Institutional device-matched FAR | — | **Not yet** (request drafted) |
| Clinical claim | false | false |

---

## Source-of-truth documents

| Artifact | Role |
|----------|------|
| `manuscript/CCTP_SDBB_manuscript.md` + `.pdf` | Final public manuscript package |
| `docs/PHASE2_INSTITUTIONAL_DATA_REQUEST.md` | Partner-facing Tier A data request |
| `docs/PHASE2_IRB_DATA_CHECKLIST.md` | IRB / partner checklist |
| `docs/EXTERNAL_VALIDATION_PHASE2_PLAN.md` | Phase 2 plan |
| `docs/EXTERNAL_VALIDATION_PHASE2_PROGRESS.md` | Progress + recommendation |
| `results/external_phase2_far.json` | Full NSRDB primary FAR |
| `results/external_phase2_summary.json` | Phase 2 summary (claims flags) |
| `results/phase2_public_control_inventory.csv` | Control inventory (Tier A empty of PHI) |
| `docs/EXTERNAL_VALIDATION_PHASE1_REPORT.md` | Phase 1 numeric baseline |

### Reproduce Phase 1 + Phase 2 interim FAR
```bash
cd /Users/johelpadilla/grok-safe/Investigaciones/Cardiac_CCTP_Pilot
python3 code/run_external_validation_phase1.py
python3 code/run_external_validation_phase2_far.py
python3 -m pytest tests/test_leadtime_detector.py tests/test_far_and_short_windows.py \
  tests/test_phase2_planning_artifacts.py tests/test_phase2_far_artifacts.py -q
```

---

## Discovery N=10 (foundational)

| Finding | Value |
|---------|--------|
| event_type | intermediate 6 (30,32,38,45,47,50); terminal 4 (31,35,36,51) |
| τ_s–excess3 concordance | **8/10 (0.8)** |
| Median lead excess3 / τ_s / var (h) | **6.86 / 5.88 / 3.90** |

Artifacts: `results/cctp_*stratified*`, `leadtime_*`, `ews_head2head*`, `figures/publication/`.  
Interpretation: `docs/JUL12_RESULTS_INTERPRETATION.md`.

---

## Next scientific step (human / institutional)

Public healthy-Holter expansion is **closed** as the main specificity arm.

1. Identify **one** real partner in classes P1–P5 (`docs/PHASE2_INSTITUTIONAL_DATA_REQUEST.md`).
2. Complete checklist §A–B with their IRB/DUA.
3. Export de-ID RR/annotations for n≈10–20 quality-first device-matched non-event controls.
4. **Do not** invent institutional inventory rows, retune frozen params, or claim S5/clinical deployability.

### Next-session prompt (copy-paste)

```
Continue External Validation Phase 2 — Tier A institutional path (post v1.2.0).
Base: public manuscript final + NSRDB FAR n=18 (τ_s ~33.7/24h, excess3 ~32.3/24h, device mismatch).
Partner brief: docs/PHASE2_INSTITUTIONAL_DATA_REQUEST.md
Checklist: docs/PHASE2_IRB_DATA_CHECKLIST.md
FROZEN: θ₃=0.08, high=0.65, W_TAU=101, abs-z≥2×3 — do not retune.
Priority: identify one real P1–P5 partner; complete checklist §A–B; no invented PHI; no S5/clinical claim.
Root: /Users/johelpadilla/grok-safe/Investigaciones/Cardiac_CCTP_Pilot
```

---

## Claims boundary

| Claim | Status |
|-------|--------|
| Discovery relational pre-VF reorganization (SDDB N=10) | **Reported** (primary message) |
| External sensitivity (Phase 1 VFDB) | **Reported** |
| Public interim control FAR (NSRDB n=18) | **Reported** with mismatch caveat |
| S5 met | **False / not claimed** |
| Clinical deployability / FDA | **None** |
| Device-matched institutional FAR | **Not yet available** |
