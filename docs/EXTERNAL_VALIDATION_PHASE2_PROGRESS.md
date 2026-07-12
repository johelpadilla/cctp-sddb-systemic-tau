# External Validation Phase 2 — Progress Report

**Date:** 2026-07-12  
**Base release:** v1.1.0 (`0a4cdbc`, DOI [10.5281/zenodo.21326738](https://doi.org/10.5281/zenodo.21326738))  
**Status:** Public interim arm **done** (full NSRDB n=18). **Tier A institutional data-request draft prepared.** **No S5 claim. No clinical claim.**  
**Primary rule:** **frozen** — θ₃=0.08, high-threshold=0.65, W_TAU=101, abs-z ≥ 2 sustained ≥ 3. **No retune.**

---

## 1. What was done

### 1.1 Public interim (completed earlier)

| Step | Result |
|------|--------|
| Download remaining NSRDB | **12/12** new records; full set **18/18** local |
| RR extract | **18/18** clean npz under `data/rr_external/` |
| Frozen primary FAR | `code/run_external_validation_phase2_far.py` (Phase-1 helpers; no retune) |
| Inventory | All NSRDB rows `present_data/nsrdb` |
| Artifacts | `results/external_phase2_far.json`, controls CSVs, summary JSON |

### 1.2 Tier A institutional path (this cycle)

| Step | Result |
|------|--------|
| Partner-facing request brief | **`docs/PHASE2_INSTITUTIONAL_DATA_REQUEST.md`** (sendable skeleton) |
| IRB checklist alignment | `docs/PHASE2_IRB_DATA_CHECKLIST.md` updated (honest N / Proposed statuses; H1 done) |
| Partner classes | Generic Tier A classes P1–P5 only — **no fabricated hospital names or PHI** |
| Inclusion / de-ID / n | Device-matched preferred; ≥4–24 h; n≈10–20 quality-first; de-ID scheme §5 |
| Frozen rule + no-claim language | Explicit in brief §§7–10 and checklist §F |
| Inventory Tier A | Templates only (`RESEARCH_ID_TBD`, `not_available`) — **no invented institutional rows** |

---

## 2. Primary FAR (frozen rule) — Phase 1 vs Phase 2 public interim

Control search window: early basal (~2 h) then search remainder, **capped at 12 h** per Holter (same as Phase 1).

| Metric | Phase 1 (n=6) | Phase 2 full NSRDB (n=18) |
|--------|---------------|---------------------------|
| Search hours | ~60.0 | ~180.0 |
| τ_s total episodes | 86 | 253 |
| τ_s **FAR / 24 h** | **~34.40** | **~33.73** |
| excess3 total episodes | 72 | 242 |
| excess3 **FAR / 24 h** | **~28.80** | **~32.27** |
| Fraction of controls with ≥1 alarm | 1.0 | 1.0 |

**Interpretation (honest):** Expanding from 6 → 18 healthy Holters **does not materially lower** primary FAR. The Phase 1 high-FAR finding was **not** an n=6 sampling fluke.

---

## 3. Limitations that remain

1. **Device mismatch (public arm):** NSRDB ≠ VFDB-like clinical telemetry; healthy Holter only.  
2. **No institutional FAR yet:** Tier A request is drafted; **zero** institutional RR/PHI received.  
3. **S5 (FAR ≤ 2/24h):** **Not met** and **not claimed** on public controls.  
4. **Ethics gates open:** checklist A1–A3, A5–A6, B6 still **N** until a real partner is engaged.  
5. **Partner not named in-repo:** by design — do not invent institutional contacts.

---

## 4. Recommendation (highest scientific value next)

| Option | Value | Recommendation |
|--------|-------|----------------|
| A. More public healthy Holter | Marginal; same mismatch class | **Deprioritize** as primary arm |
| B. Institutional / device-matched (Tier A) | Only path to realistic clinical FAR | **Active path** — request draft ready |
| C. Exploratory fusion / basal redesign | Pre-registered & labeled only | Not primary |

**Concrete next human step:** Identify **one** real partner in classes P1–P5 (hospital Holter lab, telemetry archive, EP research service, or device-core lab), complete checklist §A–B with their IRB/DUA, then export de-ID RR/annotations for n≈10–20 controls. Optionally refine de-ID fields with their privacy office before first export.

Public full-NSRDB FAR remains a **transparent interim reference** under the device-mismatch caveat — not evidence of clinical specificity.

---

## 5. Frozen parameters echoed (must match discovery / Phase 1)

```
θ₃ = 0.08
high-threshold = 0.65
W_TAU = 101
stride = 5
z_threshold = 2.0
min_consecutive = 3
alarm_rule = abs_z_from_basal >= z_threshold sustained min_consecutive
RR clean = [250, 2000] ms
```

See `results/external_phase2_far.json` → `params` and `s5_claim: false`, `clinical_claim: false`.

---

## 6. Key documents

| Artifact | Role |
|----------|------|
| `docs/PHASE2_INSTITUTIONAL_DATA_REQUEST.md` | **Partner-facing Tier A request brief** |
| `docs/PHASE2_IRB_DATA_CHECKLIST.md` | IRB / ethics / de-ID checklist |
| `docs/EXTERNAL_VALIDATION_PHASE2_PLAN.md` | Phase 2 plan (Tier A priority) |
| `results/external_phase2_far.json` | Public interim FAR |
| `results/phase2_public_control_inventory.csv` | Controls inventory (Tier A empty of PHI) |

---

## 7. Reproduce public interim FAR (unchanged)

```bash
cd /Users/johelpadilla/grok-safe/Investigaciones/Cardiac_CCTP_Pilot
python3 code/download_nsrdb_records.py --remaining
python3 code/extract_rr_external.py --db nsrdb
python3 code/run_external_validation_phase2_far.py
python3 -m pytest tests/test_phase2_planning_artifacts.py tests/test_phase2_far_artifacts.py tests/test_far_and_short_windows.py -q
```

---

## 8. Next-session prompt (copy-paste)

```
Continue External Validation Phase 2 — Tier A institutional path.
Base: v1.1.0 + public interim NSRDB FAR (n=18, τ_s ~33.7/24h, excess3 ~32.3/24h, device mismatch).
Partner brief ready: docs/PHASE2_INSTITUTIONAL_DATA_REQUEST.md
Checklist: docs/PHASE2_IRB_DATA_CHECKLIST.md
FROZEN: θ₃=0.08, high=0.65, W_TAU=101, abs-z≥2×3 — do not retune.
Priority: (1) identify one real clinical partner in classes P1–P5,
(2) complete checklist §A–B with their IRB/DUA language,
(3) refine de-ID export schema with their privacy office if needed,
(4) do NOT invent institutional inventory rows or claim S5.
Do NOT expand public healthy Holter as the main arm. No clinical deployability claims.
Root: /Users/johelpadilla/grok-safe/Investigaciones/Cardiac_CCTP_Pilot
```

---

## 9. Claims boundary

| Claim | Status |
|-------|--------|
| External sensitivity (Phase 1 VFDB) | Prior work only |
| Public interim control FAR (NSRDB n=18) | **Reported** with mismatch caveat |
| S5 met | **False / not claimed** |
| Clinical deployability / FDA | **None** |
| Device-matched institutional FAR | **Not yet available** (request drafted only) |
| Institutional partner named / DUA signed | **No** (honest empty) |
