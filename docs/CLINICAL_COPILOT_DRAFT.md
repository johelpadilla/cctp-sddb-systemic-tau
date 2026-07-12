# Clinical Copilot — First Draft (Research Use Only)

**Name (working):** CCTP Relational Risk Copilot  
**Version:** 0.1 draft specification  
**Status:** Design / offline prototype scope — **not** a diagnostic device  

---

## 1. Purpose

Provide clinicians and researchers with a **transparent, research-grade** view of relational heart-rate dynamics (Systemic Tau τ_s and RECD excess3) in the hours before a known or candidate malignant event, alongside classic early-warning statistics (variance, AR(1)).

The Copilot is a **decision-support visualization and audit trail**, not an automated VF alarm system.

---

## 2. Explicit non-diagnostic / research-use limits

- **Not cleared** by FDA, EMA, or any regulatory body.  
- **Not** for real-time bedside alarming or therapy titration.  
- Must display a persistent banner:  
  > “Research prototype. Not for clinical diagnosis or emergency decision-making. Relational metrics are experimental.”  
- No automatic ICD programming, no drug dosing recommendations, no “risk score” marketed as calibrated probability of VF without external validation (see `EXTERNAL_VALIDATION_PLAN.md`).  
- Outputs are **descriptive** (time series + deltas + lead-time under a frozen rule), not probability-calibrated.

---

## 3. Inputs

| Input | Format | Required |
|-------|--------|----------|
| RR interval series | ms, continuous or beat-indexed | Yes |
| Time base | seconds from recording start or absolute timestamps | Yes |
| Event marker (VF/VT onset) | seconds / timestamp | Optional for offline review; required for labeled lead-time audit |
| Beat quality / annotation flags | optional | Preferred |
| Clinical context | age, substrate (sinus/AF/paced), pacing mode | Optional metadata |

**Ingestion paths (v0.1):**
1. Pre-cleaned `rr_*_clean.npz` (paper pipeline)  
2. PhysioNet-style annotations via `extract_rr.py`  
3. CSV: `t_sec, rr_ms` (+ optional `event_sec`)

---

## 4. Primary computations (frozen paper stack)

1. Clean RR (250–2000 ms; report interp fraction).  
2. Proxy \(X=[z(\mathrm{RR}), z(|\Delta RR|)]\).  
3. τ_s (W=101, stride=5) via `systemictau`.  
4. excess3 (θ₃=0.08) via nested RECD.  
5. Variance & AR(1) (W=501, stride=5).  
6. Basal vs approach (last 3 h) deltas.  
7. Lead-time under abs-z ≥ 2 vs basal, ≥3 consecutive samples.

All parameters exposed as read-only “protocol card” in the UI.

---

## 5. Primary displays

### 5.1 Overview strip
- Recording duration, n beats, interp %, pacing flag, substrate.  
- Event line if known.

### 5.2 Relational panel (primary)
- τ_s(t) with basal / approach shading.  
- excess3(t) aligned.  
- Numeric Δτ_s, Δexcess3, Welch or surrogate p when available.  
- Sign concordance indicator (τ_s vs excess3).

### 5.3 Classic EWS panel (comparator)
- Rolling variance and AR(1).  
- Δvar, ΔAR(1) same windows.  
- Explicit note: univariate EWS may **disagree in sign** with relational metrics.

### 5.4 Detector / lead-time panel
- Detection time under frozen z-rule for each metric.  
- Lead-time (hours) to event.  
- Cumulative detection context (cohort reference, not patient-specific probability).

### 5.5 Risk framing language (allowed wording)

| Allowed | Forbidden |
|---------|-----------|
| “Relational reorganization relative to this recording’s basal window” | “Patient will have VF in X minutes” |
| “Δτ_s / Δexcess3 consistent with pre-event change seen in SDDB cohort (N=10)” | “High/low risk score 0–100 calibrated” |
| “Research metric crossed experimental threshold” | “Alarm — treat now” |

---

## 6. Minimal offline mock (implementation path)

A first local prototype can be a single script or notebook:

```text
python3 code/analyze_cctp_pilot.py --record 38
python3 code/run_leadtime_detector.py --records 38
# → figures/{record}/ + results JSON/CSV
```

Optional later: Streamlit/Gradio app reading one npz and rendering the five panels above. **Out of scope for this draft:** auth, multi-user EHR, cloud hosting, PHI storage.

---

## 7. Data & privacy

- Prefer de-identified Holter exports.  
- No automatic upload of clinical data in v0.1.  
- Log analysis commit hash + parameter card with every exported report.

---

## 8. Acceptance criteria for “Copilot v0.1 complete”

- [x] Spec documents inputs, displays, and research-use limits (this file).  
- [ ] Optional: one-command local HTML/PNG report for a single record.  
- [ ] Optional: side-by-side cohort reference strip from `results/cctp_cohort_stratified.csv`.  
- [ ] External validation gates from `EXTERNAL_VALIDATION_PLAN.md` before any prospective UI claim.

---

## 9. Ownership & citation

Relational metrics: Systemic Tau (`systemictau`) + nested RECD as in the CCTP/SDDB manuscript.  
Always cite the manuscript DOI / Zenodo code release when distributing Copilot screenshots or reports.
