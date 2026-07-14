# OPSP integrated detector — Holter evaluation (I0 primary)

**Status:** Exploratory multi-cohort evaluation (not production; no clinical claim)  
**Date:** 2026-07-14  
**Root:** Cardiac CCTP pilot  
**Detector:** `opsp_integrated_detect` (`code/ordinal_detectors/opc_refinements.py`)  
**Runner:** `code/run_opsp_integrated_holter_eval.py`  
**Design:** `docs/INTEGRATED_N2_PERSIST_N3_SURPLUS_DETECTOR.md`  
**Artifacts:** `results/opsp_integrated_*.csv`, `results/opsp_integrated_summary.json`

---

## 0. Question

Does the integrated detector—**Nivel-2 persistence applied to Nivel-3 synergistic surplus**, with collapse only complementary—offer a better sensitivity–specificity (FAR) balance than existing arms (OPC companion, Mode S / OPS, frozen abs-z), while keeping high-synergy reorganizations detectable?

**Primary configuration evaluated: I0** (`collapse_role="none"`, \(L=50\), \(\theta_{\Delta S}=0.08\), \(\theta_R^S=5\)).

Secondary (collapse subordinate only): **I-mod**, **I-confirm**.  
Mode-S twin: **OPS-S1** (`ops_detect` same basal-relative gates) for core-identity check.

---

## 1. Protocol (same family as prior ordinal tables)

| Axis | Definition |
|------|------------|
| Encoding | Bivariate Bandt–Pompe \(m=3\), factors \(\pi^{(1)},\pi^{(2)}\in\{0..5\}\); surplus \(S_t=\mathrm{TV}(P_{\mathrm{joint}},P_1\otimes P_2)\) |
| Sensitivity | First alarm in post-basal, pre-event search window (SDDB analytic windows; VFDB short-DB thirds) |
| Cohorts | SDDB \(n=11\), VFDB \(n=22\), total events \(n=33\) |
| FAR | NSRDB \(n=18\), Phase-2 basal \((0.25,\min(2,0.25\cdot T))\), search remainder, cap 12 h; episode refractory 0.5 h; \(\mathrm{FAR}=N_{\mathrm{ep}}/H_{\mathrm{search}}\times 24\) |
| Baselines | OPC L=50 companion and abs-z from existing same-protocol results (not retuned) |

**Proxy honesty:** \(S_t\) is a **Level-3–consistent ordinal proxy**, not continuous excess3 / nested \(\Phi_3\).

---

## 2. Results — I0 and secondary roles (this run)

### 2.1 Sensitivity

| Detector | Role | Sens SDDB | Sens VFDB | Sens all | Detected / \(n\) | Median lead (all, h) |
|----------|------|-----------|-----------|----------|------------------|----------------------|
| **I0** | **Primary** | **1.000** (11/11) | **0.682** (15/22) | **0.788** (26/33) | 26/33 | 0.285 |
| I-mod | Secondary | 1.000 (11/11) | 0.682 (15/22) | 0.788 (26/33) | 26/33 | 0.285 |
| I-confirm | Secondary | 0.545 (6/11) | 0.136 (3/22) | 0.273 (9/33) | 9/33 | 1.657 |
| OPS-S1 | Mode-S twin | 1.000 (11/11) | 0.682 (15/22) | 0.788 (26/33) | 26/33 | 0.285 |

### 2.2 FAR (NSRDB)

| Detector | FAR / 24h | Episodes | Search hours | Fraction controls alarmed |
|----------|-----------|----------|--------------|---------------------------|
| **I0** | **40.134** | 301 | 179.998 | 1.000 |
| I-mod | 40.134 | 301 | 179.998 | 1.000 |
| I-confirm | **1.733** | 13 | 179.998 | 0.278 |
| OPS-S1 | 40.134 | 301 | 179.998 | 1.000 |

### 2.3 Secondary-role impact vs I0

| Variant | \(\Delta\) sens (all) | \(\Delta\) FAR /24h | Interpretation |
|---------|----------------------|---------------------|----------------|
| I-mod | 0 (identical alarm set) | 0 | Collapse-assisted gate relaxation did **not** add Holter detections or FAR on this cohort; I0 already fires when surplus-persist holds. |
| I-confirm | **−0.515** (26→9) | **−38.4** (40.1→1.73) | Strong FAR cut via collapse co-occurrence + stricter surplus gates, at **severe** sensitivity cost; pure high-support surplus paths largely silenced. |

I0 alarms ≡ OPS-S1 alarms on all 33 event records (lead times match). That confirms the design claim: **I0 core is surplus-persist**, not a new collapse hybrid.

---

## 3. Comparison vs baselines (same protocol family)

| Arm | Sens all | FAR /24h | Notes |
|-----|----------|----------|-------|
| **OPC companion L=50** | 0.424 (14/33) | **3.733** | Collapse-persist; L3-blind by construction |
| **I0 (this run)** | **0.788** (26/33) | 40.134 | Surplus-persist; L3 path open |
| **OPS-S1 (this run)** | 0.788 (26/33) | 40.134 | = I0 core (expected) |
| **I-confirm (secondary)** | 0.273 (9/33) | 1.733 | FAR-leaning filter; not RECD default |
| **abs-z \(\tau_s\) frozen** | **0.909** (30/33) | 33.734 | Continuous production reference |

### 3.1 Does I0 improve the sens–FAR balance?

| Contrast | Sensitivity | FAR | Honest judgment |
|----------|-------------|-----|-----------------|
| I0 vs OPC | **Better** (+0.36 absolute) | **Much worse** (~10.8×) | Higher sens does **not** compensate for FAR saturation on NSRDB. OPC remains the better **specificity** ordinal arm. |
| I0 vs OPS standalone | **Equal** | **Equal** | I0 is not a new surplus detector; it *is* Mode-S S1 under another name when collapse is ignored. |
| I0 vs abs-z | **Worse** (−0.12) | **Worse** (+6.4/24h) | No win on either axis vs frozen production reference. |
| I-confirm vs OPC | Worse sens (−0.15) | Better FAR (1.73 vs 3.73) | FAR win is real but sensitivity collapses below OPC; not a free lunch. |

**Bottom line (measured, not narrative):**  
I0 **does not** currently offer a better sensitivity–specificity trade-off than prior arms. It recovers more event detections than collapse-only OPC and keeps the high-surplus path open (L3 intent), but **FAR saturates** in the same high band as OPS / abs-z-class rates. It is **not** ready to replace OPC as the project’s primary exploratory ordinal specificity arm, nor abs-z as the continuous production reference.

---

## 4. Nivel-3 respect verification

| Criterion | Result |
|-----------|--------|
| Primary predicate = surplus persistence | **Pass** (I0; Holter I0 ≡ OPS-S1) |
| Collapse-only ⇒ no alarm | **Pass** (unit suite all roles) |
| High-support + high-surplus still detectable on I0 / I-mod | **Pass** (units + Holter I0 detections without requiring collapse) |
| \(S_t\) labeled proxy ≠ excess3 | **Pass** |
| No \(\mu/\sigma\) / abs-z core | **Pass** |
| I-confirm L3 cost disclosed | **Pass** (sens 0.27; secondary only) |
| Free-OR / collapse-dominant not sold as primary | **Pass** |

Unit evidence: `PYTHONPATH=code pytest tests/test_opc_refinements.py` → **17 passed**.

### 4.1 Known gap (by design)

Pure **collapse with low surplus** may stay silent under I0. Overlap with OPC L=50:

| Pattern | Count (of 33) |
|---------|---------------|
| Both alarm | 12 |
| I0 only | 14 |
| OPC only (collapse without surplus-persist) | **2** — `vfdb/419`, `vfdb/430` |
| Neither | 5 |

Those two OPC-only records are the empirical face of the intentional hierarchy: collapse is not dominant. Keep Mode-C / OPC as a **separate** report row if collapse-only events matter.

---

## 5. Strengths, limitations, gaps

### Strengths
- Ontology matches the requested N2×N3 story: persistence *logic* on surplus *content*.
- L3 high-support path is open under I0 (and I-mod).
- Implementation invariants hold under unit tests; Holter runner is reproducible.
- I0 ≡ OPS-S1 on real data (no accidental collapse leak into the primary alarm).
- I-confirm demonstrates that collapse-as-filter can cut FAR when desired—as a **secondary** tool.

### Limitations
- **FAR ≈ 40/24h** on NSRDB with I0/I-mod/OPS-S1 — unusable as a specificity-facing primary arm without retuning surplus gates or basal control.
- I-mod adds no Holter benefit at current \(\delta_R,\delta_S\) (identical to I0).
- I-confirm’s FAR win destroys most sensitivity (including many high-surplus paths without nearby locking).
- NSRDB is not device-matched to VFDB/SDDB (same caveat as prior FAR tables).
- \(S_t\) remains a TV independence proxy, not continuous RECD excess3.

### Gaps / refinements still needed (if pursuing OPSP further)
1. **I0-strict** grid: \(\theta_{\Delta S}\in\{0.10,0.12,0.15\}\), \(\theta_R^S\in\{8,10,12\}\) with full sens/FAR tables (doc already names I0-strict as specificity-leaning).
2. Basal-relative surplus with **tighter basal estimation** or longer basal on Holter (reduce basal underestimation → false \(\Delta S\)).
3. Keep **separate Mode C** row for collapse-only events; do not force one detector to own both classes.
4. Optional: evaluate absolute \(S_t\) gates only as diagnostic (doc: high FAR risk).
5. Do **not** market free \(A_C\vee A_S\) as “integrated primary.”

---

## 6. Recommendation

| Decision | Choice |
|----------|--------|
| **Promote I0 as the project’s primary exploratory ordinal arm?** | **No** |
| **Why (measured)** | FAR ~40/24h with sens 0.79 is a worse specificity balance than OPC (sens 0.42, FAR 3.73) and does not beat abs-z on either axis. I0 is ontologically Mode-S S1, not a new superior hybrid. |
| **What to keep** | OPSP code + tests + this evaluation as the **surplus-primary design track**; I0 remains the **recommended starting integrated config for further gate tuning**, not the deployed primary comparison arm. |
| **What stays primary for tables** | OPC companion (specificity ordinal) + abs-z (continuous production) + optional separate Mode-S / I0 row for surplus narrative. |
| **Next refinement if continuing** | Run **I0-strict** / surplus-threshold grid on the same Holter protocol; promote only if a point improves sens–FAR relative to OPC without silencing high-support surplus (L3). |

**Clinical / superiority claims:** none.

---

## 7. Reproducibility

```bash
cd Investigaciones/Cardiac_CCTP_Pilot
PYTHONPATH=code python3 -m pytest tests/test_opc_refinements.py -q
python3 code/run_opsp_integrated_holter_eval.py
# outputs: results/opsp_integrated_sens_per_record.csv
#          results/opsp_integrated_nsrdb_far_per_record.csv
#          results/opsp_integrated_summary.json
#          results/opsp_integrated_comparison.csv
```

Smoke path: `python3 code/run_opsp_integrated_holter_eval.py --smoke` (subset only; does not replace full tables).
