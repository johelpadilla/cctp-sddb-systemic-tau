# Pre-registered I0 surplus-persist parameter grid vs OPC

**Status:** Exploratory grid complete (not production; no clinical claim)  
**Date:** 2026-07-14  
**Root:** Cardiac CCTP pilot  
**Runner:** `code/run_i0_surplus_param_grid.py`  
**Detector:** `opsp_integrated_detect` with `collapse_role="none"` (I0 surplus-primary only)  
**Artifacts:** `results/i0_surplus_param_grid*.csv`, `results/i0_surplus_param_grid_summary.json`  
**Prior eval:** `docs/OPSP_INTEGRATED_HOLTER_EVAL.md` (I0 default θ_ΔS=0.08, θ_R=5)

---

## 0. Question

Does any configuration of the integrated I0 detector (Nivel-2 persistence on Nivel-3 synergistic surplus) achieve a **clearly better** sensitivity–specificity balance than the OPC companion, **without** sacrificing the high-synergy (Nivel-3) path?

Prior single-point I0 (0.08, 5) had sens_all = 0.788 but FAR ≈ 40.1/24h (~10.8× OPC). This note runs a **small, pre-registered** grid and scores every cell with the **same** success rule written **before** judgment.

---

## 1. Pre-registered success criteria (locked before grid judgment)

OPC companion anchors (same protocol family; not re-estimated here):

| Anchor | Value |
|--------|-------|
| OPC sens_all | **0.424** (14/33) |
| OPC FAR /24h | **3.733** |

**Clear advance** requires **both**:

1. **FAR/24h ≤ 2 × FAR_OPC** → cap **≤ 7.467**/24h  
2. **sens_all ≥ 0.65**

| Outcome | Counts as success? |
|---------|-------------------|
| FAR ≤ 2×OPC **and** sens ≥ 0.65 | **Yes — clear advance** |
| FAR ≤ 2×OPC but sens in [OPC, 0.65) | **No** (one-sided / marginal) |
| sens ≥ 0.65 but FAR > 2×OPC | **No** (I0 default regime) |
| Any move still in the high-FAR band without both bars | **No** |
| I-confirm / free OR / collapse-filter as judged arm | **No** (out of scope for this grid) |

**Judged arm only:** I0 surplus-persist (`collapse_role=none`). Collapse-filter roles are not substitutes for I0 success.

**Rationale for the bars:**  
- FAR ≤ 2×OPC is a hard specificity ceiling so “better balance” is not mere high-sens noise.  
- sens ≥ 0.65 is well above OPC (0.42) so a FAR-controlled cell must still show a **real** sensitivity gain, not a wash.  
- Both required; no story-telling around near-misses.

Encoded in code: `PRE_REGISTERED_RULES` + `score_clear_advance` in `run_i0_surplus_param_grid.py` (unit-tested).

---

## 2. Grid specification

| Axis | Values | Fixed |
|------|--------|-------|
| θ_ΔS | {0.08, 0.10, 0.12, 0.15} | — |
| θ_R^S | {5, 8, 10, 12} | — |
| L | — | **50** |
| collapse_role | — | **none** |

**16 cells.** Product order: θ_ΔS outer, θ_R inner.

### Protocol (identical family to OPSP Holter / ordinal bake-off)

| Axis | Definition |
|------|------------|
| Encoding | Bivariate Bandt–Pompe m=3; \(S_t=\mathrm{TV}(P_{\mathrm{joint}},P_1\otimes P_2)\) |
| Sensitivity | First alarm in post-basal, pre-event window; SDDB n=11 + VFDB n=22 = **33** |
| FAR | NSRDB n=18, Phase-2 basal, search remainder, cap 12 h; episode refractory 0.5 h; FAR = episodes/search_h × 24 |
| Proxy honesty | \(S_t\) is Level-3–consistent ordinal proxy, **not** continuous excess3 |

---

## 3. Results (full non-smoke Holter run)

### 3.1 Full grid table

| θ_ΔS | θ_R | Sens SDDB | Sens VFDB | **Sens all** | Det / n | **FAR /24h** | FAR / OPC | **Clear advance?** |
|------|-----|-----------|-----------|--------------|---------|--------------|-----------|-------------------|
| 0.08 | 5 | 1.000 | 0.682 | **0.788** | 26/33 | **40.134** | 10.75 | **no** |
| 0.08 | 8 | 1.000 | 0.545 | **0.697** | 23/33 | **37.334** | 10.00 | **no** |
| 0.08 | 10 | 1.000 | 0.455 | 0.636 | 21/33 | 35.467 | 9.50 | no |
| 0.08 | 12 | 1.000 | 0.409 | 0.606 | 20/33 | 33.867 | 9.07 | no |
| 0.10 | 5 | 1.000 | 0.409 | 0.606 | 20/33 | 34.134 | 9.14 | no |
| 0.10 | 8 | 1.000 | 0.364 | 0.576 | 19/33 | 30.934 | 8.29 | no |
| 0.10 | 10 | 1.000 | 0.364 | 0.576 | 19/33 | 27.467 | 7.36 | no |
| 0.10 | 12 | 1.000 | 0.364 | 0.576 | 19/33 | 25.067 | 6.71 | no |
| 0.12 | 5 | 0.909 | 0.273 | 0.485 | 16/33 | 25.600 | 6.86 | no |
| 0.12 | 8 | 0.727 | 0.227 | 0.394 | 13/33 | 21.467 | 5.75 | no |
| 0.12 | 10 | 0.727 | 0.227 | 0.394 | 13/33 | 19.200 | 5.14 | no |
| 0.12 | 12 | 0.727 | 0.227 | 0.394 | 13/33 | 16.400 | 4.39 | no |
| 0.15 | 5 | 0.727 | 0.182 | 0.364 | 12/33 | 13.467 | 3.61 | no |
| 0.15 | 8 | 0.727 | 0.182 | 0.364 | 12/33 | 10.533 | 2.82 | no |
| 0.15 | 10 | 0.727 | 0.182 | 0.364 | 12/33 | 9.067 | 2.43 | no |
| 0.15 | 12 | 0.727 | 0.182 | 0.364 | 12/33 | **7.733** | **2.07** | **no** |

Machine-readable: `results/i0_surplus_param_grid.csv`, `results/i0_surplus_param_grid_summary.json`.

### 3.2 Comparison vs OPC

| Arm | Sens all | FAR /24h | vs clear-advance rule |
|-----|----------|----------|------------------------|
| **OPC L=50 companion** | 0.424 | **3.733** | baseline |
| I0 best sens (0.08, 5) | 0.788 | 40.134 | fails FAR (10.8×) |
| I0 lowest FAR (0.15, 12) | 0.364 | 7.733 | fails FAR (still > 2×) **and** fails sens (< OPC) |
| Cells with sens ≥ 0.65 | 2 of 16 | FAR 37–40 | both one-sided high-FAR |
| Cells with FAR ≤ 7.467 | **0 of 16** | — | no cell enters the FAR bar |
| **Clear-advance count** | — | — | **0 / 16** |

### 3.3 Trade-off geometry (honest read)

```
sens ↑
 0.79 │  ● (0.08,5)  ● (0.08,8)     ← high sens, FAR~37–40
 0.65 │ ──────── sens bar ────────
 0.60 │  ●●●● (θ_ΔS 0.08–0.10 mid)  FAR still ~25–35
 0.42 │ ── OPC ───────────────────
 0.36 │              ●●●● (0.15,*)  ← lowest FAR still ~7.7–13
      └──────────────────────────────────────── FAR →
         3.7 OPC   7.5 (2× bar)        25–40
```

Monotone pattern (as expected for stricter gates):

- Raising θ_R or θ_ΔS **reduces** both FAR and sensitivity.
- The path never crosses into the rectangle **(FAR ≤ 7.47) ∩ (sens ≥ 0.65)**.
- Best FAR cell (0.15, 12) is **just above** the FAR bar (7.733 vs 7.467) with sens **below** OPC (0.364 < 0.424).

---

## 4. Evaluation against pre-defined criteria

| Criterion | Result |
|-----------|--------|
| Rules written before judgment | **Yes** (this §1 + code constants + scorer tests) |
| Full 4×4 on real Holter protocol | **Yes** (33 events, 18 NSRDB, non-smoke) |
| Per-cell boolean clear_advance | **Yes** — all false |
| Any clear advance? | **No — 0 / 16** |
| Surplus-primary only (not confirm/OR) | **Yes** — all cells `collapse_role=none` |
| High-synergy path not sacrificed as “win” | **Yes** — unit/synthetic: strict corner (0.15, 12) still alarms on dependent high-surplus stream; Holter silence on some VFDB records is gate severity, not collapse-filter |

**Conclusion on criteria:**  
Under this proxy and protocol, **surplus-persist I0 has a practical ceiling**: no grid point is a clear advance vs OPC. High sensitivity lives only in a high-FAR band; FAR control collapses sensitivity to **at or below** OPC while still missing the 2×OPC FAR bar at the strictest corner.

---

## 5. Nivel-3 / high-synergy check

| Check | Result |
|-------|--------|
| Primary predicate = surplus persistence | Pass (I0 only) |
| Collapse never primary | Pass |
| Strict corner alarms on synthetic high-surplus stream | Pass (`tests/test_i0_surplus_param_grid.py`) |
| Proxy labeled ≠ excess3 | Pass |
| No μ/σ core | Pass |

Stricter gates reduce Holter detections (expected) but do **not** redefine the ontology around collapse. L3 path remains open when surplus is strong enough; the FAR problem is that **NSRDB also produces frequent surplus-persist runs**.

---

## 6. Strategic recommendation

| Decision | Choice |
|----------|--------|
| **Is there a clear-advance I0 cell?** | **No (0/16)** |
| **Promote any I0 grid point as primary ordinal specificity arm?** | **No** |
| **Strategic path** | **Return to OPC as the primary ordinal specificity arm** for tables; treat I0/OPS as a **surplus narrative / Mode-S track**, not a FAR-competitive replacement. If surplus work continues, prefer **structural** changes (better basal estimation, different surplus statistic, longer L, multi-scale, or a true nested Φ₃ proxy)—**not** another modest θ_ΔS/θ_R grid of the same I0 form. |
| **Do not** | Market I-confirm FAR wins as I0 success; expand this grid hoping for a free lunch; claim OPS ≡ excess3 |

### Why not “keep refining I0 thresholds”?

This grid already spans the named I0-default and I0-strict neighborhoods from the design note. The frontier is smooth: FAR and sens move together. The nearest FAR-bar approach (0.15, 12) **underperforms OPC on sensitivity** and **still fails** the FAR bar. Further micro-tuning θ_ΔS/θ_R is unlikely to open the empty (low-FAR, high-sens) rectangle without a **different signal or basal design**.

### What remains useful

- I0 / OPS code + Holter tables as the **Level-3–consistent surplus path** documentation.  
- OPC companion for **specificity-facing** ordinal comparison.  
- abs-z for continuous production reference.  
- Optional dual-mode reporting (Mode C / Mode S) rather than forcing one detector to own both collapse-only and high-synergy classes.

**Clinical / superiority claims:** none.

---

## 7. Reproducibility

```bash
cd Investigaciones/Cardiac_CCTP_Pilot
PYTHONPATH=code python3 -m pytest tests/test_i0_surplus_param_grid.py tests/test_opc_refinements.py -q
PYTHONPATH=code python3 code/run_i0_surplus_param_grid.py
# outputs:
#   results/i0_surplus_param_grid_sens_per_record.csv
#   results/i0_surplus_param_grid_nsrdb_far_per_record.csv
#   results/i0_surplus_param_grid.csv
#   results/i0_surplus_param_grid_summary.json
#   docs/I0_SURPLUS_PARAM_GRID.md  (this file)
```

Smoke (wiring only; not gating):  
`PYTHONPATH=code python3 code/run_i0_surplus_param_grid.py --smoke`
