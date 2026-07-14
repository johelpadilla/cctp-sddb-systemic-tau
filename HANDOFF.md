# CCTP — Handoff (authoritative status)

**Date**: 2026-07-14  
**Root**: `/Users/johelpadilla/grok-safe/Investigaciones/Cardiac_CCTP_Pilot/`  
**GitHub**: https://github.com/johelpadilla/cctp-sddb-systemic-tau  
**Last shipped release**: **v1.4.0** — final ordinal epilogue package (I0 + structural + ranking)  
**Zenodo concept DOI**: https://doi.org/10.5281/zenodo.21270698  
**Zenodo version DOI (v1.4.0)**: *pinned after archive* (see post-release pin commit)  
**Prior version DOI (v1.3.0)**: https://doi.org/10.5281/zenodo.21344730  

**Status snapshot**

| Workstream | State |
|------------|--------|
| SDDB discovery N=10 (stratified + lead-time + H2H) | **Done** (shipped) |
| Manuscript Phase 1 §3.9 + Phase 2 §3.10 | **Done** (shipped) |
| Manuscript **§3.11** + I0/structural + final ranking + §4.1 + §6 | **Done** — **v1.4.0** |
| External Validation Phase 1 / Phase 2 public FAR | **Done & reported** |
| Phase 2 plan / Tier A path | **Prepared** — main institutional next step |
| GitHub + Zenodo | **v1.4.0 release** |
| Native ordinal OPC/SDD bake-off + FAR + trade-off + cascade + grid | **Done** (v1.3.0+) |
| Integrated OPSP + I0 Holter + θ-grid 0/16 | **Done** — do not promote I0 |
| Structural surplus arms R0–R5 | **Done** — 0/6; **stop surplus-primary** |
| Final detector ranking (abs-z primary) | **Done** — in manuscript + `docs/FINAL_DETECTOR_RANKING.md` |
| Dual-mode C-BR / Mode-S multi-config Holter bake-off | Optional / low value after S1 stop |
| Clinical Copilot | Draft only — **do not deploy** |

**Frozen production params** (never retune on validation without pre-registration):  
θ₃=0.08, high-threshold=0.65, W_τ=101, W_EWS=501, stride=5, relative λ, detector **abs-z≥2 × 3** consecutive windows, RR clean [250, 2000] ms.

**OPC companion**: \(L=50\), \(\theta_D=0.35\), \(\theta_R=5\), joint K=36 — **keep_baseline**.

**Clinical / FDA / deployability claim**: **NONE**.

---

## Final detector ranking (locked)

| Rank | Arm | Role |
|------|-----|------|
| **1** | **abs-z \(\tau_s\)** | Preferred **primary** for pre-VF event hit rate (~0.91 sens) |
| 2 | SDD | Sensitivity ceiling (~0.97); high FAR |
| 3 | I0 surplus | Mode-S narrative only; **not promoted** |
| 4 | OPC L=50 | Specificity companion (~3.73 FAR); low sens |
| 5 | I-confirm | Secondary FAR filter only |

See: `docs/FINAL_DETECTOR_RANKING.md`, manuscript §3.11 / §4.1 / §6.

---

## Closed tracks

| Track | Outcome |
|-------|---------|
| I0 θ-grid (D′) | **0/16** clear advances |
| Structural surplus (S1) | **0/6** structural wins; stop surplus-primary |
| Promote I0 or OPC over abs-z | **Rejected** for event hit rate |

---

## Verify after clone

```bash
cd /path/to/cctp-sddb-systemic-tau
PYTHONPATH=code python3 -m pytest \
  tests/test_opc_refinements.py \
  tests/test_opsp_holter_eval_runner.py \
  tests/test_i0_surplus_param_grid.py \
  tests/test_i0_structural_arms.py \
  tests/test_ordinal_detectors.py -q
```

---

## Next scientific track

| Track | Action | Priority |
|-------|--------|----------|
| **B** | Institutional Tier A / Phase-2 device-matched non-event controls | **Main** |
| **C** | Dual-mode Holter (optional; low value after S1) | Optional |

Do **not** re-open I0 θ-grids or structural surplus micro-variants without new data.

---

## Source-of-truth documents

| Artifact | Role |
|----------|------|
| `manuscript/CCTP_SDBB_manuscript.md` (+ PDF) | Final report body |
| `docs/FINAL_DETECTOR_RANKING.md` | Ranking one-pager |
| `docs/I0_STRUCTURAL_ARMS.md` | R0–R5 + stop |
| `docs/I0_SURPLUS_PARAM_GRID.md` | 0/16 θ-grid |
| `docs/OPSP_INTEGRATED_HOLTER_EVAL.md` | I0 Holter eval |
| `docs/ORDINAL_SENSITIVITY_SPECIFICITY_TRADEOFF.md` | OPC/SDD/abs-z anchors |
