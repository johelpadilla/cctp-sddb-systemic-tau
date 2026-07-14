# I0 structural surplus arms (Track S1) — Holter results

**Status:** Complete (exploratory; no clinical claim)  
**Date:** 2026-07-14  
**Runner:** `code/run_i0_structural_arms.py`  
**Detector API:** `opsp_integrated_detect` / `ops_detect` with `basal_mode`, `hop`, optional larger `L`  
**Artifacts:** `results/i0_structural_arms_*.csv` / `_summary.json`  
**Prior:** I0 θ-grid 0/16 clear advances (`docs/I0_SURPLUS_PARAM_GRID.md`) — **no further θ-grid**

---

## 0. Question

Can **structural** changes (robust basal, hop-decoupled persistence, larger \(L_S\)) move I0 surplus-persist into a usable sens–FAR region vs OPC, without collapse as the alarm engine?

---

## 1. Pre-registered arms (fixed 6; no Cartesian product)

| Arm | \(L_S\) | Basal / thr | Persist |
|-----|---------|-------------|---------|
| **R0** I0-ref | 50 | mean + 0.08 | hop=1, θ_R=5 |
| **R1** hop | 50 | mean + 0.08 | hop=25, θ_R=3 |
| **R2** basal-q90 | 50 | thr = q90(basal) | hop=1, θ_R=5 |
| **R3** basal-MAD | 50 | med + 2.5·MAD | hop=1, θ_R=5 |
| **R4** combo | 50 | thr = q90 | hop=25, θ_R=3 |
| **R5** combo+\(L_S\) | **100** | thr = q90 | hop=50, θ_R=3 |

All: `collapse_role=none`. Protocol = OPSP Holter (SDDB+VFDB n=33; NSRDB n=18 FAR; refractory 0.5 h).

### Success rules (locked before judgment)

| Rule | FAR | sens_all |
|------|-----|----------|
| **structural_win** | ≤ 2×OPC (≤ **7.467**) | ≥ **0.55** |
| **clear_advance** (jackpot) | ≤ 2×OPC | ≥ **0.65** |
| **approaches** | ≤ 12 | ≥ 0.55 |

**Stop rule:** if no arm structural_win and none approaches → stop surplus-primary as FAR competitor.

---

## 2. Measured Holter results (full non-smoke)

| Arm | Sens SDDB | Sens VFDB | **Sens all** | Det | **FAR/24h** | FAR/OPC | Struct? | Clear? |
|-----|-----------|-----------|--------------|-----|-------------|---------|---------|--------|
| **R0** I0-ref | ~1.00 | ~0.64 | **0.758** (25/33) | 25 | **40.134** | 10.75× | no | no |
| **R1** hop | ~0.73 | ~0.14 | **0.333** (11/33) | 11 | **10.933** | 2.93× | no | no |
| **R2** q90 | ~1.00 | ~0.68 | **0.788** (26/33) | 26 | **43.334** | 11.61× | no | no |
| **R3** MAD | ~0.91 | ~0.45 | **0.606** (20/33) | 20 | **23.334** | 6.25× | no | no |
| **R4** combo | ~0.91 | ~0.32 | **0.515** (17/33) | 17 | **17.067** | 4.57× | no | no |
| **R5** combo+Ls100 | ~0.82 | ~0.18 | **0.394** (13/33) | 13 | **11.467** | 3.07× | no | no |
| OPC L=50 (anchor) | — | — | 0.424 | 14 | 3.733 | 1× | — | — |

**Counts:** structural_wins **0/6**; clear_advance **0/6**; approaches **0/6**.  
**stop_surplus_primary_recommended: true**

Machine-readable: `results/i0_structural_arms_comparison.csv`, `results/i0_structural_arms_summary.json`.

---

## 3. Interpretation (honest)

### What worked directionally

1. **Hop (R1)** is a real FAR lever: 40 → **11**/24h (~2.9×OPC). Confirms that hop=1 θ_R was mostly **pseudo-evidence** from overlapping windows.
2. **MAD basal (R3)** is a real dual lever: FAR ~halved vs R0 (40 → 23) with sens still **0.61** — best single-mechanism trade-off among structural single changes.
3. **Combo R4/R5** continue the FAR descent (17 → 11) but **sens falls through the structural bar** (0.52 → 0.39).

### What failed / surprised

1. **q90 alone (R2) raised FAR** (43 vs 40). Absolute thr = percentile(basal) is often *looser* than mean+0.08 when basal mean is moderate and the +0.08 offset sits above the empirical tail. Not a free specificity upgrade.
2. **No arm enters FAR ≤ 2×OPC with sens ≥ 0.55.** Closest FAR arms kill VFDB sensitivity (short pre-event windows cannot accumulate 3 hop-credits).
3. **L_S=100 (R5)** does not unlock a new region: smoother \(S\) + longer hop still leaves FAR ~3×OPC with sens ≤ OPC.
4. Geometry remains a **smooth frontier** (same lesson as the θ-grid): structural knobs move *along* the curve, they do not open the empty rectangle.

### Diagnosis (updated)

| Failure mode | Status after S1 |
|--------------|-----------------|
| Pseudo-persist hop=1 | **Confirmed** — hop fixes FAR, breaks short events |
| Mean basal insensitive to tails | MAD helps; percentile alone miscalibrated vs mean+Δ |
| Sparse TV at L=50 | L=100 insufficient for Holter specificity |
| NSRDB also has surplus runs | **Still true** — robust thr + hop cannot reach OPC FAR band without sens collapse |

**Practical ceiling:** surplus-primary ordinal TV on public Holter is **not** a FAR-competitive replacement for OPC. Keep as Mode-S / L3-proxy **narrative** arm only.

---

## 4. Decision

| Action | Decision |
|--------|----------|
| Promote any R* as primary ordinal | **No** |
| Further hop / basal / \(L_S\) micro-variants | **No** (stop rule fired) |
| Another θ_ΔS×θ_R grid | **No** (already closed) |
| OPC as ordinal specificity arm | **Yes** (unchanged) |
| I0/OPS Mode-S narrative | **Optional** row; not primary table |
| Institutional Tier A / packaging | Preferred next (Track B / A′) |

---

## 5. Reproduce

```bash
cd /Users/johelpadilla/grok-safe/Investigaciones/Cardiac_CCTP_Pilot
PYTHONPATH=code python3 -m pytest tests/test_i0_structural_arms.py tests/test_opc_refinements.py -q
PYTHONPATH=code python3 code/run_i0_structural_arms.py
# optional smoke: --smoke
```

---

## 6. Code surface added

| Piece | Role |
|-------|------|
| `compute_surplus_threshold` | mean_delta / percentile / mad |
| `ops_detect` / `opsp_integrated_detect` | `basal_mode`, `basal_q`, `mad_kappa`, `hop` |
| `run_i0_structural_arms.py` | Holter bake-off R0–R5 + scorers |
| `tests/test_i0_structural_arms.py` | unit coverage |

Defaults preserve legacy I0 (`basal_mode=mean_delta`, `hop=1`).
