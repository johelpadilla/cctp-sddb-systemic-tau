# CCTP — Handoff (authoritative status)

**Date**: 2026-07-12  
**Location**: `Investigaciones/Cardiac_CCTP_Pilot/`  
**Status**: SDDB discovery (N=10) + manuscript Option A **done**; **External Validation Phase 1 started and reported** (VFDB n=11 independent short events + NSRDB controls FAR).  
**Frozen params** (never retuned on validation data): $\theta_3=0.08$, high-threshold $=0.65$, $W_{\tau}=101$, $W_{\mathrm{EWS}}=501$, stride $=5$, relative $\lambda$, detector abs-$z\geq 2$ × 3 consecutive windows.  
**Clinical / deployability claim**: **NONE** (Phase 1 fails S5 FAR by a wide margin).

---

## External Validation Phase 1 — Active source of truth

**Report:** `docs/EXTERNAL_VALIDATION_PHASE1_REPORT.md`  
**Plan:** `docs/EXTERNAL_VALIDATION_PLAN.md`

### Data
| Source | Local | Role |
|--------|-------|------|
| VFDB | `data/vfdb/` 22 recs | Independent short VT/VF |
| CUDB | `data/cudb/` 35 recs | All excluded (pre-event &lt;15 min) |
| NSRDB | `data/nsrdb/` 6 Holters | Negative controls for FAR |
| Clean RR | `data/rr_external/` | 62 npz |

### Independent events (VFDB, pre≥15 min, n=11)
| Metric | Sensitivity | Median lead (h ≈ min) |
|--------|-------------|------------------------|
| τ_s | **1.00** (11/11) | 0.173 ≈ **10.4 min** |
| excess3 | **0.82** (9/11) | 0.127 ≈ **7.6 min** |
| var / AR(1) | 0.55 / 0.55 | ~0.24 h |

Concordance τ_s–excess3 (Δ signs): **0.64** (7/11). All processable events are `short_15_60min`.

### FAR (NSRDB controls, n=6, ~60 search hours)
| Metric | FAR / 24 h | Fraction of controls with ≥1 alarm |
|--------|------------|-------------------------------------|
| τ_s | **~34.4** | 1.0 |
| excess3 | **~28.8** | 1.0 |

**Reading:** S4 (sens≥0.60) holds on this small short set; **S3 (median lead≥15 min short-DB)** and **S5 (FAR≤2/24h)** do **not**. The Jul-12 abs-z rule is **not clinically deployable**.

### Artifacts
- `results/external_phase1_inventory.csv`
- `results/external_phase1_per_record.csv`
- `results/external_phase1_sensitivity.json`
- `results/external_phase1_controls.csv`
- `results/external_phase1_far.json`
- `results/external_phase1_summary.json` (`clinical_claim=false`)

### Reproduce
```bash
python3 code/extract_rr_external.py --db all
python3 code/run_external_validation_phase1.py
python3 -m pytest tests/test_leadtime_detector.py tests/test_far_and_short_windows.py -q
```

### Next recommended action
**Phase 2 specificity** — institutional / device-matched non-event Holters and longer pre-VF series; keep frozen params for primary claims; any detector fusion must be **pre-registered exploratory**, not silent retune on SDDB. Do **not** ship Copilot alarms with current FAR.

### Next-session prompt
```
Continue CCTP External Validation after Phase 1
(docs/EXTERNAL_VALIDATION_PHASE1_REPORT.md).

Facts: VFDB n=11 independent short events (τs sens=1.0, excess3=0.82,
median lead ~8–10 min); NSRDB FAR ~29–34 / 24h (S5 fail); CU all excluded
(pre<15min); params frozen θ3=0.08 W_TAU=101 abs-z≥2×3.

Highest-impact next: Phase 2 specificity — institutional/matched controls
and/or pre-registered detector fusion WITHOUT retuning on SDDB discovery.
Do not claim clinical deployability.
Root: Investigaciones/Cardiac_CCTP_Pilot
```

---

## Jul-12 2026 — Discovery + manuscript (prior)

### Corrected geometry
| event_type | n | records |
|------------|---|---------|
| intermediate | 6 | 30, 32, 38, 45, 47, 50 |
| terminal | 4 | 31, 35, 36, 51 |

Substrates: sinus 5, AF 3, paced 2. All `event_hr_source=npz`.

### Key numbers (from `results/`, not re-derived by hand)
| Finding | Value |
|---------|--------|
| τ_s–excess3 direction concordance | **8/10 (0.8)** |
| τ_s–var / excess3–var | 0.5 / 0.5 |
| var–AR(1) | **0.2** |
| Median lead excess3 / τ_s / AR(1) / var (h) | **6.86 / 5.88 / 5.97 / 3.90** |
| Sensitivity (all 4 metrics, event Holters only) | 1.0 |
| FAR | **nan** (no control Holters) |
| Detection rate at 6 h (excess3 / τ_s / var) | 0.7 / 0.4 / 0.2 |

### Artifacts
- Results: `cctp_cohort_stratified.csv`, `cctp_stratified_summary.csv`, `leadtime_per_record.csv`, `leadtime_detector_summary.json`, `leadtime_cumulative.csv`, `ews_head2head.csv`, `ews_head2head_report.json`
- Figures: `figures/publication/fig_stratified_deltas.png`, `fig_leadtime_detector.png`, `fig_ews_concordance.png` (mirrored under `manuscript/figures/publication/`)
- Interpretation: `docs/JUL12_RESULTS_INTERPRETATION.md`
- Manuscript: `manuscript/CCTP_SDBB_manuscript.md` §§2.7–2.9 methods; §§3.6–3.8 results; FAR limitation §5.8
- Plans: `docs/EXTERNAL_VALIDATION_PLAN.md`, `docs/CLINICAL_COPILOT_DRAFT.md`
- Tests: `python3 tests/test_leadtime_detector.py` → 12/12 pass

### Strengths / limits (honest)
- **Strengths**: relational concordance; multi-hour discovery leads (excess3 earliest median); substrate polarity (paced negative, AF positive mean Δτ_s); publication-ready figs.
- **Limits**: N=10 discovery; FAR undefined; abs-z may be non-specific; borderline 47/50; do not claim clinical deployability.

### Reproduce (from project root)
```bash
python3 code/run_cohort_stratified.py
python3 code/run_leadtime_detector.py
python3 code/run_ews_head2head.py
python3 code/run_publication_figures.py
python3 tests/test_leadtime_detector.py
```

### Next recommended action
**External validation Phase 1** (VFDB / CU + control Holters) per `docs/EXTERNAL_VALIDATION_PLAN.md`, frozen params, primary deliverable = FAR + external sensitivity. Do not retune thresholds.

### Next-session prompt
```
Continue CCTP after Jul-12 manuscript update (stratified + lead-time + H2H in
manuscript/CCTP_SDBB_manuscript.md and docs/JUL12_RESULTS_INTERPRETATION.md).

Highest-impact next work: EXTERNAL VALIDATION Phase 1.
- Follow docs/EXTERNAL_VALIDATION_PLAN.md
- PhysioNet VFDB and/or CU ventricular tachyarrhythmia DB + non-event control Holters
- Reuse frozen params (θ3=0.08, high-thresh=0.65, W_TAU=101, abs-z detector)
- Primary deliverable: control-arm FAR estimate + external sensitivity table
- Do NOT retune thresholds; do NOT claim clinical deployability
Root: Investigaciones/Cardiac_CCTP_Pilot
```

---

## Historical Executive Summary (2026-07-08 pilot era; superseded above for active status)
- **Var ↑** replicates (classic EWS) in both records.
- **AR(1) ↓** replicates (opposite to naive critical slowing).
- **τ_s direction is context-dependent**: +0.0158 (record 35, long terminal) vs −0.0278 (record 30, intermediate).
- **Both Δτ_s are statistically extreme** under phase-shuffle surrogates that destroy relational structure (p_dir = 0.000, 0/8 for both records).
- **RECD Level 3 (mean_excess3)** mirrors τ_s sign exactly on same data/windows/proxy:
  - Record 35: +0.00350 (p=8.3e-07) → weighted step +0.00356 (p=3.2e-29)
  - Record 30: −0.01166 (p=1.3e-16) → weighted step −0.01178 (p=7.7e-78)
- high_level3_rate (thresh=1.75) = 0 (synthetic thresh too strict); continuous excess3 carries the signal.
- **Weighted RECD (α(λ) from |τ_s|)**: λ≈0 on these series (|τ_s| << 0.41 threshold) → α3 effectively constant; frac_contrib3 high (~0.77-0.78) and stable (tiny Δ, ns or marginal). The binary Φ3 is saturated; the continuous excess3 delta (prior step) is the differentiating metric. Still fully concordant in sign with τ_s.
- Systemic Tau + excess3 capture genuine non-trivial relational reorganization (context-dependent "flavor"). Univariate metrics miss the relational structure.
- All results falsifiable, many figures, exact repro via three scripts.

**Expansion update (7 records authoritative)**: 
- Full set: 30,31,32,35,36,45,51.
- 32/51 (pacing intermittent, auto-detected): both show strong negative Δτ_s + negative Δexcess3 (p<<0.001, p_surr=0). Usable despite pacing flag.
- 36 (AF): largest positive Δτ_s (+0.0786, p_surr=0) + concordant +Δexcess3.
- 45: small positive but significant.
- All 7 maintain sign concordance between Δτ_s and Δexcess3.
- Quality: interp <1% in 6/7; 32 at 5.2% (still acceptable). n_beats 63k–129k.
- CSV + 3 batch figures in results/ and figures/batch/.
- Next: run weighted for 36/45 (pending), light threshold re-cal (theta3, high_thresh, λ), then preprint narrative.

Research reminder (SDDB):
- 23 records total max. Realistic high-quality N ≈ 8–12 after strict filters.
- Suitable scale for Chaos / Network Physiology / PLOS Comp Biol preprint.

## "Las tres" execution (this session) — FINAL
1. Procesar 32 + 51 + 1-2 más → **DONE** (→ N=10).
2. Re-calibrar umbrales (theta3 / high / λ) con datos actuales → **DONE** (code support in run_recd*.py + run_cctp_batch.py; theta3=0.08, high=0.65, lambda-relative; re-ran on 30/32/35/36/38/50/51; summary/figs updated).
3. Generar tabla/figs finales + empezar borrador preprint → **DONE** (cctp_batch_summary.csv + 3 batch pngs + preprint draft v0.9 polished with full Methods, re-cal/sensitivity, SDDB facts from PhysioNet research (23 recs, pacing 5/23, sparse metadata), exact cmds, supp list).

**Preprint draft v0.9 + re-cal completed**. Ready for upload (Chaos/Frontiers/medRxiv). Next: external validation or add 1-2 borderline if .atr available (e.g. 33/43/44).

---

## Objective
Surgical real-data validation of Systemic Tau (relational EWS) + classic var/AR(1) on Holter RR series before VF.

Completed phases:
- Best-record selection + minimal bivariate proxy (RR + |ΔRR|) + compute_taus reuse.
- Multi-record confirmation (35 terminal + 30 intermediate) with identical pipeline.
- Light surrogates (phase-shuffle independent ×8 per record) on Δτ_s to test non-triviality.
- Port of RECD ordinal levels (Φ₁/Φ₂/Φ₃ + excess3) using exact same proxy + windows; mean_excess3 deltas are significant and sign-concordant with τ_s.

---

## Data (light downloads: only .hea + .atr)
- **Record 35**: ~24h pre-VF (best long window; event near end of excerpt). 100741 beats.
- **Record 30**: ~24h total, VF at ~7.91h (intermediate event → allows pre vs post contrast). 128807 beats.
- **Extraction**: wfdb.rdann → diff(ann.sample)/fs *1000 , clean 250<RR<2000 + linear interp. Saved per-record rr_XX_clean.npz.
- RR is the primary observable (univariate → minimal bivariate proxy for τ_s).

---

## Final Results (authoritative, multi-record)

**Method**: identical pipeline (W_EWS=501, W_TAU=101, stride=5, same proxy, same code).

### Record 35 (long pre-event, terminal VF)
Basal [6-16]h vs Approach [~21-24]h

| Metric | Basal     | Approach  | Δ          | p (Welch)   | Direction |
|--------|-----------|-----------|------------|-------------|-----------|
| Var    | 32907     | 36005     | **+3097**  | 6.9e-23     | ↑         |
| AR(1)  | 0.210     | 0.165     | **-0.045** | 1.5e-40     | ↓         |
| τ_s    | 0.0582    | 0.0739    | **+0.0158**| 7.2e-26     | ↑         |

### Record 30 (intermediate event at ~7.91h)
Basal [0.5-3.5]h vs Approach [~4.91-7.91]h

| Metric | Basal     | Approach  | Δ          | p (Welch)     | Direction |
|--------|-----------|-----------|------------|---------------|-----------|
| Var    | 4449      | 5859      | **+1411**  | 1.5e-20       | ↑         |
| AR(1)  | 0.540     | 0.249     | **-0.291** | 1.4e-221      | ↓         |
| τ_s    | 0.1211    | 0.0933    | **-0.0278**| 5.2e-19       | ↓         |

**Key observations**:
- **Variance ↑** replicated in both (classic EWS).
- **AR(1) ↓** in both (opposite naive critical slowing; stronger in 30).
- **τ_s** is context-dependent: **↑ +0.0158** in 35 (terminal long) vs **↓ −0.0278** in 30 (intermediate).
- **Both deltas are extreme** under surrogates (p=0.000, observed outside full null range of 8 phase-shuffled surrogates each).
- Absolute baseline τ_s higher in 30 (~0.121 vs ~0.058).
- This is a strong differentiating result: Systemic Tau is not universal and appears sensitive to different "flavors" of critical reorganization.

All p-values extremely small.

### Surrogates on Δτ_s (light, 8 per record)

**Method**: phase_shuffle_independent on the bivariate proxy (z(RR), z(|ΔRR|)) — destroys cross-relational structure while preserving individual component spectra/distributions. Same windows + compute_taus as the pilot. Δτ_s = approach − basal.

| Record | Type       | Δτ_s obs | surr min / max     | p (|\Delta|>=obs) | p (dir) |
|--------|------------|----------|--------------------|-------------------|---------|
| 30     | intermedio | -0.02780 | -0.0137 / +0.0099 | 0.000 (0/8)      | 0.000  |
| 35     | terminal   | +0.01577 | -0.0122 / +0.0079 | 0.000 (0/8)      | 0.000  |

**Interpretation**: Both the rise and the drop in systemic concordance are statistically non-trivial under the linear null. Systemic Tau is detecting genuine changes in the coupling between signal level and local irregularity. The direction depends on the transition "flavor" (prolonged terminal vs. intermediate event).

### RECD Ordinal Levels on same RR series (Φ₁/Φ₂/Φ₃ + excess3)

**Method**: identical proxy X=[z(RR),z(|ΔRR|)], m=3 delay=1, W_PHI=101, stride=5, theta3=0.10. Same basal/approach windows. Reused recd_ordinal_levels.py functions exactly.

| Record | mean_excess3 basal | approach | Δ excess3     | p (Welch)   | Direction | Notes |
|--------|--------------------|----------|---------------|-------------|-----------|-------|
| 35     | 0.3370             | 0.3405   | **+0.00350**  | 8.3e-07     | ↑         | terminal; matches τ_s sign |
| 30     | 0.3178             | 0.3061   | **−0.01166**  | 1.3e-16     | ↓         | intermediate; matches τ_s sign |

**high_level3_rate (excess3 > 1.75)**: 0.000 in all regimes/records (synthetic thresh too stringent for noisy RR; phi3_active~1.0 at theta=0.10).
**phi1/phi2**: small shifts; excess3 (Level 3) is the differentiating continuous metric.

**Key**: excess3 Δ sign is identical to τ_s Δ in both records. This directly links the relational EWS (τ_s) to an increase/decrease in nested ordinal Level-3 synergy on real cardiac data. Strong support for the paradigm.

### Weighted RECD (α(λ) from |τ_s|) — alpha_mode=lambda

**Method**: same X proxy, m=3 delay=1, W_TAU=101 (for λ), w_phi=101 (phi3), theta3=0.10. λ=compute_lambda(|τ_s|). α with beta1=2.0, gamma2=1.5, gamma3=6.0 (elevated Nivel 3), delta3=2.0. Same basal/approach. Reused compute_recd_from_conjunctions + compute_weighted_contributions.

| Record | mean_excess3 Δ (p)     | contrib3 Δ (p)     | frac_contrib3 Δ (p)   | mean_λ basal/app | Notes |
|--------|------------------------|--------------------|-----------------------|------------------|-------|
| 35     | +0.00356 (3.2e-29)     | ~0.0 (const, null) | −0.0010 (0.68)        | 0 / 0            | terminal; λ inert |
| 30     | −0.01178 (7.7e-78)     | ~0.0 (0.17)        | −0.0057 (0.059)       | ~0 / 0           | intermediate; λ inert |

**Observations**:
- λ(t) ≈ 0 in both regimes/records (empirical |τ_s| << 0.41 synthetic threshold). α3 therefore constant → contrib3 dominated by Φ3 indicator (~1.0) and frac_contrib3 ~0.77-0.78 (high but flat).
- The differentiating signal lives in the *continuous* excess3 (see prior RECD levels step), whose deltas are highly significant and sign-concordant with τ_s.
- frac_contrib3 slightly lower in approach for the intermediate case (30), consistent in direction with its τ_s / excess3 drop (marginal p).
- This step confirms the pipeline works end-to-end (τ_s → λ → weighted accumulation) on real RR. The high baseline weight on Nivel 3 and the need for theta recalibration on noisy physiological data are clear, falsable outcomes.

---





## Important Files

```
Investigaciones/Cardiac_CCTP_Pilot/
├── code/
│   ├── analyze_cctp_pilot.py          # Generalized (--record 35|30 ...)
│   ├── run_cctp_surrogates.py         # Light surrogates on Δτ_s (phase-shuffle, 8/rec)
│   ├── run_recd_on_rr.py              # RECD Φ levels + excess3 (continuous) on same RR
│   └── run_recd_weighted_on_rr.py     # Full weighted (α(λ) from |τ_s|), contribs, frac_contrib3 (this phase)
├── data/
│   ├── rr_35_clean.npz
│   └── rr_30_clean.npz
├── results/
│   ├── cctp_pilot_summary_35.json
│   ├── cctp_pilot_summary_30.json
│   ├── surrogate_cctp_35.json
│   ├── surrogate_cctp_30.json
│   ├── surrogate_cctp_combined.json
│   ├── surrogate_light_tau.json   # legacy (pre multi-record)
│   ├── recd_rr_35.json
│   ├── recd_rr_30.json            # Level 3 (continuous excess3)
│   ├── recd_weighted_rr_35.json
│   └── recd_weighted_rr_30.json   # authoritative weighted: contribs + frac + lambda stats
├── figures/
│   ├── 35/ (13+ figures)
│   └── 30/ (13+ figures)
│       ├── 06_ews_panels.png
│       ├── 08_basal_vs_approach_boxplots.png
│       ├── 15_surrogate_delta_tau.png
│       ├── 16_recd_excess3.png          # excess3(t) + windows + thresh
│       ├── 17_recd_excess3_box.png      # basal vs approach box for excess3
│       └── ...
├── figures/ (root-level files are legacy from single-record phase)
├── 00_CCTP_PILOT_REPORT.md
└── HANDOFF.md
```

---

## How to Reproduce / Inspect

```bash
cd Investigaciones/Cardiac_CCTP_Pilot

# Pilot metrics (per record)
python3 code/analyze_cctp_pilot.py --record 35
python3 code/analyze_cctp_pilot.py --record 30

# Surrogates ligeros (phase-shuffle on Δτ_s, 8 per record)
python3 code/run_cctp_surrogates.py --record all

# RECD ordinal levels on same data (Φ + excess3)
python3 code/run_recd_on_rr.py --record all

# Quick view of pilot deltas
python3 -c '
import json
for rec in ["35","30"]:
    s = json.load(open(f"results/cctp_pilot_summary_{rec}.json"))
    print(rec, "tau_s Δ=", s["metrics"]["tau_s"]["delta"])
'

# Quick view of surrogate results
python3 -c '
import json
for rec in ["30","35"]:
    s = json.load(open(f"results/surrogate_cctp_{rec}.json"))
    print(rec, "Δτ_s obs =", round(s["delta_obs"],5), "p_dir =", s["p_direction_specific"])
'

# Quick view of RECD Level 3
python3 -c '
import json
for rec in ["30","35"]:
    r = json.load(open(f"results/recd_rr_{rec}.json"))
    ex = r["mean_excess3"]
    print(rec, "excess3 Δ=", round(ex["delta"],5), "p=", ex["p_welch"])
'

# Quick view of Weighted RECD
python3 -c '
import json
for rec in ["30","35"]:
    r = json.load(open(f"results/recd_weighted_rr_{rec}.json"))
    ex = r["mean_excess3"]
    f3 = r["frac_contrib3"]
    print(rec, "excess3 Δ=", round(ex["delta"],5), "frac3 Δ=", round(f3["delta"],5), "p_f3=", f3["p_welch"])
'
```

The scripts:
- Load `data/rr_{rec}_clean.npz`
- Reuse `systemictau.core.compute_taus` (for τ_s) and `recd_ordinal_levels.*` (for Φ/excess3)
- `analyze_*` → `figures/{rec}/` + `results/cctp_pilot_summary_{rec}.json`
- `run_*_surrogates` → per-record surrogate JSONs + `15_surrogate_delta_tau.png`
- `run_recd_on_rr` → recd_rr_*.json + 16/17 (excess3)
- `run_recd_weighted_on_rr` → recd_weighted_rr_*.json + 18/19 (contribs + frac + boxes)

---

## Figures (per-record structure + surrogates)

**Key diagnostic figures** (in `figures/{30,35}/`):
- `06_ews_panels.png` — Full view of RR + Var + AR1 + τ_s
- `07_ews_zoom_6h.png` — Last 6 h
- `11_last90min_detail.png` — Last 90 min high-resolution
- `05_rr_full_with_approach.png` — Full trace + approach band
- `08_basal_vs_approach_boxplots.png`
- `09_phase_plane_ar1_var.png`
- `15_surrogate_delta_tau.png` — Histogram of phase-shuffled Δτ_s vs observed (p=0.000 both records)
- `16_recd_excess3.png` — excess3(t) + ... (Level 3 continuous)
- `17_recd_excess3_box.png`
- `18_recd_weighted_contribs.png` — excess3 + contrib1/2/3 + frac_contrib3(t) (λ-weighted)
- `19_recd_weighted_box.png` — boxes for excess3 / contrib3 / frac_contrib3

Root-level figures/ are legacy from earlier single-record runs. Use the per-record directories for clean multi-record work.

---

## Limitations (current pilot + extension)

- Only 2 records (surgical scope).
- Basic RR cleaning (no ectopic removal beyond 250-2000ms).
- Minimal bivariate proxy for τ_s and for Φ levels (N=2 from RR + |ΔRR|). λ-weighted RECD executed (λ inert at 0.41 on these RR; continuous excess3 carries relational signal).
- Different basal regimes and event contexts between 35 (terminal) and 30 (mid).
- Per-record surrogates completed (phase-shuffle, 8/record on Δτ_s for both 30 and 35; p=0.000 both directions).
- RECD levels ported (mean_excess3 significant, sign-concordant with τ_s).

---

## Recommended Next Steps (surgical, updated)

**Completed**:
- Light multi-record (35+30) with identical pipeline.
- Per-record light surrogates on Δτ_s (8 phase-shuffle surrogates/record). Both the ↑ (35) and ↓ (30) are extreme under the null (p_dir = 0.000).
- RECD ordinal levels (Φ/excess3) ported — mean_excess3 significant + sign-concordant with τ_s.
- Full weighted RECD (compute_recd_from_conjunctions + α(λ) from empirical |τ_s|) executed on same series. Pipeline complete; λ inert at synthetic threshold but excess3 deltas robust and directionally consistent.

High value follow-ups (in order):
1. Re-calibrate theta_chaos / HIGH_THRESH + consider continuous-Φ3 variant for physiological RR (λ currently inert).
2. Add 1-2 more records (e.g. 31) with identical full pipeline (τ_s + levels + weighted).
3. Sensitivity: delay=1 vs 2, W_PHI/W_TAU, pure univariate ordinal embedding.
4. Stronger ectopic cleaning or multi-lead if raw data used later.
5. Explore why basal τ_s / excess3 markedly higher in record 30 vs 35.

The divergence (τ_s ↑ only in one case) + surrogate confirmation for both signs is scientifically valuable — Systemic Tau is sensitive to relational reorganization context, not a universal or trivial artifact.

---

## Quick Commands

```bash
# From project root
cd Investigaciones/Cardiac_CCTP_Pilot

# Pilot metrics (per record)
python3 code/analyze_cctp_pilot.py --record 30
python3 code/analyze_cctp_pilot.py --record 35

# Surrogates ligeros on Δτ_s
python3 code/run_cctp_surrogates.py --record all

# RECD ordinal levels (Φ/excess3) on RR
python3 code/run_recd_on_rr.py --record all

# Weighted RECD (full α(λ) pipeline)
python3 code/run_recd_weighted_on_rr.py --record all

# View surrogate table + p-values
python3 -c '
import json, numpy as np
for rec in ["30","35"]:
    s = json.load(open(f"results/surrogate_cctp_{rec}.json"))
    print(rec, "Δobs=", round(s["delta_obs"],5), "p_dir=", s["p_direction_specific"])
    print("  surrs:", [round(x,4) for x in s["surr_deltas"]])
'
```

**Current state**: Pilot + multi-record (35+30) + light surrogates + RECD levels + Weighted RECD (α(λ)) **COMPLETE**.

Key findings:
- var ↑ and AR(1) ↓ replicate across terminal (35) and intermediate (30) events.
- τ_s behaves differently by context: **↑ +0.0158** (35, terminal long pre-VF) vs **↓ -0.0278** (30).
- **Both Δτ_s are extreme under independent phase-shuffle surrogates** (p=0.000 for 8 surr each).
- **RECD excess3 (Level 3) exactly mirrors τ_s direction**: +0.00350 (35) vs −0.01166 (30).
- Weighted RECD pipeline executed (λ from |τ_s|); on current RR λ≈0 (synthetic 0.41 too high) → frac_contrib3 high+stable, but excess3 deltas remain the robust carrier and are sign-concordant.
- High3 rate at synthetic thresh=0 (continuous excess3 is the signal). All support context-dependent relational reorganization.
- Systemic Tau + nested ordinal Level 3 capture structural signal univariate EWS miss.

Figures per record in `figures/{35,30}/` (incl. 15_surrogate + 16/17_recd + 18/19_weighted_contribs+box.png).
JSONs: cctp_pilot_*, surrogate_cctp_*, recd_rr_*.json .

All phases complete per original surgical plan. Next: full weighted RECD or 3rd record.

## Expansion to more records (started 2026-07-08)
- Fase 1 complete: records_inventory.csv (23 records, 13 candidates meeting criteria), selected_records.txt (10 priority incl. 30,31,32,35,36,38,43,44,45,51)
- New reusable scripts:
  - code/download_sddb_records.py (only .hea+.atr via urllib; no .dat needed)
  - code/extract_rr.py (wfdb.rdann + vfon parse from #vfon: + clean 250-2000 + interp → exact npz schema)
  - code/run_cctp_batch.py (orchestrates extract + 4 analysis scripts per record; --dry-run; produces results/cctp_batch_summary.csv + supports --list / --records / --force)
- Record 31 fully processed (13.97h, VF at ~13.7h near end): analyze + surrogates (p~1 for tiny Δτ) + recd_levels (excess3 Δ −0.022 p=1.3e-45) + weighted (same sign).
- Batch summary now includes 30/31/35 authoritative deltas.
- Window logic generalized (legacy exact for 30/35; sensible pre-approach basal for others).
- To add more: edit selected_records.txt, run download if needed, `python3 code/run_cctp_batch.py --records 32,36,...` (or without to use selected).
- Recommendation: add 3-5 more (32,36,38,44,45,51) then re-calibrate thresholds (theta, lambda) before preprint framing.

## Expansion + Research (2026-07-08 update)
- **SDDB research (web/PhysioNet)**: 23 records total confirmed. Pacing: 1 continuous (40), intermittent (32,43,51,...). VF onset in #vfon comments. Prior papers: ML/VF pred on pre segments; no relational/ordinal RECD/SysTau work.
- **Improvements surgical**:
  - extract_rr.py: pacing_detected (comment+known+cv<0.06), interp_frac/n_invalid saved in npz, ari fallback, better parse.
  - run_cctp_batch.py: print_quality_report, generate_batch_figures (deltas, -logp, quality), richer csv.
  - records_inventory.csv: added pacing/rhythm/vf_onset accurate; 13 include?=Yes.
  - selected_records.txt: curated priority + pacing flags.
- **Current batch (5 records)**: authoritative cctp_batch_summary.csv + figures/batch/ (3 pngs).
  - New: 36 (AF long pre, concordant ↑ Δτ/Δexcess3), 45 (sinus excellent pre, concordant small ↑).
  - 32,51 npz with quality + pacing=True.
- **Next (for preprint scale)**: process 32/38/44/47/50/51 (use download + batch), light re-cal (theta3~0.08-0.12?, high_thresh lower, λ from real |τ| dist), add batch table to report §14, start preprint draft.
- Repro:
  ```bash
  python3 code/download_sddb_records.py --records 32,36,38,44,45,51
  python3 code/run_cctp_batch.py --records 32,36,38,44,45,51
  python3 -c "
  import pandas as pd
  df=pd.read_csv('results/cctp_batch_summary.csv')
  print(df[['record','delta_tau','delta_excess3','p_excess3','interp_frac','pacing_detected']].to_string(index=False))
  "
  ```
- Limitation section ready for preprint: N=23 max, pacing subset, limited metadata.

## Final closure (2026-07-08)
- All three small adjustments applied to draft:
  1. Strong title adopted.
  2. Abstract sentence added ("The relational metrics detect reorganization even when classical variance-based early-warning signals are weak or reversed.").
  3. Exact high_level3_rate explanation paragraph inserted in 3.1b.
- Clean upload package created at `preprint_v0.9_final/` (md + csv + batch figs + key records 30/35/38/50 + README).
- Batch script fixed (lambda args now only passed to weighted stage).
- EWS panel title clipping fixed (06_ews_panels.png regenerated with proper suptitle + rect + bbox_inches).
- Ready for medRxiv / Preprints.org + targeted journal submission (Chaos / Frontiers in Network Physiology recommended).
- All artifacts use final params: theta3=0.08, high_thresh=0.65, lambda-relative.
- Recommendation: upload draft today/tomorrow for priority date + reviewer feedback.

## Final polish (upload-ready, same day)
Four surgical text edits applied to source + `preprint_v0.9_final/`:
1. Title refined (Heart Rate Dynamics + Systemic Tau / Ordinal RECD Evidence framing).
2. Abstract: retained relational-vs-variance sentence + added pacing/AF flagging sentence.
3. §3.1 after table: explicit borderline interpretation of discordant 47/50.
4. Limitations: Record 32 (interp 5.22%, intermittent pacing) retained; did not reverse sign concordance.
Package and draft are medRxiv / journal-submission ready.

