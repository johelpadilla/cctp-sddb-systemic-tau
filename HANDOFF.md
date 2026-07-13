# CCTP — Handoff (authoritative status)

**Date**: 2026-07-13 (v1.3.0 — ordinal package + manuscript §3.11)  
**Root**: `/Users/johelpadilla/grok-safe/Investigaciones/Cardiac_CCTP_Pilot/`  
**GitHub**: https://github.com/johelpadilla/cctp-sddb-systemic-tau  
**Release**: **v1.3.0** — https://github.com/johelpadilla/cctp-sddb-systemic-tau/releases/tag/v1.3.0  
**Zenodo DOI (v1.3.0)**: https://doi.org/10.5281/zenodo.21344730  
**Zenodo concept DOI**: https://doi.org/10.5281/zenodo.21270698  
**Prior version DOI (v1.2.0)**: https://doi.org/10.5281/zenodo.21327196

**Status snapshot**

| Workstream | State |
|------------|--------|
| SDDB discovery N=10 (stratified + lead-time + H2H) | **Done** (v1.2.0+) |
| Manuscript + PDF (Phase 1 §3.9 + Phase 2 §3.10 + **§3.11 ordinal**) | **Done** (v1.3.0) |
| External Validation Phase 1 (VFDB + NSRDB FAR) | **Done & reported** |
| External Validation Phase 2 public interim (NSRDB n=18) | **Done & reported** |
| Phase 2 plan / Tier A path | **Prepared** — see `docs/EXTERNAL_VALIDATION_PHASE2_PLAN.md` |
| Tier A institutional data-request draft | **Prepared** (no partner data yet) |
| GitHub push + Zenodo archive | **v1.3.0** (this release) |
| Native ordinal formalization (OPC + SDD) | **Done & shipped** (v1.3.0) |
| Exploratory bake-off SDDB + VFDB | **Done & shipped** |
| Comparable NSRDB FAR (OPC L=50 / SDD / abs-z) | **Done & shipped** |
| Sensitivity × FAR trade-off (join-only) | **Done & shipped** |
| Light cascade SDD→OPC (±5 min, **causal**) | **Done & shipped** — exploratory; **low_priority** |
| **OPC modest param grid** (L, θ_D, θ_R; 36 cells) | **Done & shipped** — **keep_baseline** |
| Clinical Copilot | Draft only — **do not deploy** |

**Frozen production params** (never retune on validation without pre-registration):  
θ₃=0.08, high-threshold=0.65, W_τ=101, W_EWS=501, stride=5, relative λ, detector **abs-z≥2 × 3** consecutive windows, RR clean [250, 2000] ms.

**Clinical / FDA / deployability claim**: **NONE**. S5 (FAR ≤ 2/24h) is **not claimed**.

---

## Just completed (recent sessions)

### 1. Formalize two separate ordinal detectors (no fusion in singleton APIs)

| Option | Name | Core rule |
|--------|------|-----------|
| **1** | Ordinal Persistence Collapse (**OPC**) | \(D_t \le \theta_D\) **and** \(R_t \ge \theta_R\) |
| **2** | Symbolic Distribution Divergence (**SDD**) | \(\mathrm{TV}(P_t, P_{\mathrm{basal}}) \ge \theta_{\mathrm{TV}}\) (KL secondary only) |

- Spec: `docs/ORDINAL_ALARM_DETECTORS.md`
- Code: `code/ordinal_detectors/{opc,sdd}_detector.py` (pure symbol stream; **no** mean/var of continuous signal)
- Tests: `tests/test_ordinal_detectors.py`
- Singleton modules still declare “Does NOT fuse”; cascade is a **separate post-process** module

### 2. Exploratory bake-off on SDDB + VFDB

- Runner: `code/run_ordinal_exploratory_bakeoff.py`
- Write-up: `docs/ORDINAL_EXPLORATORY_BAKEOFF.md`
- Joint bivariate symbols **K=36** (Bandt–Pompe m=3, σ = s0×6+s1) from cleaned RR
- OPC **L=8 saturates** with K=36; companion **OPC L=50** is the interpretable arm

### 3. Comparable NSRDB FAR (Phase-2 episode definition)

- Runner: `code/run_ordinal_far_comparison.py` — **`fusion=False`** (singleton arms only)
- Write-up: `docs/ORDINAL_NSRDB_FAR_COMPARISON.md`
- Helper: `count_binary_alarm_episodes` in `code/cctp_metrics_core.py`
- Windowing = Phase 2: basal ~`(0.25, min(2, 0.25·H))`, search remainder, **max 12 h**, refractory **0.5 h**
- FAR = total_episodes / total_search_hours × 24

### 4. Combined sensitivity × FAR trade-off (join-only)

- Runner: `code/run_ordinal_tradeoff_analysis.py` (existing JSON/CSV only)
- Write-up: `docs/ORDINAL_SENSITIVITY_SPECIFICITY_TRADEOFF.md`
- Primary arms: **OPC L=50**, **SDD**, **abs-z frozen** (OPC L=8 excluded)

### 5. Light cascade SDD→OPC (±5 min) — causal; low_priority

**Rule (fixed params, not retuned):**
1. **Candidate:** SDD (`L_c=50`, `θ_TV=0.35`, `θ_S=1`)
2. **Confirm:** OPC L=50 (`θ_D=0.35`, `θ_R=5`) within closed window **[t−W, t+W]**, **W = ±5 min**
3. **Causality (events):** only OPC with `t_opc < event_hr`;  
   `decision_time = max(t_SDD, t_OPC_confirm) < event_hr`  
   (**no post-event look-ahead**)

| Artifact | Path |
|----------|------|
| Pure merger | `code/ordinal_detectors/cascade_fusion.py` |
| Evaluation entry | `code/run_ordinal_cascade_fusion.py` |
| Tests | `tests/test_ordinal_cascade_fusion.py` (**20 passed**) |
| Write-up | `docs/ORDINAL_CASCADE_FUSION.md` |
| Results | `results/ordinal_cascade_{comparison,per_record,nsrdb_far_per_record,gain_loss,summary}.*` |

**Head-to-head (causal cascade vs singletons):**

| Detector | Sens SDDB (11) | Sens VFDB (22) | Sens all (33) | FAR /24h (NSRDB 18) | Episodes |
|----------|----------------|----------------|---------------|---------------------|----------|
| **Cascade SDD→OPC (causal)** | 0.545 (6/11) | 0.318 (7/22) | **0.394** (13/33) | **~3.87** | 29 |
| OPC L=50 alone | 0.545 (6/11) | 0.364 (8/22) | **0.424** (14/33) | **~3.73** | 28 |
| SDD alone | 1.000 (11/11) | 0.955 (21/22) | **0.970** (32/33) | **~46.3** | 347 |
| abs-z τ_s frozen | 1.000 (11/11) | 0.864 (19/22) | **0.909** (30/33) | **~33.73** | 253 |

**Recommendation:** `low_priority` — cascade does **not** beat OPC L=50 on both sens and FAR. Do **not** promote. Do **not** fish windows post-hoc.

### 6. Modest OPC (L, θ_D, θ_R) parameter exploration — **just finished**

**Goal:** See if modest, justified knobs improve OPC sensitivity without wrecking NSRDB FAR. Exploratory only; **not** clinical optimization; **abs-z not retuned**.

**Grid (36 cells):**

| Axis | Values | Justification |
|------|--------|---------------|
| L | {40, 50, 60, 70} | Around L=50 companion for K=36 (need L > θ_D·K ≈ 12.6) |
| θ_D | {0.30, 0.35, 0.40} | Around 0.35 diversity threshold |
| θ_R | {4, 5, 6} | Around 5 consecutive low-div windows |

**Methodology (same as bake-off + FAR):**
- Event sens: first OPC alarm post-basal / pre-event (SDDB n=11, VFDB n=22)
- FAR: binary-alarm episodes, refractory 0.5 h, Phase-2 basal/search, cap 12 h (NSRDB n=18)
- Symbols built **once per record**, then OPC re-run across all cells
- Slack for “FAR not much worse”: ≤ max(1.5× baseline FAR, baseline+2) ≈ **5.73 /24h**

| Artifact | Path |
|----------|------|
| Runner | `code/run_ordinal_opc_param_explore.py` |
| Tests | `tests/test_ordinal_opc_param_explore.py` (**16 passed**) |
| Grid table | `results/ordinal_opc_param_explore_grid.csv` |
| Summary JSON | `results/ordinal_opc_param_explore_summary.json` |
| Report | `results/ordinal_opc_param_explore_report.md` |

**Baseline row (reproduced; matches prior bake-off/FAR):**

| Metric | Value |
|--------|-------|
| L, θ_D, θ_R | 50, 0.35, 5 |
| sens_sddb | 6/11 ≈ 0.545 |
| sens_vfdb | 8/22 ≈ 0.364 |
| **sens_all** | **14/33 ≈ 0.424** |
| **FAR** | **3.733 /24h** (28 ep) |

Dual full-grid runs: identical baseline counts (reproducible).

**Recommendation: `keep_baseline`**

No cell raised all-event sensitivity with FAR within exploratory slack. Sensitivity gains only with substantially higher FAR:

| L | θ_D | θ_R | sens_all | FAR /24h | Class |
|---|-----|-----|----------|----------|--------|
| **50** | **0.35** | **5** | **0.424** | **3.73** | **baseline** |
| 50 | 0.35 | 4 | 0.424 | 3.87 | same_sens_far_ok |
| 50 | 0.35 | 6 | 0.394 | 3.33 | lower_sens_lower_far |
| 40 | 0.35 | 6 | 0.545 | 9.73 | higher_sens_far_worse |
| 60 | 0.40 | 4 | 0.515 | 8.00 | higher_sens_far_worse |
| 50 | 0.40 | 6 | 0.636 | 14.5 | higher_sens_far_worse |
| 40 | 0.40 | 6 | 0.788 | 28.1 | higher_sens_far_worse |

**Qualitative effects (grid margins, descriptive only):**
- **L↑** → lower mean sens and lower FAR (longer window damps brief collapses)
- **θ_D↑** → main lever: higher sens **and** much higher FAR
- **θ_R↑** → small effect (Δmean_sens ≈ −0.01; mild FAR drop)

**Do not adopt new OPC defaults** on this grid. Keep **(L=50, θ_D=0.35, θ_R=5)**. Do not expand to exhaustive search without a priori design + hold-out; n is small.

**Reproduce:**
```bash
python3 code/run_ordinal_opc_param_explore.py --write-report   # ~3–4 min full grid
python3 -m pytest tests/test_ordinal_opc_param_explore.py -q
```

---

## Observational joint table (singleton trade-off; still valid)

| Detector | Sens SDDB (n=11) | Sens VFDB (n=22) | Sens all (n=33) | FAR /24h (NSRDB n=18) | Control episodes |
|----------|------------------|------------------|-----------------|------------------------|------------------|
| **OPC L=50** | 0.545 (6/11) | 0.364 (8/22) | **0.424** (14/33) | **~3.73** | 28 |
| **SDD** | 1.000 (11/11) | 0.955 (21/22) | **0.970** (32/33) | **~46.3** | 347 |
| **abs-z τ_s** | 1.000 (11/11) | 0.864 (19/22) | **0.909** (30/33) | **~33.73** | 253 |

**Objective-conditioned balance (observational only — not clinical ranking):**

| Objective | Preferred arm | Rationale (short) |
|-----------|---------------|-------------------|
| Lowest false-alarm burden | OPC L=50 | FAR ~9× lower than abs-z; accepts lower hit rate |
| Maximum event hit rate | SDD | 32/33 events; highest control FAR |
| Frozen pilot baseline continuity | abs-z τ_s | Production-frozen continuous detector |
| Sens per unit FAR | OPC L=50 | Highest sens_all / FAR among fixed-param arms |
| Light cascade SDD→OPC | **Not preferred** | Sens ≤ OPC, FAR ≈ OPC |
| OPC param retune (36-cell grid) | **Keep baseline** | No Pareto-ish gain within FAR slack |

---

## Shipped package (v1.3.0)

```
docs/ORDINAL_ALARM_DETECTORS.md
docs/ORDINAL_EXPLORATORY_BAKEOFF.md
docs/ORDINAL_NSRDB_FAR_COMPARISON.md
docs/ORDINAL_SENSITIVITY_SPECIFICITY_TRADEOFF.md
docs/ORDINAL_CASCADE_FUSION.md
code/ordinal_detectors/                 # opc + sdd + cascade_fusion (post-process only)
code/run_ordinal_exploratory_bakeoff.py
code/run_ordinal_far_comparison.py      # fusion=False guards intact
code/run_ordinal_tradeoff_analysis.py
code/run_ordinal_cascade_fusion.py      # exploratory cascade arm
code/run_ordinal_opc_param_explore.py   # modest L/θ_D/θ_R grid → keep_baseline
code/cctp_metrics_core.py               # + count_binary_alarm_episodes only (abs-z untouched)
tests/test_ordinal_detectors.py
tests/test_ordinal_exploratory_bakeoff.py
tests/test_ordinal_far_comparison.py
tests/test_ordinal_tradeoff_analysis.py
tests/test_ordinal_cascade_fusion.py
tests/test_ordinal_opc_param_explore.py
results/ordinal_*                       # bake-off + FAR + trade-off + cascade + opc_param_explore
HANDOFF.md
```

**Verify:**
```bash
cd /Users/johelpadilla/grok-safe/Investigaciones/Cardiac_CCTP_Pilot
python3 -m pytest tests/test_ordinal_detectors.py \
  tests/test_ordinal_exploratory_bakeoff.py \
  tests/test_ordinal_far_comparison.py \
  tests/test_ordinal_tradeoff_analysis.py \
  tests/test_ordinal_cascade_fusion.py \
  tests/test_ordinal_opc_param_explore.py \
  tests/test_far_and_short_windows.py -q
python3 code/run_ordinal_tradeoff_analysis.py              # join-only; instant
python3 code/run_ordinal_cascade_fusion.py --write-doc     # ~40–60 s; causal cascade
# optional re-run OPC grid (~3–4 min):
# python3 code/run_ordinal_opc_param_explore.py --write-report
# optional re-run FAR (~20–30 s):
# python3 code/run_ordinal_far_comparison.py --write-doc
```

---

## v1.2.0 public close-out (already on GitHub + Zenodo)

| Fact | Phase 1 | Phase 2 interim |
|------|---------|-----------------|
| NSRDB n | 6 | **18** |
| Search hours (capped 12 h/record) | ~60 | **~180** |
| τ_s FAR / 24 h | ~34.4 | **~33.7** |
| excess3 FAR / 24 h | ~28.8 | **~32.3** |
| Device-matched to VFDB? | **No** | **No** |
| S5 (FAR ≤2/24h) | Failed | **Still failed / not claimed** |

---

## Source-of-truth documents

| Artifact | Role |
|----------|------|
| `manuscript/CCTP_SDBB_manuscript.md` + `.pdf` | Public manuscript package (v1.3.0; includes §3.11) |
| `docs/EXTERNAL_VALIDATION_PHASE2_PLAN.md` | Phase 2 plan (Tier A priority; still active next step) |
| `docs/ORDINAL_ALARM_DETECTORS.md` | Formal OPC + SDD (singleton APIs) |
| `docs/ORDINAL_EXPLORATORY_BAKEOFF.md` | SDDB/VFDB exploratory sensitivity |
| `docs/ORDINAL_NSRDB_FAR_COMPARISON.md` | Comparable FAR methodology + tables |
| `docs/ORDINAL_SENSITIVITY_SPECIFICITY_TRADEOFF.md` | Joint sens × FAR trade-off (join-only) |
| `docs/ORDINAL_CASCADE_FUSION.md` | Causal SDD→OPC ±5 min exploratory cascade |
| `results/ordinal_opc_param_explore_report.md` | OPC L/θ_D/θ_R modest grid + **keep_baseline** |
| `docs/PHASE2_INSTITUTIONAL_DATA_REQUEST.md` | Tier A partner brief |
| `docs/PHASE2_IRB_DATA_CHECKLIST.md` | IRB / partner checklist |
| `results/external_phase2_far.json` | Phase 2 abs-z primary FAR |
| `results/ordinal_nsrdb_far_summary.json` | Ordinal vs abs-z FAR aggregates |
| `results/ordinal_tradeoff_summary.json` | Full trade-off narrative + rankings |
| `results/ordinal_cascade_summary.json` | Cascade head-to-head + recommendation |
| `results/ordinal_opc_param_explore_summary.json` | OPC param grid + recommendation |
| `code/cctp_metrics_core.py` → `detect_lead_time` / `count_alarm_episodes` | **Frozen** abs-z path |

### Reproduce Phase 1 + Phase 2 + ordinal package
```bash
cd /Users/johelpadilla/grok-safe/Investigaciones/Cardiac_CCTP_Pilot
python3 code/run_external_validation_phase1.py
python3 code/run_external_validation_phase2_far.py
python3 code/run_ordinal_exploratory_bakeoff.py
python3 code/run_ordinal_far_comparison.py --write-doc
python3 code/run_ordinal_tradeoff_analysis.py
python3 code/run_ordinal_cascade_fusion.py --write-doc
python3 code/run_ordinal_opc_param_explore.py --write-report
```

---

## Discovery N=10 (foundational)

| Finding | Value |
|---------|--------|
| event_type | intermediate 6 (30,32,38,45,47,50); terminal 4 (31,35,36,51) |
| τ_s–excess3 concordance | **8/10 (0.8)** |
| Median lead excess3 / τ_s / var (h) | **6.86 / 5.88 / 3.90** |

Artifacts: `results/cctp_*stratified*`, `leadtime_*`, `ews_head2head*`, `figures/publication/`.

---

## Claims boundary

| Claim | Status |
|-------|--------|
| Discovery relational pre-VF reorganization (SDDB N=10) | **Reported** |
| External sensitivity (Phase 1 VFDB) | **Reported** |
| Public interim control FAR abs-z (NSRDB n=18) | **Reported** (~33.7/24h) |
| Ordinal FAR comparison (OPC/SDD vs abs-z) | **Shipped (exploratory)** — observational |
| Ordinal sens × FAR trade-off | **Shipped (exploratory)** — exploratory join; no clinical winner |
| Causal cascade SDD→OPC ±5 min | **Shipped (exploratory)** — exploratory; **not preferred** |
| OPC L/θ_D/θ_R modest grid (36 cells) | **Shipped (exploratory)** — **keep_baseline**; no better balance within FAR slack |
| S5 met | **False / not claimed** |
| Clinical deployability / FDA | **None** |
| Device-matched institutional FAR | **Not yet available** |
| OPC/SDD clinical superiority vs abs-z | **Not claimed** |
| Production fused detector | **Not shipped** — cascade is exploratory post-process only |
| New OPC production defaults | **Not adopted** — keep L=50, θ_D=0.35, θ_R=5 |

---

## Next steps (pick one track; do not mix casually)

### Track A — Ship ordinal package (**done in v1.3.0**)
1. Review formalization + bake-off + FAR + trade-off + cascade + **OPC param explore** docs/results.
2. Commit + push full uncommitted ordinal set (all ORDINAL docs, detectors, runners, tests, results, HANDOFF).
3. Suggested message theme:  
   `ordinal: OPC/SDD formalization, bake-off, NSRDB FAR, sens×FAR trade-off, causal cascade exploratory, OPC L/θ_D/θ_R grid keep_baseline (no abs-z retune)`.
4. Keep abs-z as reported primary detector; cascade **low_priority**; OPC params **unchanged** at (50, 0.35, 5).

### Track B — Institutional Phase 2 (specificity that can support S5)
1. Identify **one** real partner P1–P5 (`docs/PHASE2_INSTITUTIONAL_DATA_REQUEST.md`).
2. Checklist §A–B with IRB/DUA.
3. Device-matched non-event controls (n≈10–20).
4. **Do not** invent inventory rows, retune frozen abs-z, or claim S5 from NSRDB alone.

### Track C — Ordinal research (only if deliberate)
- OPC L/θ_D/θ_R **already modest-grid explored** → **keep_baseline**; do **not** re-grid the same axes without a priori design + hold-out.
- Pre-register any **SDD** basal/threshold policy changes before re-running FAR (next open lever if any).
- Cascade SDD→OPC ±5 min **already evaluated (causal)** → **low_priority**; do **not** re-grid confirm window without a priori justification.
- Do **not** retune abs-z production constants in place.

### Explicit non-goals (still true)
- No promotion of cascade (or any fusion) to production without a new design decision + partner data
- No swap of production `detect_lead_time` / abs-z
- No adopt of looser OPC θ_D / smaller L as “defaults” from this grid (overfits small n)
- No S5 / clinical / FDA claims from public Holters
- No invented partners or PHI
- No post-hoc window fishing on cascade (±10 min optional robustness only if pre-registered)
- No reprocessing of RR solely to re-derive the join-only trade-off table

---

## Next-session prompts (copy-paste)

**Ship ordinal work (recommended):**
```
Commit and push the full ordinal package (formalization + bake-off + NSRDB FAR +
sens×FAR trade-off + causal cascade SDD→OPC exploratory + OPC L/θ_D/θ_R param
explore with keep_baseline).
Root: /Users/johelpadilla/grok-safe/Investigaciones/Cardiac_CCTP_Pilot
Include: docs/ORDINAL_*.md, code/ordinal_detectors/, code/run_ordinal_*.py,
  code/cctp_metrics_core.py (count_binary_alarm_episodes only),
  tests/test_ordinal_*.py, results/ordinal_*, HANDOFF.md
Do NOT retune frozen abs-z; no clinical/S5 superiority claims.
OPC defaults remain L=50, θ_D=0.35, θ_R=5 (grid found no better balance).
Cascade is exploratory/low_priority (causal confirm, sens 13/33, FAR ~3.87).
FAR runner remains fusion=False for singleton arms.
```

**Institutional Tier A:**
```
Continue External Validation Phase 2 — Tier A institutional path (post v1.2.0).
Base: public manuscript final + NSRDB FAR n=18 (τ_s ~33.7/24h, device mismatch).
Partner brief: docs/PHASE2_INSTITUTIONAL_DATA_REQUEST.md
Checklist: docs/PHASE2_IRB_DATA_CHECKLIST.md
FROZEN: θ₃=0.08, high=0.65, W_TAU=101, abs-z≥2×3 — do not retune.
Priority: identify one real P1–P5 partner; complete checklist §A–B; no invented PHI; no S5/clinical claim.
Root: /Users/johelpadilla/grok-safe/Investigaciones/Cardiac_CCTP_Pilot
```

**Ordinal research (only if deliberate; OPC grid + cascade already done):**
```
Do NOT re-open OPC L/θ_D/θ_R grid (36 cells → keep_baseline) or cascade ±5 min
(causal → low_priority) without a priori design + hold-out.
If continuing ordinal work: pre-register one SDD-only lever (basal / θ_TV / θ_S)
or wait for Tier A partner data. Do not retune abs-z.
Baseline: results/ordinal_opc_param_explore_report.md + docs/ORDINAL_CASCADE_FUSION.md
Root: /Users/johelpadilla/grok-safe/Investigaciones/Cardiac_CCTP_Pilot
```
