# CCTP — Handoff (authoritative status)

**Date**: 2026-07-14  
**Root**: `/Users/johelpadilla/grok-safe/Investigaciones/Cardiac_CCTP_Pilot/`  
**GitHub**: https://github.com/johelpadilla/cctp-sddb-systemic-tau  
**Last shipped release**: **v1.4.0** — https://github.com/johelpadilla/cctp-sddb-systemic-tau/releases/tag/v1.4.0  
**Zenodo DOI (v1.4.0)**: https://doi.org/10.5281/zenodo.21348295  
**Zenodo concept DOI**: https://doi.org/10.5281/zenodo.21270698  
**Prior version DOI (v1.3.0)**: https://doi.org/10.5281/zenodo.21344730  
**HEAD (local `main`)**: `563f0cf` — archive Surplus/P2/Surrogate diagnostic line (closed; no detector)  
**origin/main (shipped)**: `6f19efb` (v1.4.0 DOI pin + HANDOFF); release content `46b097b`  
**Local vs origin:** `main` **ahead 1** (not pushed).  

**Surplus line:** closed and archived on local `main`. Authoritative close: **`docs/SURPLUS_LINE_CLOSURE.md`**. Not a detector; not I0.

**Clinical / FDA / deployability claim**: **NONE**.

---

## Status snapshot

| Workstream | State |
|------------|--------|
| SDDB discovery + manuscript Phase 1/2 + §3.11 ranking | **Done** — **v1.4.0** shipped (GitHub + Zenodo) |
| External Validation Phase 1 / Phase 2 public FAR | **Done & reported** |
| Phase 2 plan / Tier A path | **Prepared** — main institutional next step |
| Native ordinal OPC/SDD + cascade + trade-off | **Done** (v1.3.0+) |
| Integrated OPSP + I0 Holter + θ-grid 0/16 | **Done** — do **not** promote I0 |
| Structural surplus arms R0–R5 | **Done** — 0/6; **stop surplus-primary** |
| Final detector ranking (abs-z primary) | **Done** — manuscript + `docs/FINAL_DETECTOR_RANKING.md` |
| **Surplus / P2 / Surrogate line (A→E)** | **Closed locally 2026-07-14** — see `docs/SURPLUS_LINE_CLOSURE.md` |
| Dual-mode C-BR / Mode-S multi-config Holter bake-off | Optional / low value after S1 stop |
| Clinical Copilot | Draft only — **do not deploy** |

**Frozen production params** (never retune on validation without pre-registration):  
θ₃=0.08, high-threshold=0.65, W_τ=101, W_EWS=501, stride=5, relative λ, detector **abs-z≥2 × 3** consecutive windows, RR clean [250, 2000] ms.

**OPC companion**: \(L=50\), \(\theta_D=0.35\), \(\theta_R=5\), joint K=36 — **keep_baseline**.

---

## Final detector ranking (locked — v1.4.0)

| Rank | Arm | Role |
|------|-----|------|
| **1** | **abs-z \(\tau_s\)** | Preferred **primary** for pre-VF event hit rate (~0.91 sens) |
| 2 | SDD | Sensitivity ceiling (~0.97); high FAR |
| 3 | I0 surplus | Mode-S narrative only; **not promoted** |
| 4 | OPC L=50 | Specificity companion (~3.73 FAR); low sens |
| 5 | I-confirm | Secondary FAR filter only |

See: `docs/FINAL_DETECTOR_RANKING.md`, manuscript §3.11 / §4.1 / §6.

---

## Closed tracks (do not reopen without new data / new statistic)

| Track | Outcome |
|-------|---------|
| I0 θ-grid (D′) | **0/16** clear advances |
| Structural surplus (S1) R0–R5 | **0/6** structural wins; **stop surplus-primary** |
| Promote I0 or OPC over abs-z (event hit rate) | **Rejected** |
| Another I0-like θ / micro-structural grid on current proxy | **Forbidden** — floor is structural |
| **Surplus proxy redesign + P2 + Phase β \(S^{\mathrm{ex}}\)** | **Closed 2026-07-14** — floor improved; no approach rise; **not interesting**; archive only |

---

## Local research (uncommitted — independent of detector ranking)

**Line status: CLOSED.** Full narrative close → **`docs/SURPLUS_LINE_CLOSURE.md`**.

Do **not** merge these into I0 / ranking docs. Optional separate commit/release only if archiving (Track D).

### A. Surplus NSRDB diagnostic (closed analysis)

**Purpose:** Why \(S_t\) is noisy (high FAR) on healthy NSRDB — **not** a detector.

| Artifact | Path |
|----------|------|
| Report | `docs/SURPLUS_NSRDB_DIAGNOSTIC.md` |
| Runner | `code/run_surplus_nsrdb_diagnostic.py` |
| Results | `results/surplus_diag_*.csv` + `surplus_diag_summary.json` |

| Contrast | Value |
|----------|-------|
| NSRDB mean of mean \(S\) (18 × 12 h) | **≈ 0.373** |
| NSRDB mean of p90 \(S\) | **≈ 0.441** |
| Fraction windows \(S \ge 0.40\) | **≈ 30%** |
| Event basal mean \(S\) (32 records) | **≈ 0.448** (higher than NSRDB; \(d \approx -1.83\)) |
| Approach − basal \(\Delta S\) | **≈ 0** (~41% records rise) |
| Dominant co-factor of high \(S\) | **MI_bits** (\(r \approx 0.96\)) |

**Root cause:** Proxy \([z(\mathrm{RR}),\, z(|\Delta\mathrm{RR}|)]\) is structurally non-independent in healthy HRV; \(S_t\) tracks ordinal MI, not a rare pre-VF synergy state.

```bash
cd /path/to/cctp-sddb-systemic-tau
PYTHONPATH=code python3 code/run_surplus_nsrdb_diagnostic.py
# smoke: ... --smoke
```

### B. Proxy redesign (design only — closed)

**Purpose:** Propose alternative bivariate proxies so healthy null is nearer independence, without optimizing I0 gates.

| Artifact | Path |
|----------|------|
| Design note | `docs/SURPLUS_PROXY_REDESIGN.md` |

### C. Phase-α proxy bake-off (closed diagnostic)

**Purpose:** Compare five frozen proxies on NSRDB floor + event basal/approach contrast. **Not** a detector; **not** I0.

| Artifact | Path |
|----------|------|
| Report | `docs/SURPLUS_PROXY_BAKEOFF.md` |
| Proxies | `code/proxy_variants.py` (does **not** mutate `build_bivariate_proxy`) |
| Runner | `code/run_surplus_proxy_bakeoff.py` |
| Results | `results/surplus_bakeoff_*.csv` + `surplus_bakeoff_summary.json` |

| Proxy | NSRDB mean \(S\) | frac \(S\ge0.40\) | Event basal mean \(S\) | \(d\) (event vs NSRDB) | Promising? |
|-------|------------------|-------------------|------------------------|------------------------|------------|
| R0 | 0.373 | 0.301 | 0.448 | 1.83 | baseline |
| P1 AR(1) residual | **0.524** (worse) | 0.963 | 0.619 | 2.47 | **No** |
| **P2 lag-10** | **0.269** (−28%) | **0.005** (−98%) | 0.289 | 0.78 | **Yes** |
| P3 signed Δ | 0.492 (worse) | 0.897 | 0.563 | 2.64 | **No** |
| P4 \|Δ\| residual | 0.386 | 0.395 | 0.447 | 1.36 | **No** |

**Key takeaways:** Only **P2** passed pre-registered floor bars with contrast still non-inverted. Approach−basal remains ≈0 on all arms. Naive contemporaneous AR residual (P1) **raises** the floor (algebraic collinearity of \(\mathrm{RR}_t\) and \(e_t\)).

```bash
PYTHONPATH=code python3 code/run_surplus_proxy_bakeoff.py
# smoke: ... --smoke
```

### D. P2 lag screen (Phase β mini — closed)

**Purpose:** Discrete lags only on the promising Phase-α arm. Floor vs contrast trade-off. **Not** lag optimization.

| Artifact | Path |
|----------|------|
| Report | `docs/SURPLUS_P2_LAG_SCREEN.md` |
| Runner | `code/run_surplus_p2_lag_screen.py` |
| Results | `results/surplus_p2lag_*.csv` + `surplus_p2lag_summary.json` |

| ID | NSRDB mean \(S\) | frac \(S\ge0.40\) | Event basal | \(d\) (ev vs NSRDB) |
|----|------------------|-------------------|-------------|---------------------|
| P2-5 | 0.270 | 0.007 | 0.296 | 0.64 |
| **P2-10** | **0.269** | **0.005** | 0.289 | **0.78** (best \(d\)) |
| P2-20 | 0.268 | 0.004 | 0.289 | 0.70 |

**Takeaway:** All three cut the R0 floor ~28%; differences among lags are small. **P2-10** remains the balanced default. Approach−basal still ≈0. No detector promotion.

```bash
PYTHONPATH=code python3 code/run_surplus_p2_lag_screen.py
# smoke: ... --smoke
```

### E. Phase β — Surplus-over-surrogate on P2-10 (**closed**)

**Purpose:** Record-specific null via phase-randomization of channel 2 only:

\[
S^{\mathrm{ex}}_t = S_t - \mathrm{median}_k(S^{\mathrm{surr}}_{k,t}),\quad
S^{z}_t = (S_t - \mu_k)/\sigma_k
\]

Diagnostic only — **not** a detector, **not** I0.

| Artifact | Path |
|----------|------|
| Report | `docs/SURPLUS_P2_SURROGATE_PHASEB.md` |
| Helpers | `code/surplus_surrogates.py` |
| Runner | `code/run_surplus_p2_surrogate_phaseb.py` |
| Tests | `tests/test_surplus_surrogates.py` (22 pass w/ proxy tests) |
| Results | `results/surplus_p2_surrogate_*` (`"smoke": false`, full 18+32) |

**Methods frozen:** proxy **P2-10**; BP \(m=3\), \(L=50\), TV; \(n_{\mathrm{surr}}=19\) phase-randomize ch2 (IAAFT optional `--method iaaft`).

**Full cohort (means of per-record stats):**

| Family | NSRDB mean | Event basal | \(\Delta\) app−bas | \(d\) |
|--------|------------|-------------|---------------------|-------|
| raw \(S\) | 0.269 | 0.289 | −0.009 | **0.78** |
| \(S^{\mathrm{ex}}\) | **−0.002** | 0.013 | −0.009 | **0.48** |
| \(S^{z}\) | −0.088 | 0.229 | −0.216 | 0.54 |

| Screening | Result |
|-----------|--------|
| Floor drop \(S^{\mathrm{ex}}\) vs \(S\) | **1.008** (≥0.20) **pass** — NSRDB mean \(S^{\mathrm{ex}}\approx 0\) |
| \(d\ge 0.6\) on \(S^{\mathrm{ex}}\) | **fail** (0.48) |
| Approach rise on \(S^{\mathrm{ex}}\) | **fail** (negative) |
| **`interesting_screening`** | **False** |

**Takeaway:** Median \(S_{\mathrm{surr}}\approx S\) on NSRDB → residual P2-10 floor is largely **null-explainable** (finite-sample independence bias). Subtracting it zeros the healthy mean but **weakens** event contrast and does **not** create pre-VF rise. **Stop** surplus-over-surrogate / \(S^{\mathrm{ex}}\) promotion on this design. No detector; ranking unchanged.

```bash
PYTHONPATH=code python3 code/run_surplus_p2_surrogate_phaseb.py          # full ~10–20 min
PYTHONPATH=code python3 code/run_surplus_p2_surrogate_phaseb.py --smoke
PYTHONPATH=code python3 -m pytest tests/test_surplus_surrogates.py -q
```

---

## Verify after clone (shipped detector suite)

```bash
cd /path/to/cctp-sddb-systemic-tau
PYTHONPATH=code python3 -m pytest \
  tests/test_opc_refinements.py \
  tests/test_opsp_holter_eval_runner.py \
  tests/test_i0_surplus_param_grid.py \
  tests/test_i0_structural_arms.py \
  tests/test_ordinal_detectors.py -q
```

Local proxy / Phase-β unit tests (optional):

```bash
PYTHONPATH=code python3 -m pytest tests/test_proxy_variants.py tests/test_surplus_surrogates.py -q
```

Diagnostic runners remain exploratory (no promotion criteria).

---

## Next tracks

| Track | Action | Priority |
|-------|--------|----------|
| **B** | Institutional Tier A / Phase-2 device-matched non-event controls | **Main product** |
| **D** | Commit/ship local surplus research as **separate** archive — not mixed into detector ranking | If archiving |
| **C** | Dual-mode Holter bake-off | Low value after S1 stop |

**Closed research (do not reopen without new design / new data):** entire Surplus/P2/Surrogate line (see `docs/SURPLUS_LINE_CLOSURE.md`); I0 θ-grids 0/16; structural R0–R5 0/6; P1/P3/P4; continuous lag grids; \(S^{\mathrm{ex}}\) / P2 as detectors; mix local surplus into locked ranking.

**Secondary research (only with explicit new design):** more orthogonal bivariate channels (RR vs local RMSSD / entropy / ectopy count); RECD Nivel-3 excess3 / three-channel synergy — not more null micro-variants on P2-10+TV.

---

## Source-of-truth documents

| Artifact | Role |
|----------|------|
| `manuscript/CCTP_SDBB_manuscript.md` (+ PDF) | Final report body (v1.4.0) |
| `docs/FINAL_DETECTOR_RANKING.md` | Ranking one-pager |
| `docs/I0_STRUCTURAL_ARMS.md` | R0–R5 + stop |
| `docs/I0_SURPLUS_PARAM_GRID.md` | 0/16 θ-grid |
| `docs/OPSP_INTEGRATED_HOLTER_EVAL.md` | I0 Holter eval |
| `docs/ORDINAL_SENSITIVITY_SPECIFICITY_TRADEOFF.md` | OPC/SDD/abs-z anchors |
| **`docs/SURPLUS_LINE_CLOSURE.md`** | **Authoritative close** — Surplus/P2/Surrogate line |
| `docs/SURPLUS_NSRDB_DIAGNOSTIC.md` | Surplus noise diagnosis |
| `docs/SURPLUS_PROXY_REDESIGN.md` | Proxy proposals + Phase plan |
| `docs/SURPLUS_PROXY_BAKEOFF.md` | Phase-α bake-off |
| `docs/SURPLUS_P2_LAG_SCREEN.md` | P2 lag ℓ∈{5,10,20} |
| `docs/SURPLUS_P2_SURROGATE_PHASEB.md` | Phase β \(S^{\mathrm{ex}}\) (not interesting) |
| `HANDOFF.md` (this file) | Session restart |

---

## Leave-off (for next agent)

**Shipped:** `origin/main` @ `6f19efb` = **v1.4.0** (abs-z primary; I0/surplus-primary stopped; GitHub + Zenodo).  
**Local archive:** `563f0cf` on `main` (ahead of origin; not pushed).

**Clinical claim:** NONE.

### Research line CLOSED — Surplus / P2 / Surrogate (2026-07-14)

Authoritative summary: **`docs/SURPLUS_LINE_CLOSURE.md`**.

| Finding | Detail |
|---------|--------|
| Floor (healthy) | Largely geometric (RR–\|ΔRR\| coupling); **P2-10** cut NSRDB mean \(S\) ~28%; Phase β zeros NSRDB \(S^{\mathrm{ex}}\) |
| Pre-VF approach | **No** systematic rise under this bivariate BP+TV design on public data |
| Phase β screening | **Not interesting** (\(d\) 0.78→0.48; approach flat/neg) |
| Promotion | **None** — not detector, not I0 |

| Exercise | Doc |
|----------|-----|
| Line close | `docs/SURPLUS_LINE_CLOSURE.md` |
| NSRDB diagnostic | `docs/SURPLUS_NSRDB_DIAGNOSTIC.md` |
| Redesign | `docs/SURPLUS_PROXY_REDESIGN.md` |
| Phase α | `docs/SURPLUS_PROXY_BAKEOFF.md` |
| Lag screen | `docs/SURPLUS_P2_LAG_SCREEN.md` |
| Phase β | `docs/SURPLUS_P2_SURROGATE_PHASEB.md` |

**Track D (archive):** local commit of surplus diagnostic stack on `main` (separate from detector ranking; **no** promotion). Not part of v1.4.0 release tag. Push/tag only if user asks.

**Locked:** detector ranking (abs-z primary, OPC companion). Do **not** reopen I0 on R0. Do **not** mutate `build_bivariate_proxy`. Do **not** promote \(S^{\mathrm{ex}}\) / P2 as detectors. Do **not** re-run lag/null micro-grids on this pipeline.

**Preferred next (product):** Track **B** institutional Tier A / Phase-2 device-matched non-event controls.  
**Optional:** push Track D archive to origin / tag if desired.  
**Only with new design:** orthogonal channels or RECD Nivel-3 — not more P2/null variants.

**Restart phrase:** read this file + `docs/SURPLUS_LINE_CLOSURE.md`. Surplus/P2/Surrogate line is **closed**. Prefer Track B.
