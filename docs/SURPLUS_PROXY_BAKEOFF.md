# Phase α — Proxy bake-off diagnostic

**Date:** 2026-07-14  
**Status:** Closed diagnostic exercise (independent of I0 / detector ranking)  
**Parent design:** `docs/SURPLUS_PROXY_REDESIGN.md`  
**Prior diagnostic:** `docs/SURPLUS_NSRDB_DIAGNOSTIC.md`  
**Runner:** `code/run_surplus_proxy_bakeoff.py`  
**Proxy module:** `code/proxy_variants.py`  
**Artifacts:** `results/surplus_bakeoff_*`

**Clinical claim:** NONE. Not a detector. No threshold optimization.

---

## 1. Objective

Compare a **frozen set of five bivariate proxies** under identical Bandt–Pompe / TV surplus settings, and measure only:

1. **Healthy surplus floor** on NSRDB (18 × 12 h).  
2. **Contrast** between NSRDB and basals/approach of SDDB+VFDB event records.

Question: *Which proxies lower the healthy \(S\) floor without destroying all event-vs-control contrast?*

This is screening, not optimization and not promotion to I0.

---

## 2. Methods (fixed across arms)

| Item | Setting |
|------|---------|
| Encoding | Bandt–Pompe \(m=3\), delay \(=1\), factors \(\{0..5\}\), joint \(K=36\) |
| Surplus | \(S_t=\mathrm{TV}(P_{\mathrm{joint}},\,P_1\otimes P_2)\) on last \(L=50\) |
| NSRDB | 18 controls, search cap 12 h; Phase-2 basal \((0.25,\min(2,0.25T))\) h |
| Events | SDDB analytic windows + VFDB thirds (`short_db_windows`); \(n=32\) usable (`vfdb/419` no_basal) |
| High-\(S\) episode | \(S\ge 0.40\), run ≥ 5 samples, merge gap ≤ 10 (search region only on NSRDB) |
| MI | Ordinal \(I(\pi_1;\pi_2)\) on same \(L\)-windows, stride 25 → \(\mathrm{corr}(S,\mathrm{MI})\) |

**Production `build_bivariate_proxy` was not modified.** Variants live in `proxy_variants.py`.

### Frozen proxy set

| ID | Definition | Role |
|----|------------|------|
| **R0** | \([z(\mathrm{RR}),\, z(\lvert\Delta\mathrm{RR}\rvert)]\) | Current baseline |
| **P1** | \([z(\mathrm{RR}),\, z(e)]\), \(e=\) global AR(1) residual of RR | Main redesign candidate (P1a) |
| **P2** | \([z(\mathrm{RR}_t),\, z(\lvert\Delta\mathrm{RR}\rvert_{t-10})]\), \(\ell=10\) frozen | Lag screen |
| **P3** | \([z(\mathrm{RR}),\, z(\Delta\mathrm{RR})]\) signed | Simple control |
| **P4** | \([z(\mathrm{RR}),\, z(r)]\), \(r=\lvert\Delta\mathrm{RR}\rvert - \mathrm{OLS}(\lvert\Delta\mathrm{RR}\rvert\sim\mathrm{RR})\) | Residualize second channel |

Reproduce:

```bash
cd /path/to/cctp-sddb-systemic-tau
PYTHONPATH=code python3 code/run_surplus_proxy_bakeoff.py
# smoke: ... --smoke
```

### Pre-registered “promising” bars (exploratory only)

A non-R0 arm is **promising** if **(A or B)** and **C**:

| | Criterion |
|--|-----------|
| **A** | NSRDB mean \(S\) drops ≥ **25%** vs R0 |
| **B** | Fraction of windows with \(S\ge 0.40\) drops ≥ **30%** vs R0 |
| **C** | Event basals still have mean \(S\) **≥** NSRDB mean \(S\) (contrast not inverted) |

These are **not** go/no-go for detectors — only which arms deserve further research.

---

## 3. Main comparative table

Values are **means of per-record statistics** (NSRDB \(n=18\); event basals \(n=32\)).

| Proxy | NSRDB mean \(S\) | median \(S\) | p90 \(S\) | frac \(S\ge0.40\) | high-\(S\) eps / 24 h | \(r(S,\mathrm{MI})\) | Event basal mean \(S\) | Approach−basal \(\Delta S\) | Cohen’s \(d\) (event basal vs NSRDB) | Promising? |
|-------|------------------|--------------|-----------|-------------------|------------------------|----------------------|------------------------|-----------------------------|-------------------------------------|------------|
| **R0** | **0.373** | 0.371 | 0.441 | **0.301** | ~859 | 0.956 | **0.448** | −0.005 | **+1.83** | (baseline) |
| **P1** | 0.524 | 0.526 | 0.603 | 0.963 | ~137 | 0.975 | 0.619 | −0.004 | +2.47 | **No** (floor ↑) |
| **P2** | **0.269** | 0.266 | 0.325 | **0.005** | **~27** | 0.939 | **0.289** | −0.009 | **+0.78** | **Yes** |
| **P3** | 0.492 | 0.493 | 0.577 | 0.897 | ~333 | 0.976 | 0.563 | −0.007 | +2.64 | **No** (floor ↑) |
| **P4** | 0.386 | 0.384 | 0.459 | 0.395 | ~953 | 0.968 | 0.447 | −0.007 | +1.36 | **No** (floor ~same / slightly ↑) |

### Relative to R0 (NSRDB floor)

| Proxy | \(\Delta\) mean \(S\) | Drop fraction (mean) | Drop fraction (frac≥0.40) | Floor bar |
|-------|----------------------|----------------------|---------------------------|-----------|
| P1 | +0.152 | **−41%** (worse) | −220% (worse) | Fail |
| **P2** | **−0.104** | **+28%** | **+98%** | **Pass (A and B)** |
| P3 | +0.119 | −32% (worse) | −198% (worse) | Fail |
| P4 | +0.014 | −3.7% (worse) | −31% (worse) | Fail |

R0 numbers match the prior surplus diagnostic within rounding (mean ≈ 0.373, p90 ≈ 0.441, frac≥0.40 ≈ 30%, event basal ≈ 0.448, \(d\approx 1.83\)).

---

## 4. Narrative findings

### 4.1 Only P2 lowers the healthy floor

**P2 (lag \(\ell=10\))** is the sole arm that clearly moves the healthy null:

- Mean \(S\): 0.373 → **0.269** (−28% ≥ 25% bar).  
- Tail mass at absolute 0.40: 30% → **0.5%** (−98%).  
- Diagnostic episode rate: ~860 → **~27 / 24 h**.  
- Mean ordinal MI also drops (0.74 → 0.43 bits), consistent with weaker simultaneous coupling.

Interpretation: delaying the absolute-increment channel breaks **beat-synchronous** algebraic coupling that dominates R0. Healthy RSA-like level–\|Δ\| locking is partially destroyed by a 10-beat lag.

### 4.2 P2 preserves direction of event contrast, but weakens effect size

- Event basal mean \(S\) (0.289) remains **above** NSRDB (0.269) → criterion **C** holds.  
- Cohen’s \(d\) falls from **1.83 → 0.78** — still positive, no longer a large separation.  
- Approach − basal stays **≈ 0** (mean −0.009; only 38% of records rise) — same qualitative failure mode as R0: **no systematic pre-event surplus elevation**.

So P2 is a better *null geometry*, not a free pre-VF surplus detector.

### 4.3 P1 (global AR(1) residual) **raises** the floor — important negative result

A priori P1 was the “main” redesign. Empirically:

- NSRDB mean \(S\) **0.52**, frac \(S\ge0.40\) **96%**.  
- Event basals also high (0.62); \(d\) still large, but the null is *worse*.

**Why (structural, not a coding glitch):**  
With \(e_t = \mathrm{RR}_t - \hat c - \hat\phi\,\mathrm{RR}_{t-1}\), the pair

\[
X_t = \big[z(\mathrm{RR}_t),\; z(e_t)\big]
\]

is **contemporaneously algebraically dependent** (both channels contain \(\mathrm{RR}_t\)). Removing linear AR structure does **not** orthogonalize channel 2 from channel 1 at the same time index; Bandt–Pompe on that pair can *increase* joint–product TV. R0 match checks for shapes; residualization target was wrong for independence-TV.

Episode rate *falls* under P1 only because almost every window is above 0.40 and runs **merge** — not because the signal is sparse.

**Implication for redesign:** future residual arms should avoid co-time collinearity, e.g. \([z(\mathrm{RR}_{t-1}), z(e_t)]\), \([z(\widehat{\mathrm{RR}}_t), z(e_t)]\), rolling residual of **channel 2 only** without feeding \(\mathrm{RR}_t\) into both axes, or true Gram–Schmidt in a lagged basis — not the naive \([z(\mathrm{RR}_t), z(e_t)]\).

### 4.4 P3 (signed Δ) also worsens the floor

Signed increments tighten level–increment ordinal coupling (baroreflex / negative feedback). Mean \(S\) rises to ~0.49; frac≥0.40 ~90%. Useful as a **negative control**: “just drop the abs” is not a floor fix.

### 4.5 P4 (linear residual of \|ΔRR\| on RR) is nearly R0

OLS residualization of \(\lvert\Delta\mathrm{RR}\rvert\) on RR barely moves ordinal TV (mean 0.386 vs 0.373). Continuous linear partialling is almost invisible to Bandt–Pompe ranks. Residualizing the second channel **linearly** is insufficient; lag (P2) mattered more than linear orthogonalization (P4).

### 4.6 \(S\) remains an ordinal-dependence twin under all arms

Across all five proxies, mean \(r(S,\mathrm{MI})\) stays **0.94–0.98**. Changing the proxy rescales the dependence field; it does not turn TV-surplus into something other than an independence deficit.

---

## 5. Promising-arm evaluation (pre-registered)

| Arm | Floor A/B | Contrast C | Promising |
|-----|-----------|------------|-----------|
| R0 | — | Yes (event > NSRDB) | Baseline |
| P1 | Fail (floor ↑) | Yes | **No** |
| **P2** | **Pass** (−28% mean; −98% frac40) | **Yes** (event 0.289 ≥ NSRDB 0.269) | **Yes** |
| P3 | Fail | Yes | **No** |
| P4 | Fail | Yes | **No** |

**Stop / pivot notes from redesign doc, applied honestly:**

- Residualization *as implemented (P1/P4)* did **not** deliver a healthy floor drop → pure AR residualization of this form is **not** the next primary bet.  
- **Lag (P2)** did deliver floor drop with preserved (weaker) direction of contrast → **primary Phase-β candidate**.  
- No arm produces a reliable approach rise → any future detector still cannot assume “surplus climbs into VF.”

---

## 6. Recommendation

### Continue researching

1. **P2 (lagged second channel)** — only promising arm.  
   - Next (still diagnostic): small lag screen \(\ell\in\{5,10,20\}\) for floor vs \(d\) trade-off (frozen a priori, no event-tuned thr).  
   - Then Phase β: surplus-over-surrogate \(S^{\mathrm{ex}}\) **on P2**, not on R0.  
2. **Reformulated residual (P1′)** — only if redesigned to avoid contemporaneous \(\mathrm{RR}_t\)–\(e_t\) collinearity; do **not** re-run the current P1 definition.

### Do not pursue (for floor reduction under this statistic)

- **P3** signed alone.  
- **P4** global linear residual of \(\lvert\Delta\mathrm{RR}\rvert\) on RR.  
- **P1** as currently defined.

### Still out of scope until β clears

- I0-like θ-grids on any new proxy.  
- Replacing abs-z primary or OPC companion in the locked ranking.  
- Clinical claims.

---

## 7. Limitations

1. Single frozen lag \(\ell=10\) for P2 — not a lag optimization study.  
2. Global AR(1) / global OLS only — rolling residual variants not tested.  
3. Absolute high cut 0.40 for episodes is R0-centric; under P2 the distribution sits lower, so episode rates are not comparable as “alarm quality,” only as diagnostic density of absolute high states.  
4. Event \(n=32\) public records; VFDB short windows dominate approach definitions.  
5. \(S_t\) remains TV independence proxy, not nested \(\Phi_3\).  
6. No multi-testing correction; criteria were exploratory screens.

---

## 8. Artifacts

| File | Content |
|------|---------|
| `code/proxy_variants.py` | R0–P4 builders; R0 delegates to `build_bivariate_proxy` |
| `code/run_surplus_proxy_bakeoff.py` | Full bake-off runner |
| `results/surplus_bakeoff_nsrdb_per_record.csv` | Per control × proxy |
| `results/surplus_bakeoff_event_per_record.csv` | Per event × proxy |
| `results/surplus_bakeoff_proxy_summary.csv` | Aggregated table |
| `results/surplus_bakeoff_summary.json` | Methods + promising evaluation |

---

## 9. One-paragraph summary

Under identical Bandt–Pompe \(m=3\), \(L=50\), TV surplus, **only the lagged proxy P2** (\(\ell=10\)) reduced the NSRDB surplus floor (mean \(S\) 0.373→0.269; frac \(S\ge0.40\) 30%→0.5%) while keeping event basals slightly above controls (\(d\) 1.83→0.78). Global AR(1) residual (P1) and signed Δ (P3) **raised** the floor; linear residualization of \(\lvert\Delta\mathrm{RR}\rvert\) (P4) was near R0. No proxy produced a systematic approach rise. Phase α therefore nominates **P2** (and lag/surrogate follow-ups) for Phase β, and **rejects** naive contemporaneous AR residualization as a floor fix.

---

*Observational diagnostic. No clinical claims. No device claims.*
