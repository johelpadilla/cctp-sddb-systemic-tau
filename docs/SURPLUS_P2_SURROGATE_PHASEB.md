# Phase β — Surplus-over-surrogate on P2-10 (diagnostic)

**Date:** 2026-07-14  
**Status:** **Closed** full-cohort diagnostic (not a detector; not I0)  
**Parent:** `docs/SURPLUS_P2_LAG_SCREEN.md` (P2-10 default) · `docs/SURPLUS_PROXY_BAKEOFF.md` (Phase α)  
**Helpers:** `code/surplus_surrogates.py`  
**Runner:** `code/run_surplus_p2_surrogate_phaseb.py`  
**Artifacts:** `results/surplus_p2_surrogate_*` (`"smoke": false`)

**Clinical claim:** NONE.

---

## 1. Objective

On the only Phase-α promising proxy (**P2-10**), ask whether a **record-specific independence null** removes residual healthy floor while preserving (or improving) event contrast.

Define pointwise surplus-over-surrogate:

\[
S^{\mathrm{ex}}_t = S_t - \mathrm{median}_{k=1\ldots n_{\mathrm{surr}}}\bigl(S^{\mathrm{surr}}_{k,t}\bigr),
\qquad
S^{z}_t = \frac{S_t - \mu_k}{\sigma_k}
\]

where each surrogate keeps channel 1 fixed and **phase-randomizes channel 2 only** (lagged \(|\Delta\mathrm{RR}|\) after proxy z-scoring).

**Screening bars only** (pre-registered in runner; not validation):

| Criterion | Bar |
|-----------|-----|
| Floor drop | mean \(S^{\mathrm{ex}}\) NSRDB ≥ **20%** lower than mean raw \(S\) (as fraction of \|mean \(S\)\|; code: `floor_drop_frac_Sex_vs_S ≥ 0.2`) |
| Event contrast | Cohen’s \(d\) (event basal vs NSRDB) on \(S^{\mathrm{ex}}\) ≥ **0.6**, **or** mean(approach − basal) on \(S^{\mathrm{ex}}\) **> 0** |
| Logic | floor_ok **AND** (d_ok **OR** positive approach rise) → `interesting_screening` |

Not a detector; not I0; production `build_bivariate_proxy` untouched.

---

## 2. Methods (frozen)

| Item | Setting |
|------|---------|
| Proxy | **P2-10**: \([z(\mathrm{RR}_t),\, z(|\Delta\mathrm{RR}|_{t-10})]\) |
| Encoding | Bandt–Pompe \(m=3\), delay \(=1\), joint \(K=36\) |
| Surplus | \(S_t=\mathrm{TV}(P_{\mathrm{joint}},\,P_1\otimes P_2)\), \(L=50\) |
| Null | Phase-randomize **channel 2 only**; \(n_{\mathrm{surr}}=19\) (odd) |
| Optional | `--method iaaft` (not used for primary numbers) |
| NSRDB | 18 controls, search cap 12 h |
| Events | SDDB + VFDB; \(n=32\) usable basals (`vfdb/419` no_basal) |
| High windows (diagnostic only) | raw \(S\ge 0.40\); \(S^{\mathrm{ex}}\ge 0.05\); \(S^{z}\ge 2\); run ≥5, merge ≤10 |
| Seed | 20260714 (stable hash offset for surrogates) |

Reproduce:

```bash
cd /path/to/cctp-sddb-systemic-tau
PYTHONPATH=code python3 code/run_surplus_p2_surrogate_phaseb.py          # full ~10–20 min
PYTHONPATH=code python3 code/run_surplus_p2_surrogate_phaseb.py --smoke
PYTHONPATH=code python3 -m pytest tests/test_surplus_surrogates.py -q
```

Confirm full cohort: `results/surplus_p2_surrogate_summary.json` must have `"smoke": false`.

---

## 3. Full-cohort results

Values are **means of per-record statistics** (NSRDB \(n=18\); event basals \(n=32\)).  
Source: `results/surplus_p2_surrogate_summary.json` (timestamp 2026-07-14T17:09:43Z).

### 3.1 Comparison table

| Family | NSRDB mean | NSRDB p90 | frac “high”\* | high eps / 24 h | Event basal mean | Approach−basal \(\Delta\) | Cohen’s \(d\) (ev basal vs NSRDB) | frac app > bas |
|--------|------------|-----------|---------------|-----------------|------------------|---------------------------|-----------------------------------|----------------|
| **raw \(S\)** | **0.269** | 0.325 | 0.0047 (\(S\ge0.40\)) | ~27 | **0.289** | **−0.009** | **0.78** | 0.375 |
| **\(S^{\mathrm{ex}}\)** | **−0.002** | 0.054 | 0.117 (\(S^{\mathrm{ex}}\ge0.05\)) | ~599 | **0.013** | **−0.009** | **0.48** | 0.344 |
| **\(S^{z}\)** | −0.088 | — | — | ~252 | 0.229 | −0.216 | 0.54 | — |

\*High-fraction and episode rates use **different thresholds** per family; they are **not** comparable as a floor metric. Floor screening uses **mean** levels.

**Surrogate floor check (NSRDB):** mean of mean \(S_{\mathrm{surr}}^{\mathrm{med}} \approx 0.271\) vs mean of mean \(S \approx 0.269\) → on healthy records, the phase-null median nearly **matches** raw \(S\). Residual P2-10 floor is largely **finite-sample / null independence bias**, not excess coupling beyond the ch2-phase null.

### 3.2 Screening decision

| Check | Value | Pass? |
|-------|-------|-------|
| Floor drop \(S^{\mathrm{ex}}\) vs raw \(S\) | **1.008** (≥ 0.20) | **Yes** — NSRDB mean \(S^{\mathrm{ex}}\) ≈ 0 |
| Cohen’s \(d\) on \(S^{\mathrm{ex}}\) | **0.48** (≥ 0.60) | **No** |
| Approach − basal on \(S^{\mathrm{ex}}\) | **−0.009** (> 0) | **No** |
| **`interesting_screening`** | | **False** |

Floor criterion alone is not enough. Excess surplus removes the healthy mean floor but **weakens** event–NSRDB separation (\(d\): 0.78 → 0.48) and does **not** create a pre-event rise.

### 3.3 Sanity vs prior P2-10

Raw-\(S\) row matches the lag-screen / bake-off P2-10 path within rounding (mean \(S\) 0.269, frac≥0.40 ~0.005, event basal 0.289, \(d\) 0.78). Unit test `test_p2_10_raw_matches_bakeoff_path` enforces path parity.

---

## 4. Narrative answers (honest)

### Does \(S^{\mathrm{ex}}\) lower the healthy floor beyond raw P2-10?

**Yes, on the mean.** NSRDB mean \(S^{\mathrm{ex}} \approx 0\) while mean raw \(S \approx 0.27\). Phase-randomizing channel 2 alone accounts for essentially the entire residual P2-10 floor. That is scientifically useful: it identifies the residual as **null-like under ch2 phase scramble**, not as a rare healthy synergy state.

**Caveat on “high” rates:** At the diagnostic threshold \(S^{\mathrm{ex}}\ge 0.05\), a large fraction of NSRDB windows still exceed the bar (mean frac ~0.12; ~599 episodes/24 h). That threshold was a provisional screen cut, **not** a calibrated specificity target. Do not read episode rates as a failed floor drop; read **mean recenter** as the floor result.

### Does event contrast improve?

**No.** Cohen’s \(d\) (event basal vs NSRDB) falls from **0.78** (raw) to **0.48** (\(S^{\mathrm{ex}}\)) and **0.54** (\(S^{z}\)). Event basal \(S^{\mathrm{ex}}\) is only ~0.013 above a ~0 NSRDB mean — a small absolute excess that does not clear the \(d \ge 0.6\) bar.

### Is there an approach − basal rise?

**No.** Mean \(\Delta\) on raw \(S\), \(S^{\mathrm{ex}}\), and \(S^{z}\) is **negative** (≈ −0.009 for \(S\) and \(S^{\mathrm{ex}}\); more negative for \(S^{z}\)). Fraction of records with approach > basal remains ~⅓–⅜. Same flat/negative pattern as Phase α and the lag screen.

### Screening verdict

**Not interesting** under pre-registered bars. Do **not** promote \(S^{\mathrm{ex}}\) or \(S^{z}\) to detector / I0 work on this pilot.

---

## 5. Interpretation

1. **P2-10 residual floor is mostly null-explainable.** Median surrogate \(S\) ≈ raw \(S\) on NSRDB ⇒ finite-sample independence bias (and/or marginal structure preserved under ch2 phase scramble) dominates healthy TV surplus.
2. **Subtracting that null does not create a pre-VF signal.** Approach windows are not systematically above basal on excess metrics.
3. **Weak residual event elevation** (basal \(S^{\mathrm{ex}}\) slightly above NSRDB) is consistent with slightly higher average coupling or nonstationarity on event records, not with a lead-time structure usable as a detector on this design.
4. **Consistent with structural stop (R0–R5, I0 grids):** further micro-nulls on the same BP/TV pipeline are unlikely to yield a free pre-event rise without a **new** statistic, windowing, or data design.

---

## 6. What this does *not* claim

- No clinical, FDA, or deployability claim.
- No change to locked detector ranking (abs-z primary; OPC companion).
- No I0 promotion; no threshold promotion for \(S^{\mathrm{ex}}\) / \(S^{z}\).
- No mutation of production `build_bivariate_proxy`.
- IAAFT was implemented but **not** required for the primary phase-null conclusion; optional sensitivity only.
- Episode rates at thr 0.05 / \(z=2\) are diagnostic instrumentation, not FAR claims.

---

## 7. Decision / stop rules

| Action | Decision |
|--------|----------|
| Archive Phase β results + this doc | **Yes** (local research; optional commit) |
| Promote \(S^{\mathrm{ex}}\) / \(S^{z}\) detector or I0 arm | **No** |
| Re-run I0 θ-grids or R0–R5 structural arms | **No** |
| Continuous lag grids / more P1–P4 variants | **No** without new design hypothesis |
| Institutional Tier A non-event controls | **Still** the main product path (Track B) |

**Research close for surplus-null line on P2-10:** floor problem is largely diagnosed and null-subtractable; **pre-event temporal contrast is not recovered.** Further surplus-primary work needs a new design, not another null flavor on the same proxy.

---

## 8. Artifacts

| Path | Role |
|------|------|
| `code/surplus_surrogates.py` | phase / IAAFT ch2; `surplus_from_proxy_matrix` |
| `code/run_surplus_p2_surrogate_phaseb.py` | Full Phase-β runner |
| `tests/test_surplus_surrogates.py` | Shapes + raw P2-10 parity |
| `results/surplus_p2_surrogate_summary.json` | Aggregate + screening |
| `results/surplus_p2_surrogate_summary.csv` | Compact table |
| `results/surplus_p2_surrogate_nsrdb_per_record.csv` | Per-control |
| `results/surplus_p2_surrogate_event_per_record.csv` | Per-event basal / approach |
