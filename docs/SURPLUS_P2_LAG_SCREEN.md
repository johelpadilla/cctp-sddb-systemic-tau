# Phase β mini — P2 lag screen (diagnostic)

**Date:** 2026-07-14  
**Status:** Closed mini-diagnostic (not a detector; not lag optimization)  
**Parent:** `docs/SURPLUS_PROXY_BAKEOFF.md` (Phase α)  
**Runner:** `code/run_surplus_p2_lag_screen.py`  
**Proxy module:** `code/proxy_variants.py` (`P2-5`, `P2-10`, `P2-20`)  
**Artifacts:** `results/surplus_p2lag_*`

**Clinical claim:** NONE.

---

## 1. Objective

Within the only Phase-α promising arm (**P2**), compare three **frozen** lags

\[
X_t = \bigl[z(\mathrm{RR}_t),\; z\bigl(|\Delta\mathrm{RR}|_{t-\ell}\bigr)\bigr],
\quad \ell \in \{5,\,10,\,20\}
\]

and measure the trade-off between:

1. **Healthy surplus floor** (NSRDB), and  
2. **Event contrast** (basal vs NSRDB; approach−basal).

Not a continuous grid, not I0, not detector thresholds.

---

## 2. Methods (identical to Phase α except lag)

| Item | Setting |
|------|---------|
| Encoding | Bandt–Pompe \(m=3\), delay \(=1\), joint \(K=36\) |
| Surplus | \(S_t=\mathrm{TV}(P_{\mathrm{joint}},\,P_1\otimes P_2)\), \(L=50\) |
| NSRDB | 18 controls, search cap 12 h |
| Events | SDDB + VFDB; \(n=32\) usable (`vfdb/419` no_basal) |
| High-\(S\) episode | \(S\ge 0.40\), run ≥ 5, merge ≤ 10 |
| Reference within screen | **P2-10** (Phase-α lag) |

Reproduce:

```bash
cd /path/to/cctp-sddb-systemic-tau
PYTHONPATH=code python3 code/run_surplus_p2_lag_screen.py
# smoke: ... --smoke
```

`build_bivariate_proxy` was **not** modified. Lag variants are registered as `P2-5` / `P2-10` / `P2-20`; production `P2` remains ℓ=10.

---

## 3. Comparative table

Values are **means of per-record statistics** (NSRDB \(n=18\); event basals \(n=32\)).

| ID | \(\ell\) | NSRDB mean \(S\) | median \(S\) | p90 \(S\) | frac \(S\ge0.40\) | high-\(S\) eps / 24 h | \(r(S,\mathrm{MI})\) | Event basal mean \(S\) | Approach−basal \(\Delta S\) | Cohen’s \(d\) (event basal vs NSRDB) |
|----|----------|------------------|--------------|-----------|-------------------|------------------------|----------------------|------------------------|-----------------------------|-------------------------------------|
| **P2-5** | 5 | 0.270 | 0.268 | 0.328 | 0.0071 | ~38 | 0.941 | **0.296** | −0.011 | 0.64 |
| **P2-10** | 10 | 0.269 | 0.266 | 0.325 | 0.0047 | ~27 | 0.939 | 0.289 | −0.009 | **0.78** |
| **P2-20** | 20 | **0.268** | 0.266 | **0.324** | **0.0041** | **~25** | 0.938 | 0.289 | −0.006 | 0.70 |

### Relative to Phase-α R0 (from bake-off; not re-run here)

R0: mean \(S\approx 0.373\), frac≥0.40 \(\approx 0.301\), \(d\approx 1.83\), eps/24h \(\sim 860\).

| ID | Drop mean \(S\) vs R0 | Drop frac≥0.40 vs R0 | \(d\) retained vs R0’s 1.83 |
|----|----------------------|----------------------|------------------------------|
| P2-5 | −28% | −98% | 0.64 (weaker) |
| P2-10 | −28% | −98% | 0.78 |
| P2-20 | −28% | −99% | 0.70 |

### Relative to P2-10 (within-screen reference)

| ID | \(\Delta\) mean \(S\) | \(\Delta\) frac≥0.40 | \(\Delta\,d\) |
|----|----------------------|----------------------|---------------|
| P2-5 | +0.6% (slightly higher floor) | +50% relative (0.007 vs 0.005) | −0.14 |
| P2-10 | 0 | 0 | 0 |
| P2-20 | −0.3% (slightly lower floor) | −14% relative | −0.08 |

**Sanity check:** P2-10 numbers match Phase-α arm **P2** to machine precision (same code path, same lag).

---

## 4. Narrative answers (honest)

### Is there a clearly superior lag?

**No.** Differences among \(\ell=5,10,20\) are **small** compared with the P2 vs R0 gap. All three:

- Cut NSRDB mean \(S\) by ~28% vs R0  
- Nearly eliminate windows with \(S\ge 0.40\) (~0.4–0.7% vs ~30%)  
- Keep event basal mean \(S\) ≥ NSRDB mean (contrast not inverted)  
- Show **approach−basal ≈ 0** (no free pre-VF rise)

### Does longer lag (20) keep reducing the floor or lose contrast?

- Floor: **marginal further reduction** (mean \(S\) 0.268 vs 0.269; frac≥0.40 0.004 vs 0.005; episodes slightly fewer).  
- Contrast: **does not improve** — Cohen’s \(d\) **falls** slightly (0.70 vs 0.78). Event basal level is essentially flat (0.289).  
- So ℓ=20 is not harmful, but it is **not a free lunch**: a tiny quieter floor with a slightly weaker separation statistic.

### Is short lag (5) insufficient to lower the floor?

**No.** ℓ=5 already delivers nearly the full P2 floor reduction vs R0. It is **not** “insufficient.” Relative to ℓ=10 it is only a bit noisier (higher frac≥0.40 and episode rate, slightly higher mean \(S\)) and has the **weakest** \(d\) in this set (0.64).

### Recommendation (balance)

| Choice | Role |
|--------|------|
| **P2-10 (ℓ=10)** | **Default reference** — best \(d\) in the discrete set, mid floor, matches Phase α |
| P2-20 | Acceptable if one prioritizes minimal high-\(S\) rate; no clear contrast gain |
| P2-5 | Acceptable; slightly noisier floor, weaker \(d\) — no reason to prefer over 10 |

**Do not** treat this as optimized lag selection. Three points only; no continuous search; no detector thresholds; production ranking unchanged.

---

## 5. Interpretation (why lags look similar)

Lagging \(|\Delta\mathrm{RR}|\) breaks **contemporaneous** ordinal coupling that drives the healthy MI/surplus floor under R0. Once that contemporaneous joint structure is removed, further delay (5→10→20 beats) mostly reshuffles which increment window is paired with \(\mathrm{RR}_t\); Bandt–Pompe ranks over \(L=50\) absorb much of that change, so cohort-level \(S\) statistics stay close.

The residual event-vs-NSRDB gap under all P2 lags is real but **modest** (\(d\sim 0.6\text{–}0.8\) vs R0’s \(\sim 1.8\)) and is **level/contrast**, not a rising pre-event approach signal.

---

## 6. What this does *not* authorize

- Promoting P2 (any lag) into I0 / detector ranking  
- A continuous lag grid or threshold optimization  
- Claims of pre-VF “detection” from surplus under P2  
- Reopening structural R0–R5 or I0 θ-grids on these proxies without a new, pre-registered design

Optional next research (still diagnostic): surplus-over-surrogate / \(S^{\mathrm{ex}}\) on **P2-10 only**, or institutional non-event controls (Track B).
