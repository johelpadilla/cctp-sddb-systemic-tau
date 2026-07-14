# Diagnostic analysis: synergistic surplus \(S_t\) on NSRDB vs event basals

**Date:** 2026-07-14  
**Status:** Independent exploratory exercise (closed)  
**Not a detector. Not I0 / θ-grid / structural-arm continuation.**  
**Runner:** `code/run_surplus_nsrdb_diagnostic.py`  
**Artifacts:** `results/surplus_diag_*`

---

## 1. Objective

Understand **why synergistic surplus generates so much noise (high FAR) in healthy NSRDB controls**, and how that noise compares to **basal periods** of VFDB/SDDB records that contain ventricular events.

Questions:

1. Why does surplus rise frequently in healthy Holter?
2. What factors co-occur with high-surplus windows?
3. How does NSRDB surplus differ from basals of records with VF/VT (or SDDB events)?

This exercise is **diagnostic only**. No alarms, no promotion criteria, no parameter optimization.

---

## 2. Methods (fixed, inherited)

| Item | Setting |
|------|---------|
| Encoding | Bivariate Bandt–Pompe \(m=3\), delay \(=1\), factors in \(\{0..5\}\); joint \(K=36\) |
| Proxy channels | \([z(\mathrm{RR}),\; z(\lvert\Delta\mathrm{RR}\rvert)]\) (`build_bivariate_proxy`) |
| Surplus | \(S_t = \mathrm{TV}(P_{\mathrm{joint}},\, P_1\otimes P_2)\) on last \(L=50\) samples |
| NSRDB | 18 controls, search cap 12 h; basal \((0.25,\, \min(2,\,0.25\cdot T))\) h (Phase-2 family) |
| Event basals | SDDB analytic windows (`get_event_and_windows`); VFDB thirds / short-DB rules (`short_db_windows`) |
| High-\(S\) episode (diagnostic) | Absolute \(S_t \ge 0.40\), run ≥ 5 samples, merge gap ≤ 10 |
| Covariates | Rolling windows at stride 25: RR stats, diversity, entropy, MI, joint concentration |

**Honesty:** \(S_t\) is the same Level-3–consistent **ordinal independence-TV proxy** used in OPS/I0 — **not** continuous excess3 / nested \(\Phi_3\).

Reproduce:

```bash
PYTHONPATH=code python3 code/run_surplus_nsrdb_diagnostic.py
# smoke: ... --smoke
```

---

## 3. Main findings

### 3.1 Healthy NSRDB lives on a **high surplus floor**

Across **18/18** NSRDB records (12 h cap):

| Metric (per-record, then mean over records) | Value |
|---------------------------------------------|-------|
| Mean of mean \(S\) | **0.373** (range 0.339–0.407) |
| Mean of median \(S\) | ≈ 0.37 |
| Mean of p90 \(S\) | **0.441** |
| Mean of std \(S\) | 0.053 |
| Mean fraction \(S \ge 0.40\) | **30.1%** of windows |
| Mean fraction \(S \ge 0.50\) | **1.8%** of windows |
| Basal mean \(S\) vs search mean \(S\) | 0.373 vs 0.372 (**no drift**) |

**Interpretation:** In healthy Holter, surplus is **not a rare spike**. It fluctuates around a floor near **0.35–0.40**, with ~30% of windows already above the absolute diagnostic cut 0.40. That alone explains why absolute or weakly basal-relative surplus gates fire almost continuously.

Diagnostic episode counts at \(S\ge 0.40\) yield ~**860 episodes/24h** (median ~866) — this is **not** an alarm design; it is a stress test showing that “high surplus” is the **typical** state under this definition, not an anomaly.

### 3.2 High-\(S\) windows co-occur with **ordinal dependence / mild locking**, not with “noisy RR” per se

Pooled high-vs-low contrast (top 10% vs bottom 10% of \(S\) within each NSRDB record; mean Cohen’s \(d\) across 18 records):

| Feature | Mean Cohen’s \(d\) (high−low \(S\)) | Mean Pearson \(r\) with \(S\) | Direction |
|---------|-------------------------------------|-------------------------------|-----------|
| **MI_bits** (ordinal mutual information) | **+6.2** | **+0.96** | Strongly co-elevated |
| **NMI** | **+5.2** | **+0.91** | Strongly co-elevated |
| Diversity / joint support | **−1.6** | −0.43 | High \(S\) → **slightly fewer** joint symbols |
| \(H_{\mathrm{joint}}\) | −1.6 | −0.43 | Mild concentration |
| RR lag-1 autocorrelation | −1.3 | −0.36 | High \(S\) → **less smooth** RR |
| top3 joint mass | +1.0 | +0.31 | Mild pattern concentration |
| RR CV / std / range | −0.5 to −0.7 | weak | High \(S\) is **not** “high HRV chaos” |
| RMSSD | ~+0.3 | weak | Inconsistent / small |
| `interp_frac` (record-level) | — | ~0 with mean \(S\) | **Not** quality artifact |

Record-level: \(\mathrm{corr}(\overline{S},\, \overline{\mathrm{MI}})\approx 0.998\).  
So on this encoding, **\(S_t\) is essentially an ordinal dependence / MI twin**. Peaks are windows where joint Bandt–Pompe symbols **deviate from independence**, often with **mild support concentration** (partial locking), not pure high-entropy noise.

**What does *not* explain high surplus in NSRDB:**

- Interpolation / invalid-beat fraction  
- A large basal→search rise (means are flat)  
- Extreme RR range alone  

**What does:**

- Structural ordinal coupling between the two proxy channels  
- Frequent short runs of dependent symbol patterns in normal HRV  
- Mild repertoire concentration when coupling is strong  

### 3.3 Event basals have **higher** surplus than healthy controls — approach does **not** systematically raise it further

| Group | \(n\) records | Mean of mean \(S\) | Mean of p90 \(S\) |
|-------|---------------|--------------------|-------------------|
| NSRDB (all / basal / search) | 18 | **0.373** | 0.441 |
| Event basals (SDDB+VFDB) | 32\* | **0.448** | 0.505 |
| SDDB basals | 11 | 0.423 | — |
| VFDB basals | 21 | 0.461 | — |
| Event approach windows | 32 | 0.443 | — |

\*`vfdb/419` skipped (`no_basal`).

**Contrast NSRDB mean \(S\) vs event basal mean \(S\):**  
difference ≈ **−0.076**, Cohen’s \(d \approx -1.83\) (event basals **higher**).

**Approach − basal \(\Delta S\):** mean **−0.005**, median −0.009; only **41%** of event records have approach mean \(>\) basal mean (SDDB 55%, VFDB 33%).

Substrate notes:

| Substrate (record means) | NSRDB | Event basal |
|--------------------------|-------|-------------|
| Mean diversity | 0.568 | 0.542 |
| Mean MI (bits) | 0.736 | **0.939** |
| Mean RR CV | 0.061 | **0.132** |

Event basals sit on a **more dependent, higher-MI, higher RR-variability** substrate (ectopy / pathology / different short-DB sampling), **not** on a quiet low-surplus baseline that later “lights up” before VF.

---

## 4. Why surplus is so noisy on NSRDB (root-cause synthesis)

### Cause A — Structural non-independence of the bivariate proxy (primary)

The channels are \(z(\mathrm{RR})\) and \(z(\lvert\Delta\mathrm{RR}\rvert)\). These are **functionally related** on almost any RR series: large successive RR changes co-occur with local RR structure. Bandt–Pompe factors of the two channels therefore **cannot** be expected to be independent under a healthy null.  
**\(S_t \approx 0.35\) is closer to a baseline coupling bias than to “pathological synergy.”**

### Cause B — High floor + modest variance ⇒ continuous threshold crossings

With mean ≈ 0.37, std ≈ 0.05, p90 ≈ 0.44:

- Absolute gates near 0.15–0.40 are almost always on.  
- Basal-relative gates of the form \(\mathrm{mean}(S_{\mathrm{basal}})+\Delta S\) with \(\Delta S\sim 0.08\) set thr ≈ 0.45, still crossed in a large tail of **healthy** search (p90 sits near thr).  
- Persistence for a few samples does not rescue specificity when high states are common and re-enter often.

### Cause C — Basal control cannot invent contrast that is not in the series

NSRDB basal and search means are **equal**. There is no healthy “quiet basal vs active search” story for this \(S_t\).  
On event records, basals are **higher still**, so basal-relative thresholds become **stricter** on the sick cohort while healthy controls keep a lower thr — a recipe for **high FAR without selective pre-event rise**.

### Cause D — Pre-event approach is not a reliable surplus elevation regime

Mean approach \(S\) is statistically **flat to slightly lower** than basal on VFDB. Surplus-persist detectors that expect “synergy rises before VF” are fighting the empirical shape of this proxy on public Holter.

### Cause E — \(S_t\) tracks MI, not a rare synergistic collapse event

High-\(S\) windows co-track MI / NMI and mild joint concentration. That is **normal ordinal dependence**, not a sparse Level-3 “synergy storm.” Detector FAR is then the price of treating a dense dependence field as an alarm process.

---

## 5. Hypotheses for future structural work (do **not** implement here)

Ordered by diagnostic priority:

1. **Proxy redesign (highest leverage)**  
   Replace or orthogonalize \([\mathrm{RR},\,\lvert\Delta\mathrm{RR}\rvert]\) so the healthy null is nearer independence (e.g. RR vs a delayed/orthogonal feature, amplitude-free phase features, or surrogate-calibrated residual surplus).

2. **Null calibration via surrogates**  
   For each window or record, compare \(S_t\) to a distribution under phase-randomized / IAAFT RR surrogates; alarm on **excess over surrogate null**, not raw TV.

3. **True nested Level-3 statistic**  
   Continuous excess3 / \(\Phi_3\) (or a closer discrete nested mutual-info residual) may separate “synergy beyond pairwise” from the baseline pairwise dependence that dominates TV-to-product-of-margins.

4. **Multi-scale / slower aggregation**  
   \(L=50\) beat windows resolve fast HRV dependence; longer \(L\) or hop-aggregated credits may reduce micro-burst noise (structural arms already hinted FAR–sens trade-offs; not reopened here).

5. **Do not use raw \(S_t\) as FAR-primary**  
   Given floor ~0.37 and flat approach−basal, surplus alone is a poor specificity arm. Keep OPC (collapse) and abs-z \(\tau_s\) as the primary table arms; treat surplus as **substrate / narrative**, not production gate.

6. **Conditioning**  
   If surplus is revisited, require co-occurrence with **independent** structure (e.g. true collapse of repertoire *and* excess-over-null \(S\)), not \(\Delta S\) alone — knowing I-confirm already showed extreme sens cost on the old gates.

7. **Cohort caveat**  
   Event basals (esp. short VFDB) are not “healthy-like”; higher basal \(S\) and RR CV may reflect ectopy load and DB length. Any future surplus claim needs **matched control segments** of similar duration and ectopy density.

---

## 6. Implications for the closed I0 / surplus-primary line

These diagnostics **retrospectively explain** the closed surplus-primary results without reopening them:

| Prior observation | Diagnostic explanation |
|-------------------|------------------------|
| I0 / OPS FAR ≈ 40/24h on NSRDB | Healthy \(S\) floor ~0.37 + thr ~ mean+0.08 → frequent persist runs |
| No θ-grid clear advance | Parameterization cannot invent contrast; floor and variance dominate |
| Structural arms FAR-sensitive | MAD/percentile/hop reshape tails of a dense field, not a sparse event process |
| Approach not a clean surplus rise | Event approach \(\Delta S\) vs basal ≈ 0 on average |
| Event basals “noisy” too | Event basal \(S\) **higher** than NSRDB |

**Strategic takeaway:** The FAR problem is **mostly measurement / encoding / null-model**, not “need a slightly better θ.” Further micro-grids of I0-like detectors are low EV without a new statistic or proxy.

---

## 7. Deliverables checklist

| Deliverable | Location |
|-------------|----------|
| Per-record NSRDB \(S\) stats | `results/surplus_diag_nsrdb_per_record.csv` |
| High-\(S\) episodes + peak covariates | `results/surplus_diag_nsrdb_episodes.csv` |
| High vs low feature contrasts | `results/surplus_diag_nsrdb_high_vs_low.csv` (+ `_pooled.csv`) |
| Event basal / approach stats | `results/surplus_diag_event_basal_per_record.csv` |
| Cohort comparison table | `results/surplus_diag_cohort_comparison.csv` |
| Machine-readable summary | `results/surplus_diag_summary.json` |
| This report | `docs/SURPLUS_NSRDB_DIAGNOSTIC.md` |

---

## 8. Bottom line (one paragraph)

Under the project’s fixed bivariate Bandt–Pompe encoding, synergistic surplus \(S_t\) in healthy NSRDB sits on a **high floor (~0.37)** with frequent excursions above 0.40, driven almost 1:1 by **ordinal mutual dependence** (MI) and mild joint concentration — not by sparse pathology or signal-quality artifacts. Basals of public event records show **even higher** mean \(S\), and pre-event approach windows do **not** systematically elevate surplus above basal. Therefore the high FAR of surplus-primary detectors is expected: the statistic measures **common ordinal coupling** in the RR–|ΔRR| proxy, not a rare pre-VF synergistic state. Future progress requires **structural** changes (proxy, null, or true Level-3 residual), not further I0-style gate tuning.

---

*Observational research software. No clinical claims. No device / S5 claims.*
