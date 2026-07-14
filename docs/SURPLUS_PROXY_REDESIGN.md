# Proxy redesign for synergistic surplus \(S_t\)

**Date:** 2026-07-14  
**Status:** Design / analysis only (no detector, no θ-grid, no implementation required)  
**Parent exercise:** `docs/SURPLUS_NSRDB_DIAGNOSTIC.md`  
**Scope:** Alternatives to the bivariate proxy that feeds ordinal surplus — **not** retuning I0/OPS gates.

**Clinical claim:** NONE.

---

## 0. Why this document exists

The NSRDB diagnostic established:

| Fact | Implication |
|------|-------------|
| Mean healthy \(S \approx 0.37\), p90 \(\approx 0.44\) | High floor, not sparse events |
| \(r(S,\mathrm{MI}) \approx 0.96\) | \(S_t\) ≈ ordinal dependence twin |
| Event basals **higher** than NSRDB; approach−basal \(\approx 0\) | No free pre-VF surplus rise |
| I0 θ-grid 0/16; structural arms 0/6 | Gates cannot invent contrast |

**Root cause (one line):**  
Under \(X = [z(\mathrm{RR}),\, z(|\Delta\mathrm{RR}|)]\), Bandt–Pompe factors are **structurally non-independent** in normal HRV.  
\(S_t = \mathrm{TV}(P_{\mathrm{joint}},\, P_1\otimes P_2)\) therefore measures **common ordinal coupling**, not rare pre-VF synergy.

This document asks: **what proxy redesigns could push the healthy null nearer independence, while preserving (or improving) sensitivity to pathologically relevant reorganizations?**

---

## 1. Anatomy of the current proxy

### 1.1 Construction (frozen code)

```text
rr_t
drr_t = |RR_t − RR_{t−1}|          # prepend RR_0
X = [ z_global(rr), z_global(drr) ]
π_i = Bandt–Pompe(X[:,i], m=3, delay=1)   # factors ∈ {0..5}
S_t = TV(P_joint, P1⊗P2) on last L=50
```

Source: `build_bivariate_proxy` in `code/cctp_metrics_core.py`.

### 1.2 Why healthy RR almost *must* couple these channels

Think in continuous time first, then in ordinal patterns.

**A. Algebraic / dynamical coupling**

- \(|\Delta\mathrm{RR}|\) is a **deterministic functional** of the RR path:  
  \(d_t = |rr_t - rr_{t-1}|\).  
  Given local level and successive differences, the two series share the same generating process.
- In a random-walk-like or AR(1) HRV regime, large \(|\Delta\mathrm{RR}|\) clusters co-occur with excursions of RR itself (slow wander + bursty increments).  
  Ordinal patterns on both channels therefore **co-align** more often than chance.

**B. Shared timescale of Bandt–Pompe \(m=3\), delay=1**

- Patterns are formed on **adjacent beats**.  
  Channel 1 patterns encode “RR went up/down/flat over 3 beats.”  
  Channel 2 patterns encode “|\(\Delta\mathrm{RR}\)| went up/down/flat over the same 3 beats.”  
- A sinusoid or respiratory sinus arrhythmia (RSA) produces **phase-locked** rank patterns between level and absolute increment (e.g. steepest \(|\Delta\mathrm{RR}|\) near zero-crossings of the oscillation).  
  That is **normal physiology**, not pathology.

**C. Absolute value folds signed increments**

- Using \(|\Delta\mathrm{RR}|\) destroys sign but **preserves magnitude coupling** to local RR dynamics.  
  RSA and baroreflex produce regular envelopes of \(|\Delta\mathrm{RR}|\) that are **not** independent of RR ordinal structure.

**D. Global z-score does not orthogonalize**

- \(z(\cdot)\) is an affine transform per channel; it **does not** remove cross-channel dependence.  
  Bandt–Pompe is rank-based and invariant to monotone transforms *within* a channel — so z-scoring is almost irrelevant to ordinal dependence between channels.

**E. Statistic identity**

- \(\mathrm{TV}(P_{12},\, P_1\otimes P_2)\) is a dependence measure (related to total-variation independence).  
  When healthy MI is already high, \(S\) **must** sit high. Diagnostic: \(r \approx 0.96\).

### 1.3 What we would *want* from a redesigned proxy

| Design goal | Rationale |
|-------------|-----------|
| **G1.** Healthy null: \(\mathbb{E}[S]\) nearer 0 (or a tight, stable floor ≪ 0.37) | Room for absolute / relative gates |
| **G2.** Variance of healthy \(S\) not dominated by dense short runs of normal dependence | Reduce episode rate |
| **G3.** Pre-VF / high-risk substrate still moves \(S\) (or excess-\(S\)) | Keep scientific + eventual detector value |
| **G4.** Still computable from **RR only** (or RR + cheap derived series) on public Holter | Same data constraints as pilot |
| **G5.** Compatible with discrete ordinal pipeline (\(m=3\), factors, joint \(K\)) if possible | Reuse OPS machinery; lower implementation risk |
| **G6.** Physiologically interpretable | Avoid opaque ML feature soup |

**Hard truth:** G1 and G3 can trade off. A fully independent healthy pair may also be **blind** to pathology that only reorganizes *level–increment* coupling. Redesign must aim for **orthogonalization of normal structure**, not destruction of all coupling.

---

## 2. Design space (map)

Proxy redesign can touch one or more of:

```text
  RR series
     │
     ├─► Channel selection / construction     ← main leverage (this doc)
     ├─► Pre-whitening / residualization
     ├─► Symbolization (BP vs amplitude bins vs phase)
     ├─► Delay / embedding (m, τ)
     ├─► Null model (raw S vs surplus-over-surrogate)
     └─► Statistic (TV vs residual MI vs true nested Φ₃)
```

This document prioritizes **channel construction + symbolization**, with brief notes on null/statistic (already flagged in the diagnostic).

---

## 3. Concrete proxy proposals

Each proposal: definition, healthy-null hypothesis, pre-VF hypothesis, pros/cons, feasibility, risk.

---

### P0 — Baseline (current)

\[
X = \big[z(\mathrm{RR}),\; z(|\Delta\mathrm{RR}|)\big]
\]

| | |
|--|--|
| **Healthy null** | Strong ordinal dependence (observed floor ~0.37) |
| **Pre-VF hope** | Further locking / synergy rise — **empirically weak** (flat approach−basal) |
| **Verdict** | Keep for τ_s continuity and narrative; **do not** treat raw \(S\) as FAR-primary |

---

### P1 — Orthogonal residual increment (recommended pilot #1)

**Idea:** Channel 2 should be the **innovation residual**, not raw \(|\Delta\mathrm{RR}|\).

\[
\begin{aligned}
e_t &= \mathrm{RR}_t - \widehat{\mathrm{RR}}_t \quad\text{(e.g. AR(1) or local linear predict from RR history)} \\
X &= \big[z(\mathrm{RR}),\; z(e)\big]
\quad\text{or}\quad
X &= \big[z(\mathrm{RR}),\; z(|e|)\big]
\end{aligned}
\]

**Variants:**

| ID | Residual | Notes |
|----|----------|-------|
| P1a | Global AR(1) residual of RR | Simple; one φ per record or per long basal |
| P1b | Rolling AR(1) residual (window ~50–100 beats) | Adaptive to nonstationarity |
| P1c | Residual of \(|\Delta\mathrm{RR}|\) after regressing on RR (and lag-RR) | Directly attacks structural coupling of current pair |
| P1d | Gram–Schmidt / partial out linear projection of channel-2 onto channel-1 in z-space | Continuous orthogonalization before BP |

**Healthy-null hypothesis:**  
Linear (or AR) structure of normal HRV is removed from channel 2 → ordinal patterns of residual nearer independent of RR ranks → \(\mathbb{E}[S]\) drops.

**Pre-VF hypothesis:**  
Ectopy clusters, runs of abnormal coupling, or nonlinear reorganization leave **non-AR residual structure** coordinated with RR level (or with collapse of repertoire) → residual dependence rises selectively.

| Pros | Cons |
|------|------|
| Directly targets Cause A (functional dependence) | Residualization is continuous pre-processing — mildly leaves “pure discrete” purity |
| Cheap; reuses BP + TV pipeline | Choice of AR order / window is a new hyperparameter (must freeze once) |
| Interpretable: “dependence beyond AR predictability” | Over-whitening may erase pre-VF signal if pathology is mostly AR-like rate change |
| Natural bridge to surrogates (P-null) | Global AR may fail on nonstationary Holter |

**Feasibility:** High (hours–days of code).  
**Impact potential:** High for lowering healthy floor.  
**Risk to pre-VF info:** Medium — needs empirical check on approach vs NSRDB.

---

### P2 — Delayed / multi-lag second channel (timescale decoupling)

**Idea:** Break simultaneous algebraic coupling by **time lag**.

\[
X_t = \big[z(\mathrm{RR}_t),\; z(|\Delta\mathrm{RR}|_{t-\ell})\big], \quad \ell \in \{5,10,20\}\ \text{beats}
\]

or asymmetric delays in Bandt–Pompe (channel-wise delay).

**Healthy-null hypothesis:**  
RSA and beat-to-beat dependence decorrelate over a few respiratory cycles; lag \(\ell\) near half respiratory period (~3–6 s ≈ 4–8 beats at 60–80 bpm) may **minimize** healthy MI.

**Pre-VF hypothesis:**  
Pathological coupling may be **longer-range** or **scale-free** (ectopy trains, ischemia-related HRV reorganization) and survive lag.

| Pros | Cons |
|------|------|
| Minimal code change | Lag is a free parameter; respiratory rate varies |
| Keeps both familiar physiological quantities | May only *shift* dependence, not remove it |
| Easy grid of 3–4 fixed lags for a **pilot bake-off** | Absolute value still folds sign |

**Feasibility:** Very high.  
**Impact potential:** Medium (likely partial floor reduction).  
**Best use:** Cheap **screening** arm before heavier redesigns.

---

### P3 — Level vs *signed* increment (restore direction)

\[
X = \big[z(\mathrm{RR}),\; z(\Delta\mathrm{RR})\big]
\quad\text{with}\quad
\Delta\mathrm{RR}_t = \mathrm{RR}_t - \mathrm{RR}_{t-1}
\]

**Healthy-null hypothesis:**  
Signed increments still couple to RR (negative feedback: high RR tends to be followed by negative Δ), but ordinal patterns may differ from absolute-value envelope; floor might change.

**Pre-VF hypothesis:**  
Asymmetry of accelerations/decelerations (e.g. before VF) is clinically studied; signed channel may capture **asymmetric** reorganization that \(|\cdot|\) erases.

| Pros | Cons |
|------|------|
| One-line change from P0 | Structural dependence likely **remains** (maybe even stronger via baroreflex) |
| Access to acceleration/deceleration asymmetry | May *increase* healthy MI if feedback is tight |
| Physiologically classic | Does not orthogonalize |

**Feasibility:** Trivial.  
**Impact potential:** Low–medium for floor; **medium for narrative** if asymmetry is the pre-VF story.  
**Verdict:** Include as **control arm** in any bake-off, not as sole hope for G1.

---

### P4 — Amplitude–phase (Hilbert) split

**Idea:** Separate **slow envelope** from **instantaneous phase** of RR (or of band-limited RR).

\[
\begin{aligned}
u_t &= \mathrm{RR}_t - \mathrm{trend}_t \quad\text{(band-pass or HP)} \\
a_t &= | \mathcal{H}(u)_t |,\quad
\phi_t = \arg \mathcal{H}(u)_t \\
X &= \big[z(a),\; z(\phi)\big]
\quad\text{or}\quad
X = \big[z(\mathrm{RR}),\; z(\phi)\big]
\end{aligned}
\]

For ordinal pipeline: phase can be wrapped to \([-\pi,\pi]\) then Bandt–Pompe on unwrapped local phase increments, or **phase bins** (circular quantization) instead of BP on φ.

**Healthy-null hypothesis:**  
In healthy RSA-dominated HRV, amplitude and phase of a narrow band can be **nearly independent** over short windows (or weakly dependent); phase of respiration-driven oscillation is the “clock,” amplitude the “gain.”

**Pre-VF hypothesis:**  
Loss of phase stability, phase slips, or amplitude–phase locking changes before malignant arrhythmia (autonomic decoupling / ectopy interrupting the oscillator).

| Pros | Cons |
|------|------|
| Classic signal-processing split; good physiology story | Needs band choice (e.g. 0.15–0.4 Hz HF); Holter RR is unevenly sampled → interpolate carefully |
| May achieve lower healthy dependence than level–\|Δ\| | Interpolation + Hilbert edge effects |
| Connects to known HRV phase literature | Phase of broad-band RR is noisy |

**Feasibility:** Medium.  
**Impact potential:** Medium–high if RSA is the main healthy coupling driver.  
**Risk:** Over-tuned band → dataset-specific; freeze band *a priori*.

---

### P5 — Two weakly related HRV *domains* (spectral / multiscale features as channels)

**Idea:** Replace “level vs derivative” with **two features that are not algebraically twins**.

Examples (each z-scored, short rolling estimates ~30–60 s, then BP or coarse bins):

| Channel A | Channel B | Physiological reading |
|-----------|-----------|------------------------|
| Local mean RR (or HF power) | Local RMSSD or pNN50 | Tone vs short-term variability |
| HF power | LF power | Respiratory vs vasomotor bands |
| Sample entropy (short) | Local RR | Complexity vs level |
| SD1 (Poincaré) | SD2 | Short vs long axis — **but** SD1≈RMSSD-related; still coupled |

**Preferred concrete pair for pilot:**

\[
X = \big[z(\overline{\mathrm{RR}}_{w}),\; z(\mathrm{RMSSD}_{w})\big]
\quad\text{with non-overlapping or lightly overlapping blocks}
\]

or beat-synchronous:

\[
X_t = \big[z(\mathrm{RR}_t),\; z(\mathrm{RMSSD}_{t,w})\big]
\]

with \(w\) large enough that RMSSD is not a near-deterministic function of the last 2–3 RR only.

**Healthy-null hypothesis:**  
Mean rate and short-term variability are **related but not deterministic**; ordinal independence TV may sit lower than RR–\|ΔRR\|.

**Pre-VF hypothesis:**  
Decoupling of rate and variability (or joint collapse) is a known autonomic theme pre-event.

| Pros | Cons |
|------|------|
| Channels are not strict functionals of each other at lag 0 | Rolling stats induce **strong serial dependence** and smoothing bias |
| Good clinical interpretability | Estimation noise on short windows inflates spurious joint symbols |
| Easy to explain to clinical collaborators | May blur short pre-event windows (VFDB) |

**Feasibility:** Medium.  
**Impact potential:** Medium.  
**Note:** LF/HF as pair is more independent-looking but **requires careful RR→PSD** pipeline; higher implementation cost.

---

### P6 — Ordinal pattern of RR vs **pattern of a surrogate-orthogonal complement**

**Idea:** Keep RR as channel 1; construct channel 2 as a **random projection / phase-randomized companion** constrained to match spectrum of \(|\Delta\mathrm{RR}|\) but destroy **cross**-dependence, then… wait — that would make healthy \(S\approx 0\) by construction and kill signal.

**Useful form (not a new X for raw S, but a calibration):**

Define channel pair as now (P0), but report

\[
S^{\mathrm{ex}}_t = S_t - \mathbb{E}\big[S_t^{\mathrm{surr}}\big]
\quad\text{or}\quad
\frac{S_t - \mu_{\mathrm{surr}}}{\sigma_{\mathrm{surr}}}
\]

where surrogates phase-randomize **one** channel (or IAAFT both while destroying cross-spectrum).

This is **null redesign** more than proxy redesign, but it **interacts** with proxy choice: better proxies → smaller \(\mu_{\mathrm{surr}}\) variance → cleaner excess.

| Pros | Cons |
|------|------|
| Directly answers “is this more dependent than chance?” | Costly if per-window surrogates on full Holter |
| Works with *any* proxy including P0 | Still high FAR if true healthy dependence >> surrogate (structural coupling is real, not chance) |
| Already listed in diagnostic as priority #2 | Does not fix structural non-independence of P0 — only **re-centers** it |

**Feasibility:** Medium (vectorized record-level surrogates first).  
**Impact potential:** High **if** combined with P1/P2; alone on P0 may only shift thr, not create sparse events.

**Important:** On P0, healthy dependence is **real**, not sampling noise. Surrogates that destroy that real structure will make \(S^{\mathrm{ex}}\) large **everywhere** in NSRDB — same FAR problem in new units. Surrogates help when residual dependence after a good proxy is near zero.

---

### P7 — Different symbolization (not different continuous channels)

Even with fixed \(X = [z(\mathrm{RR}), z(|\Delta\mathrm{RR}|)]\):

| ID | Change | Hope |
|----|--------|------|
| P7a | Larger embedding delay τ on channel 2 only | Decorrelate simultaneous patterns |
| P7b | \(m=2\) factors (smaller alphabet) | Coarser dependence; may lower TV floor **or** saturate |
| P7c | Amplitude-aware ordinal (Bandt–Pompe with ties / weighted BP) | Less pure rank coupling |
| P7d | Independent quantization: equal-frequency bins on each channel (not BP) | Patterns = joint histogram of levels, not ranks |
| P7e | **Permutation of increments only**: channel1 = BP(RR), channel2 = BP of *residual after removing shared permutation structure* | Advanced; research-grade |

**P7a + P2** are natural siblings.  
**P7d** changes ontology: from “ordinal pattern coupling” to “state co-occurrence” — still valid as dependence, different physiology.

| Pros | Cons |
|------|------|
| No new continuous features | May not fix Cause A if channels remain twins |
| Easy A/B in existing symbol generators | Alphabet / m changes break comparison to published OPC K=36 |

**Feasibility:** High for P7a–d.  
**Impact potential:** Low–medium alone; good **combo** with P1.

---

### P8 — Three-channel / reduced-synergy construction (true Level-3 path)

**Idea:** Move closer to nested synergy rather than 2-margin TV.

\[
Y = \big[z(\mathrm{RR}),\; z(\Delta\mathrm{RR}),\; z(\mathrm{local\ context})\big]
\]

e.g. local context = slow trend, or respiratory proxy from RR, or time-of-day surrogate.

Then use **excess3 / Φ₃** or residual mutual information:

\[
I(1;2;3) \quad\text{or}\quad
\mathrm{MI}(X_1,X_2 \mid X_3)
\]

instead of bivariate TV-to-product.

**Healthy-null hypothesis:**  
Much of RR–\|ΔRR\| dependence is **explained by a third slow factor** (respiratory drive, activity). Conditioning removes healthy floor; residual synergy may be rare.

**Pre-VF hypothesis:**  
Irreducible three-way reorganization is the original RECD Level-3 story; may finally separate narrative from pairwise MI.

| Pros | Cons |
|------|------|
| Aligns with original CCTP / RECD ontology | Higher dimension, more data per window, estimation noise |
| Directly addresses “S tracks pairwise MI” | excess3 on RR noise was historically “always on” at soft θ — needs careful design |
| Best long-term scientific fit | Not a quick pilot |

**Feasibility:** Medium–low for a clean Holter bake-off.  
**Impact potential:** High if estimation works; historically fragile on this data.  
**Verdict:** **Strategic** track, not first pilot.

---

### P9 — “Asymmetric information” channels: RR vs binary event process

**Idea:** Channel 2 is not a continuous twin but a **sparse marker**:

- Premature beat flag (if annotations available)  
- Large-deviation indicator: \(1\{|\Delta\mathrm{RR}| > \kappa\cdot\mathrm{RMSSD}\}\)  
- Run-length or local ectopy density

\[
X = \big[z(\mathrm{RR}),\; \text{binary or low-cardinality process}\big]
\]

Then joint symbols live on a **product of rich × poor alphabets**.

**Healthy-null hypothesis:**  
In clean NSRDB, binary large-deviation process is sparse → joint near independence most of the time → low mean \(S\), spikes only when large deviations cluster with specific RR patterns.

**Pre-VF hypothesis:**  
Ectopy density and coupling of ectopy to rate **do** rise before many ventricular events.

| Pros | Cons |
|------|------|
| Builds sparsity into the proxy (good for FAR) | Annotation quality varies; pure RR large-deviation threshold is a free κ |
| Clinical face validity | May collapse to “ectopy detector,” losing ordinal synergy claim |
| Orthogonal to smooth HRV dependence | Misses reorganizations without ectopy |

**Feasibility:** Medium (RR-only version high).  
**Impact potential:** Medium–high for **specificity**; different scientific claim.

---

### P10 — Poincaré-section / return-map coordinates

\[
X_t = \big[z(\mathrm{RR}_t),\; z(\mathrm{RR}_{t+1})\big]
\quad\text{or}\quad
\big[z(\mathrm{RR}_t),\; z(\mathrm{RR}_{t+1}-\mathrm{RR}_t)\big]
\]

Wait: \((\mathrm{RR}_t, \mathrm{RR}_{t+1})\) is **maximally dependent** by construction (lag-1 plot). Bad for independence null.

**Useful variant:**  
Compare **cloud shape** features as two channels only after a transform, e.g.

\[
X = \big[z(\mathrm{SD1}_{local}),\; z(\mathrm{SD2}_{local})\big]
\]

Same family as P5.  
**Avoid** raw \((\mathrm{RR}_t, \mathrm{RR}_{t+1})\) for TV-to-product — floor will be extreme.

---

## 4. Comparative matrix

Rough scores (1–5). **Not measured** — reasoned priors for prioritization.

| ID | Proposal | Feasibility | Floor↓ potential | Pre-VF info keep | Physiology clarity | Overall priority |
|----|----------|-------------|------------------|------------------|--------------------|------------------|
| **P1** | Residual / orthogonal increment | 5 | 5 | 3–4 | 4 | **A (first pilot)** |
| **P2** | Lagged \|ΔRR\| | 5 | 3 | 3 | 4 | **A (cheap screen)** |
| **P6** | Excess over surrogate (with P1) | 4 | 4* | 4 | 4 | **A (with P1)** |
| **P3** | Signed ΔRR | 5 | 2 | 3 | 5 | B (control arm) |
| **P7a** | Channel-wise BP delay | 5 | 2–3 | 3 | 3 | B (combo with P1/P2) |
| **P4** | Amplitude–phase Hilbert | 3 | 4 | 3–4 | 4 | **B (second wave)** |
| **P5** | RR mean vs RMSSD domains | 3–4 | 3 | 3 | 5 | B |
| **P9** | RR vs sparse large-dev / ectopy | 3–4 | 4 | 3 | 4 | B (different claim) |
| **P8** | 3-way / true Level-3 residual | 2 | 4 | 4 | 5 | **C (strategic)** |
| **P0** | Current | — | 1 | 2 | 4 | Baseline only |
| **P10** raw Poincaré pair | — | 0 | ? | 3 | Avoid for independence-TV |

\*P6 alone on P0 is weaker (see §3 P6).

---

## 5. Recommended program (without implementing now)

### Phase α — Proxy bake-off diagnostics (mirror NSRDB diagnostic)

**Goal:** Measure healthy floor and MI coupling under alternative \(X\), **not** alarm rates.

For each candidate proxy in a **small frozen set**:

1. NSRDB (18 × 12 h): mean/median/p90 \(S\), frac \(S\ge\) absolute cuts, \(\mathrm{corr}(S,\mathrm{MI})\), episode rate at a fixed absolute cut  
2. Event basals + approach: mean \(S\), \(\Delta S\) approach−basal  
3. Optional: Cohen’s \(d\) NSRDB vs event basal  

**Suggested first set (max 6 arms — avoid combinatorial explosion):**

| Arm | Definition |
|-----|------------|
| **R0** | P0 current \([z(\mathrm{RR}), z(|\Delta\mathrm{RR}|)]\) |
| **R1** | P3 signed \([z(\mathrm{RR}), z(\Delta\mathrm{RR})]\) |
| **R2** | P2 lag \(\ell=5\): \([z(\mathrm{RR}_t), z(|\Delta\mathrm{RR}|_{t-5})]\) |
| **R3** | P2 lag \(\ell=10\) |
| **R4** | P1b rolling AR(1) residual: \([z(\mathrm{RR}), z(e^{\mathrm{AR1}})]\) |
| **R5** | P1c residualize \(|\Delta\mathrm{RR}|\) on RR: \([z(\mathrm{RR}), z(r)]\) |

**Success criteria for “better proxy” (diagnostic, pre-registered):**

| Criterion | Pass if |
|-----------|---------|
| Floor reduction | NSRDB mean \(S\) ≤ **0.20** (or ≥ 40% relative drop vs R0) |
| Dependence deflation | mean \(r(S,\mathrm{MI})\) still high is OK; absolute MI should drop |
| Not flat pathology | Event approach or basal shows **some** separation vs NSRDB (any of: higher p90, higher mean, or approach rise in ≥60% records) |
| No quality artifact | corr with interp_frac still ~0 |

**Stop / pivot rules:**

- If **no arm** reduces NSRDB mean \(S\) by ≥25% → residualization alone insufficient; escalate to P4 or P6+P1.  
- If floor drops but **all** event contrast dies → proxy over-whitened; try milder residual (P1c with partial regression) or P9.  
- Do **not** turn a winning proxy into an I0 θ-grid until diagnostic criteria pass.

### Phase β — Null calibration on the winning proxy

Only after Phase α picks a winner:

- Record-level phase-randomization of channel 2  
- Report \(S^{\mathrm{ex}} = S - \mathrm{median}(S_{\mathrm{surr}})\)  
- Re-run NSRDB floor + event approach contrast  

### Phase γ — Optional detector (only if β looks sparse)

- Single pre-registered surplus-persist arm on \(S^{\mathrm{ex}}\)  
- Compare FAR/sens to **abs-z primary** and **OPC** anchors  
- Success bar: same structural_win region as I0 doc (FAR ≤ 2×OPC, sens ≥ 0.55) — or abandon as primary again

### Explicitly out of scope until α/β done

- Another I0-like θ Cartesian grid on P0  
- Claiming clinical utility  
- Replacing abs-z primary in the ranking table  

---

## 6. Physiological intuition cheat-sheet

| Situation | What healthy heart does | What we want proxy to do |
|-----------|-------------------------|--------------------------|
| RSA / HF oscillation | RR and \|ΔRR\| phase-locked | **Not** flag as synergy |
| Slow rate wander | Level and increments co-move | **Not** flag |
| Isolated ectopy | Local spike in \|ΔRR\| | Optional flag (P9) or mild residual bump (P1) |
| Ectopy trains / instability | Nonlinear, non-AR coupling; repertoire may collapse | **Flag** residual dependence + optional OPC co-condition |
| Pre-VF autonomic shift | Possible rate–variability decoupling or new locking | **Flag** if channels chosen to see that axis |

The design target is **not** “zero dependence forever,” but:

> **Independence under normal autonomic stationarity; dependence under non-normal residual coupling.**

---

## 7. Implementation notes (for whoever codes Phase α)

1. **Do not modify** `build_bivariate_proxy` in place — add  
   `build_bivariate_proxy_variant(rr, kind=...)` or a small `proxy_variants.py`  
   so τ_s / published detectors remain bitwise stable.  
2. Keep Bandt–Pompe \(m=3\), \(L=50\), TV definition **fixed** across arms so differences are **proxy-only**.  
3. Reuse `ordinal_synergy_surplus` and the NSRDB diagnostic runner structure (`run_surplus_nsrdb_diagnostic.py`) with a `--proxy` flag.  
4. Freeze residual windows / lags **before** looking at event sens if a detector is ever built.  
5. Document each arm’s healthy floor in a table next to R0 — same style as `SURPLUS_NSRDB_DIAGNOSTIC.md`.

---

## 8. Bottom line

1. **Current proxy is structurally guilty:** RR and \(|\Delta\mathrm{RR}|\) are not two independent “systems”; they are two readouts of one beat process. High NSRDB \(S\) is expected.  
2. **Highest leverage next step is channel redesign**, not thresholds — especially **residual / orthogonal increment (P1)** and **lagged second channel (P2)** as a cheap screen.  
3. **Surrogate excess (P6)** is powerful **after** the healthy continuous coupling is reduced; alone on P0 it mostly renames the same dense field.  
4. **Amplitude–phase (P4)** and **true Level-3 residual (P8)** are the best **second-wave** scientific bets if P1/P2 fail.  
5. **abs-z \(\tau_s\) + OPC** remain the production narrative; any new proxy is a **research track** until Phase α/β diagnostics clear pre-registered bars.

---

## 9. One-paragraph summary for handoff

The surplus FAR problem is largely a **proxy geometry** problem: \([z(\mathrm{RR}), z(|\Delta\mathrm{RR}|)]\) forces ordinal dependence in healthy Holter, so \(S_t\approx\mathrm{MI}\) sits on a ~0.37 floor. Promising redesigns prioritize **orthogonal residual increments (P1)**, **lagged second channels (P2)**, and **excess-over-surrogate on the winning proxy (P6)**; signed increments (P3) and symbolization tweaks (P7) are cheap controls; Hilbert amplitude–phase (P4) and nested three-way residuals (P8) are deeper follow-ons. Next work should be a **small frozen proxy bake-off** measuring NSRDB floor and event basal/approach contrast — not another I0 gate grid.

---

*Observational research design note. No clinical claims. No device claims.*
