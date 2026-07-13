# Native Ordinal Alarm Detectors for Systemic Tau / RECD

**Status:** Methodological proposal (formalization only)  
**Version:** 0.1.0-draft  
**Date:** 2026-07-12  
**Scope:** Conceptual alternatives to the frozen continuous abs-z detector  
**Not in scope:** Empirical bake-off, parameter retuning of production CCTP, fusion of the two options, clinical claims

---

## 0. Motivation and ontological constraint

### 0.1 Current detector (frozen production baseline)

The CCTP pilot currently alarms on a **continuous scalar metric** \(m(t)\) (e.g. \(\tau_s(t)\) or RECD excess) via a basal z-score:

\[
z(t) = \frac{m(t) - \mu_{\mathrm{basal}}}{\sigma_{\mathrm{basal}} + \varepsilon},
\qquad
A_{\mathrm{abs\text{-}z}}(t) = \mathbf{1}\{|z(t)| \ge 2\}
\]

with a sustained-run requirement of \(\ge 3\) consecutive alarmed samples (`detect_lead_time` in `code/cctp_metrics_core.py`).

This rule is empirically useful and is **frozen** for discovery/validation reproducibility (\(\theta_3=0.08\), high-threshold \(0.65\), \(W_{\tau}=101\), abs-z \(\ge 2 \times 3\)). It is **not** ontologically native to the ordinal construction of Systemic Tau and RECD: the decision depends on the first two moments of a continuous trajectory of an already-aggregated metric, not on the discrete symbol stream that generates ordinal structure.

### 0.2 Design principles for alternatives

Both options developed below must:

1. Operate on a **discrete ordinal symbol stream** \(\sigma_t \in \Sigma\), not on \(\mu,\sigma\) of the raw series or of \(m(t)\).
2. Respect Bandt–Pompe / RECD **topological** information (ranks, patterns, joint configurations).
3. Remain **completely independent** of each other (no fused rule in this document).
4. Stay free of clinical, FDA, or deployability claims.

### 0.3 Input: ordinal symbol stream

**Definition 0.1 (Symbol alphabet).**  
Let \(x_t \in \mathbb{R}^d\) be a (uni- or multi-variate) observation stream already mapped to ordinal symbols by a fixed embedding \((m,\tau)\):

\[
\sigma_t = \pi\bigl(x_t, x_{t+\tau}, \ldots, x_{t+(m-1)\tau}\bigr) \in \Sigma,
\qquad
|\Sigma| = K < \infty.
\]

Typical cardiac settings in this project:

| Encoding | \(K\) | Notes |
|----------|------:|-------|
| Univariate Bandt–Pompe, \(m=3\) | \(3! = 6\) | Standard BP |
| Bivariate independent pair of BP symbols | \(6^2 = 36\) | Joint code of two channels |
| RECD relation codes on pairs | \(3\) per pair (EQ/GT/LT) | Aligns with \(\Phi_2\) machinery |
| Joint \(N\)-variable BP rows | \((m!)^N\) | Full multivariate symbol |

**Assumption (fixed encoding).** The alphabet \(\Sigma\) and the map \(x \mapsto \sigma\) are **frozen** before detection. Detectors never re-estimate moments of \(x\) or of continuous \(\tau_s\); they only consume \(\{\sigma_t\}\).

**Notation.** Time indices are discrete (beat index or embedding step). Windows are left-closed, right-closed index sets of consecutive symbols. Empty-window edge cases return “no alarm”.

---

# Option 1 — Ordinal Persistence Collapse Detector

**Working name (manuscript):** *Ordinal Persistence Collapse* (OPC)  
**RECD alignment:** Collapse of ordinal support + sustained persistence \(\leftrightarrow\) \(\Phi_2\)-style locking and path toward deeper conjunction (\(\Phi_3\))

## 1.1 Variables

### 1.1.1 Observation window

Fix an observation length \(L \in \mathbb{N}\), \(L \ge 2\).

\[
W_t := \{\sigma_{t-L+1}, \sigma_{t-L+2}, \ldots, \sigma_t\}
\quad\text{for } t \ge L-1.
\]

For \(t < L-1\), all quantities below are undefined and the detector outputs no alarm.

### 1.1.2 Ordinal support and diversity

**Definition 1.1 (Support).**  
\[
\operatorname{supp}(W_t) := \{\, s \in \Sigma : \exists\, u \in W_t,\; \sigma_u = s \,\}.
\]

**Definition 1.2 (Support diversity — primary).**  
\[
D_t := \frac{\bigl|\operatorname{supp}(W_t)\bigr|}{K} \in \Bigl\{\tfrac{1}{K},\tfrac{2}{K},\ldots,1\Bigr\}
\quad\text{(or \(0\) if \(W_t\) empty)}.
\]

\(D_t\) is purely combinatorial: it counts **how many distinct ordinal states** are active, normalized by alphabet size. It does **not** use frequencies beyond presence, nor amplitude.

**Definition 1.3 (Shannon diversity — optional diagnostic).**  
Let \(p_t(s) = L^{-1}\#\{u\in W_t:\sigma_u=s\}\). Then

\[
H_t := -\sum_{s\in\Sigma} p_t(s)\,\log K\, p_t(s)
\quad\text{(base-\(K\) entropy; \(H_t\in[0,1]\))}.
\]

The **alarm rule of Option 1 uses \(D_t\) only**. \(H_t\) may be reported for diagnostics but is not required for the predicate.

### 1.1.3 Persistence (sustained low-diversity run)

**Definition 1.4 (Low-diversity indicator).**  
Given a diversity threshold \(\theta_D \in (0,1)\),

\[
\ell_t := \mathbf{1}\{D_t \le \theta_D\}.
\]

**Definition 1.5 (Persistence length — primary).**  
\[
R_t := \max\Bigl\{ r \in \mathbb{N}_0 :
  t-r+1 \ge L-1,\;
  \ell_{t-r+1}=\ell_{t-r+2}=\cdots=\ell_t=1
\Bigr\},
\]

with \(R_t=0\) if \(\ell_t=0\).

Thus \(R_t\) is the number of **consecutive observation endpoints** ending at \(t\) for which the ordinal support has collapsed below \(\theta_D\). This is a discrete run-length on the collapse indicator, not a z-score run.

### 1.1.4 Auxiliary: dominant-symbol run (Φ₂-adjacent diagnostic)

**Definition 1.6 (Dominant-symbol run).**  
Let \(s^\star_t \in \arg\max_{s\in\Sigma} \#\{u\in W_t:\sigma_u=s\}\) (ties broken by fixed total order on \(\Sigma\)). Define the raw symbol run length

\[
\rho_t := \max\{ r \ge 1 : \sigma_{t-r+1}=\cdots=\sigma_t \},
\]

or \(0\) if \(t\) invalid. \(\rho_t\) is **not** part of the primary alarm predicate; it is available for RECD interpretation (locking of a single ordinal state).

## 1.2 Alarm rule

**Definition 1.7 (OPC alarm).**  
Declare an alarm at time \(t\) if and only if

\[
\boxed{
A_{\mathrm{OPC}}(t) \;=\;
\mathbf{1}\{ D_t \le \theta_D \}
\;\wedge\;
\mathbf{1}\{ R_t \ge \theta_R \}
}
\]

Equivalently: the support diversity has been **jointly low and sustained** for at least \(\theta_R\) consecutive window endpoints.

**First-alarm time** (pre-event search, if an event time \(t_E\) and search start \(t_0\) are given):

\[
t^\star_{\mathrm{OPC}}
=
\min\bigl\{ t : t_0 < t < t_E,\; A_{\mathrm{OPC}}(t)=1 \bigr\}
\]

(lead-time \(t_E - t^\star\) if the min exists; otherwise no detection). This timing wrapper is identical in form to the production detector’s “first sustained alarm,” but the **predicate** is ordinal.

## 1.3 Parameters and suggested initial ranges

| Symbol | Role | Suggested initial range | Notes |
|--------|------|-------------------------|--------|
| \(L\) | Observation window length (symbols) | \(5,\,8,\,10\) (try \(L=8\) first) | Must be \(\ll\) basal length; large enough to estimate support |
| \(\theta_D\) | Diversity collapse threshold | \(0.25\)–\(0.45\) (try \(0.35\)) | Fraction of alphabet; for \(K=6\), \(D\le 0.33\) means \(\le 2\) symbols |
| \(\theta_R\) | Minimum persistence (consecutive low-\(D\) windows) | \(4\)–\(6\) (try \(5\)) | Analogous spirit to “sustained” abs-z runs, but on \(\ell_t\) |
| \(K\) | Alphabet size | Fixed by encoding | Not a free tuning knob once encoding is frozen |
| \((m,\tau)\) | Embedding | Project defaults (e.g. \(m=3,\tau=1\)) | Inherited from RECD/τ_s pipeline; not re-optimized here |

**Calibration note (honest).** Ranges follow the seed discussion in `Contexto.docx`. They are **starting proposals** for later empirical study, not optimized thresholds and not replacements for the frozen abs-z configuration.

### 1.3.1 Edge cases

1. \(t < L-1\): \(A_{\mathrm{OPC}}(t)=0\).
2. \(K=1\): degenerate alphabet; detector disabled.
3. All symbols identical in \(W_t\): \(D_t=1/K\), \(\ell_t=1\) if \(1/K\le\theta_D\).
4. Ties in dominant symbol (diagnostic only): fixed total order on \(\Sigma\).

## 1.4 Algorithm (pseudocode)

```
Algorithm OPC_DETECT(σ[0..T-1], L, θ_D, θ_R, K):
  Input:  symbol stream σ_t ∈ {0,...,K-1}
  Output: alarm[0..T-1] ∈ {0,1}, D[·], R[·]

  alarm ← zeros(T); D ← NaN(T); R ← zeros(T)
  run ← 0

  for t ← 0 to T-1:
    if t < L-1:
      continue
    W ← σ[t-L+1 .. t]
    n_unique ← |set(W)|
    D[t] ← n_unique / K
    if D[t] ≤ θ_D:
      run ← run + 1
    else:
      run ← 0
    R[t] ← run
    if (D[t] ≤ θ_D) and (R[t] ≥ θ_R):
      alarm[t] ← 1

  return alarm, D, R
```

**Complexity.** \(O(T L)\) naïve, or \(O(T)\) with a sliding multiset / frequency table of size \(K\).

**Reference implementation (optional, non-production):**  
`code/ordinal_detectors/opc_detector.py` — pure functions on symbol arrays; **not** wired into `detect_lead_time`.

## 1.5 Advantages vs abs-z ≥ 2 (sustained)

| Aspect | abs-z on continuous \(m(t)\) | OPC |
|--------|----------------------------|-----|
| Ontology | Moments of an aggregated real series | Discrete support of ordinal states |
| Scale / units | Sensitive to basal \(\sigma\); device gain affects \(m\) | Invariant to monotone transforms of raw \(x\) (via ordinal encoding) |
| What is “abnormal” | \(|m-\mu|/\sigma\) large | Few ordinal patterns + that state **persists** |
| RECD link | Indirect (alarm on τ_s / excess amplitude) | Direct: collapse + persistence \(\approx\) \(\Phi_2\) locking |
| Interpretability | “Metric left basal band” | “System trapped in a small ordinal repertoire for too long” |

## 1.6 Limitations and open points

1. **Threshold sensitivity.** \(\theta_D\) interacts with \(K\) and \(L\); a fixed \(\theta_D\) is not automatically portable across encodings without re-expressing as absolute support size \(n_{\mathrm{unique}}\le n^\star\).
2. **Gradual reorganization.** If the system reorganizes while keeping moderate diversity, OPC may not fire (by design: it targets **collapse**, not distribution drift).
3. **Window alignment.** Sliding every symbol step can yield highly autocorrelated \(R_t\); reporting may prefer non-overlapping evaluation or hop size \(h\ge 1\).
4. **No amplitude cue.** Rare but high-impact amplitude events that do not change ranks will not appear in \(\sigma_t\) (feature of ordinal methods, not a bug of OPC alone).
5. **Empirical performance unknown** on SDDB/VFDB/NSRDB under this goal; no superiority claim is made vs frozen abs-z.

## 1.7 Coherence with RECD / Systemic Tau

- **Ordinal paradigm:** All inputs are Bandt–Pompe / RECD symbols.
- **\(\Phi_2\) (persistent relations):** \(R_t\) is a global, windowed persistence of **collapsed repertoire**, the system-level counterpart of pairwise relation locking in `compute_phi2`.
- **Path toward \(\Phi_3\):** Prolonged low diversity concentrates mass on few joint configurations—necessary (not sufficient) for deeper conjunction structure.
- **τ_s:** OPC does **not** replace τ_s as a relational score; it proposes a **native alarm layer** on the same symbol substrate that underlies τ_s / RECD, avoiding a second-order continuous z-score on the score itself.

**Option 1 is complete as a standalone detector. No coupling to Option 2 is defined.**

---

# Option 2 — Symbolic Distribution Divergence Detector

**Working name (manuscript):** *Symbolic Distribution Divergence* (SDD)  
**RECD alignment:** Structural reorganization of the ordinal pattern law (change of empirical measure on \(\Sigma\)), independent of continuous amplitude

## 2.1 Variables and distributions

### 2.1.1 Current-window empirical law

Fix current-window length \(L_c \ge 2\).

\[
W^{\mathrm{cur}}_t := \{\sigma_{t-L_c+1},\ldots,\sigma_t\},
\qquad
P_t(s) := \frac{1}{L_c}\#\{u\in W^{\mathrm{cur}}_t:\sigma_u=s\},
\quad s\in\Sigma.
\]

\(P_t\) is a probability vector on the **finite** alphabet \(\Sigma\).

### 2.1.2 Basal empirical law

**Definition 2.1 (Basal segment).**  
Let \(B = \{t_b^{\mathrm{start}},\ldots,t_b^{\mathrm{end}}\}\) be a pre-specified basal index set (or time window mapped to indices), with length \(L_b = |B| \ge L_c\), and **no overlap** with the detection search region.

\[
P_{\mathrm{basal}}(s)
:=
\frac{1}{L_b}\#\{u\in B:\sigma_u=s\},
\quad s\in\Sigma.
\]

**Estimation policies (choose one; freeze before evaluation):**

| Policy | Definition | Use when |
|--------|------------|----------|
| **Fixed basal (recommended default)** | Single \(P_{\mathrm{basal}}\) from the first eligible basal block (e.g. first 1–2 h of Holter / early pre-event segment) | Stable controls; matches CCTP basal-window spirit without using \(\mu,\sigma\) of \(m\) |
| Sliding basal | \(P_{\mathrm{basal},t}\) from \([t-L_b-G,\,t-G]\) with gap \(G\ge 0\) | Slow nonstationarity; more parameters |
| Pooled multi-block basal | Average of \(M\) early blocks | Longer recordings with early artifacts |

**This formalization adopts Fixed basal as the primary rule.** Sliding basal is noted only as a robustness variant, not as a second detector.

### 2.1.3 Support handling

Let \(S_{\mathrm{basal}} = \{s: P_{\mathrm{basal}}(s)>0\}\), \(S_t=\{s:P_t(s)>0\}\).

- **Total variation** (recommended) needs no smoothing.
- **KL** (secondary) needs additive smoothing when supports differ (see §2.2.2).

## 2.2 Divergence measures

### 2.2.1 Recommended: Total Variation (TV)

**Definition 2.2 (Total variation distance).**  

\[
\mathrm{TV}\bigl(P_t,\,P_{\mathrm{basal}}\bigr)
:=
\frac{1}{2}\sum_{s\in\Sigma}\bigl|P_t(s)-P_{\mathrm{basal}}(s)\bigr|
\in [0,1].
\]

**Properties used here:**

- Metric on probability simplices; symmetric.
- Bounded in \([0,1]\) → interpretable threshold \(\theta_{\mathrm{TV}}\in(0,1)\).
- No logarithm; **robust** when \(P_{\mathrm{basal}}(s)=0\) but \(P_t(s)>0\) (no infinite divergence).
- Purely a function of two discrete distributions on \(\Sigma\).

### 2.2.2 Secondary (alternative only): Kullback–Leibler

**Definition 2.3 (Smoothed KL — not the primary alarm).**  
For \(\varepsilon>0\) small (e.g. \(\varepsilon=1/L_b\)),

\[
\tilde P(s)=\frac{P(s)+\varepsilon}{1+K\varepsilon},
\qquad
\mathrm{KL}\bigl(\tilde P_t \,\|\, \tilde P_{\mathrm{basal}}\bigr)
=
\sum_{s\in\Sigma}\tilde P_t(s)\,\log\frac{\tilde P_t(s)}{\tilde P_{\mathrm{basal}}(s)}.
\]

KL is **asymmetric**, unbounded, and more fragile under sparse counts. It may be reported as a diagnostic companion; **the recommended alarm uses TV only**.

## 2.3 Alarm rule

**Definition 2.4 (SDD alarm — TV).**  

\[
\boxed{
A_{\mathrm{SDD}}(t)
=
\mathbf{1}\Bigl\{
  \mathrm{TV}\bigl(P_t,\,P_{\mathrm{basal}}\bigr)
  \ge
  \theta_{\mathrm{TV}}
\Bigr\}
}
\]

Optional **sustainment** (still Option 2 only; not a fusion with OPC): require the TV predicate for \(\theta_S\) consecutive endpoints,

\[
A_{\mathrm{SDD}}^{\mathrm{sust}}(t)
=
\mathbf{1}\Bigl\{
  \min_{0\le j<\theta_S}
  \mathrm{TV}\bigl(P_{t-j},\,P_{\mathrm{basal}}\bigr)
  \ge \theta_{\mathrm{TV}}
\Bigr\},
\]

with default \(\theta_S=1\) (instantaneous) or \(\theta_S\in\{2,3\}\) if short-window noise is a concern. Sustainment here still uses only TV on ordinal distributions—not abs-z.

**First-alarm time** (if event timing is defined): same min-time construction as §1.2 with \(A_{\mathrm{SDD}}\) or \(A_{\mathrm{SDD}}^{\mathrm{sust}}\).

## 2.4 Parameters and basal estimation

| Symbol | Role | Suggested initial range | Notes |
|--------|------|-------------------------|--------|
| \(L_c\) | Current window length | \(30\)–\(100\) symbols (try \(50\)) | Larger than OPC’s \(L\): need stable \(P_t\) |
| \(L_b\) | Basal length | \(\ge 3 L_c\), ideally \(\ge 200\) symbols | Better estimate of \(P_{\mathrm{basal}}\) |
| \(\theta_{\mathrm{TV}}\) | TV alarm threshold | \(0.25\)–\(0.45\) (try \(0.35\)) | On \([0,1]\); tune later only in a dedicated study |
| \(\theta_S\) | Optional sustainment | \(1\) (default) or \(2\)–\(3\) | Not required for formal definition |
| \(\varepsilon\) | KL smoothing (secondary only) | \(1/L_b\) | Unused if TV-only pipeline |
| Basal policy | Fixed vs sliding | **Fixed** primary | Matches honest “baseline repertoire” |

**How to estimate \(P_{\mathrm{basal}}\) (operational):**

1. Map the continuous basal time window used in CCTP (e.g. early Holter hours) to symbol indices after embedding.
2. Count frequencies of \(\sigma_u\) on that index set; normalize by \(L_b\).
3. Freeze \(P_{\mathrm{basal}}\) for the entire pre-event search (Fixed policy).
4. Never recompute basal from the search region (no look-ahead leakage).

## 2.5 Algorithm (pseudocode)

```
Algorithm SDD_DETECT(σ[0..T-1], B_start, B_end, L_c, θ_TV, K, θ_S=1):
  Input:  symbol stream σ; basal index range [B_start, B_end]
  Output: alarm[0..T-1], TV[·]

  # --- basal distribution (fixed) ---
  counts_b ← zeros(K)
  for u ← B_start to B_end:
    counts_b[σ[u]] ← counts_b[σ[u]] + 1
  L_b ← B_end - B_start + 1
  P_basal ← counts_b / L_b

  alarm ← zeros(T); TV ← NaN(T)
  run ← 0

  for t ← 0 to T-1:
    if t < L_c - 1:
      continue
    # skip or mask t inside basal if required by protocol
    counts ← zeros(K)
    for u ← t-L_c+1 to t:
      counts[σ[u]] ← counts[σ[u]] + 1
    P ← counts / L_c
    TV[t] ← 0.5 * sum_s |P[s] - P_basal[s]|

    if TV[t] ≥ θ_TV:
      run ← run + 1
    else:
      run ← 0

    if run ≥ θ_S:
      alarm[t] ← 1

  return alarm, TV, P_basal
```

**Complexity.** \(O(T L_c + L_b)\) naïve; \(O(T+K)\) with sliding frequency updates.

**Reference implementation (optional, non-production):**  
`code/ordinal_detectors/sdd_detector.py` — pure functions; **not** coupled to OPC; **not** wired into `detect_lead_time`.

## 2.6 Advantages vs abs-z ≥ 2 (sustained)

| Aspect | abs-z on continuous \(m(t)\) | SDD (TV) |
|--------|----------------------------|----------|
| Ontology | Moments of aggregated real series | Law of ordinal patterns on \(\Sigma\) |
| What is “abnormal” | Amplitude excursion of \(m\) | Change of **which patterns** are used |
| Gradual drift | Needs large \(|z|\) | Detects redistribution even if a scalar score stays near basal mean |
| Scale invariance | No (depends on \(\sigma_{\mathrm{basal}}\)) | Yes (distributions of ranks/patterns) |
| Device / gain mismatch | Can inflate FAR when \(m\) variance changes | Orthogonal to continuous gain if encoding is rank-based |

## 2.7 Limitations and robustness

1. **Basal quality.** Contaminated basal (early pathology, artifacts) poisons \(P_{\mathrm{basal}}\). Fixed basal is only as honest as the basal window.
2. **Sample size.** Small \(L_c\) → high variance of \(P_t\) → spurious TV spikes; mitigate with larger \(L_c\) or \(\theta_S>1\).
3. **Alphabet size.** Large \(K\) (e.g. full multivariate) needs longer windows for reliable frequencies.
4. **Not collapse-specific.** A redistribution that **increases** diversity can alarm SDD while OPC stays silent (and vice versa)—expected, because the options are separate.
5. **KL misuse.** Using unsmoothed KL as primary would be fragile; hence TV is recommended.
6. **No empirical superiority claim** vs abs-z under this goal.

## 2.8 Robustness notes (ordinal-native)

- TV is Lipschitz in \(\ell_1\) of count vectors: single-symbol flips change TV by at most \(1/L_c\).
- Invariance: any strictly increasing transform of raw amplitudes leaves Bandt–Pompe symbols unchanged → SDD unchanged.
- No dependence on \(\mu,\sigma\) of raw ECG/RR or of continuous \(\tau_s\).

**Option 2 is complete as a standalone detector. No coupling to Option 1 is defined.**

---

# 3. Comparative summary (no fusion)

Both options replace the **ontological role** of abs-z as an alarm layer. Neither is declared the production winner here; neither is combined.

| Criterion | Option 1 — OPC (Persistence Collapse) | Option 2 — SDD (Symbolic Divergence) |
|-----------|----------------------------------------|--------------------------------------|
| Core idea | Low ordinal diversity **and** sustained persistence | Current pattern law diverges from basal law |
| Primary statistic | \(D_t\), \(R_t\) | \(\mathrm{TV}(P_t,P_{\mathrm{basal}})\) |
| Sensitivity profile | Abrupt **collapse / trapping** | Gradual or abrupt **reorganization** |
| RECD alignment | Very high (\(\Phi_2\) persistence, path to \(\Phi_3\)) | High (structural change of ordinal measure) |
| Basal moments of continuous series | None | None |
| Basal object | Thresholds only (no basal distribution required) | Requires \(P_{\mathrm{basal}}\) |
| Interpretability | “Trapped in few ordinal states too long” | “Pattern repertoire reorganized vs basal” |
| Robustness to scale/gain | High (ordinal) | High (ordinal) |
| Main weakness | Misses pure redistribution without collapse | Needs good basal; weaker “collapse” story |
| Implementation complexity | Low–medium | Medium (basal estimation) |
| Relation to frozen abs-z | Conceptual alternative only | Conceptual alternative only |

### Explicit non-goals of this document

- No rule of the form \(A_{\mathrm{OPC}}\wedge A_{\mathrm{SDD}}\) or \(A_{\mathrm{OPC}}\vee A_{\mathrm{SDD}}\).
- No replacement of frozen CCTP abs-z parameters in discovery/validation tables.
- No claim that either option reduces NSRDB FAR or improves VFDB sensitivity without a future dedicated study.

### Suggested manuscript placement (later)

A methods subsection **“Toward native ordinal detectors”** can lift §1 and §2 almost verbatim as two independent proposals, with §3 as a short comparison, and the production abs-z rule retained as the **reported empirical detector** until a bake-off exists.

---

## 4. Consistency with project decisions

| Source | Decision reflected here |
|--------|-------------------------|
| `Contexto.docx` | Two separate options first; no fusion yet; OPC diversity 0.25–0.45 and persistence 4–6; SDD prefers TV over KL |
| CCTP frozen pipeline | abs-z remains production baseline; not retuned |
| RECD (`recd_ordinal_levels.py`) | Symbol stream + persistence language aligned with \(\Phi_2\); no change to \(\Phi_1,\Phi_2,\Phi_3\) definitions |
| Ontological critique | Detectors avoid \(\mu,\sigma\) of original continuous signal and of continuous metric as the alarm trigger |

---

## 5. Reference code map (optional)

| Path | Role |
|------|------|
| `docs/ORDINAL_ALARM_DETECTORS.md` | This formalization (primary deliverable) |
| `code/ordinal_detectors/opc_detector.py` | Option 1 pure reference |
| `code/ordinal_detectors/sdd_detector.py` | Option 2 pure reference |
| `code/ordinal_detectors/__init__.py` | Exports both **separately** (no fused API) |
| `tests/test_ordinal_detectors.py` | Synthetic symbol-stream unit tests per option |

Production entry points (`detect_lead_time`, Phase 1/2 FAR scripts) are **unchanged**.
