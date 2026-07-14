# Mode-C / Mode-S specificity refinement (OPC-G, OPC-BR, OPS)

**Status:** Analysis + design (not a production retune; not a full multi-cohort bake-off)  
**Date:** 2026-07-13  
**Root:** Cardiac CCTP pilot  
**Builds on:** `docs/OPC_RECD_LEVEL3_REFINEMENT_PROPOSALS.md` (proposals A–F), Contexto.docx (separate-first OPC/SDD), `code/ordinal_detectors/opc_refinements.py`  
**Companion baseline:** OPC \((L,\theta_D,\theta_R)=(50,0.35,5)\), sens \(\approx 0.424\) (14/33), FAR \(\approx 3.733\)/24h, decision `keep_baseline`  
**Does not retune:** frozen abs-z; production OPC companion parameters

---

## 0. Goal restatement and the honest tension

**Objective.** Improve **specificity** of the ordinal family so false alarms fall, **without** sacrificing the RECD ability to represent **high-synergy (Nivel 3 / excess3)** reorganizations.

**Hard tension.** Companion OPC is already the specificity-leaning arm (FAR \(\approx 3.73\)/24h vs abs-z \(\approx 33.7\), SDD \(\approx 46.3\)). Its low sensitivity (\(\approx 0.42\)) is not mainly a threshold-tuning failure: the absolute collapse grid returned **`keep_baseline`**—no cell raised all-event sens while keeping FAR in exploratory slack. Further pure collapse tightening can lower FAR further but **deepens Level-3 blindness**.

Therefore:

| What we can do | What we cannot do |
|----------------|-------------------|
| Make **Mode C** (collapse) more specific via ordinal gates that cut *non-informative* collapses | Make Mode C “see” high-synergy Nivel 3 by tightening \(\theta_D\) / \(L\) / \(\theta_R\) |
| Give **Mode S** (ordinal synergy surplus / OPS) specificity so it is usable without FAR explosion | Treat free \(A_C \vee A_S\) as the primary FAR solution |
| Score arms **separately first** (Contexto.docx spirit) | Claim OPS ≡ continuous excess3, or rebrand abs-z(excess3) as ordinal |

**Central honesty statement (non-negotiable):**  
Mode C (OPC / OPC-G / OPC-BR) remains **collapse-based and Level-3–blind** unless Mode S is present in the family. OPC-G and OPC-BR can improve Mode-C FAR without *worsening* L3 blindness relative to absolute-only tightening, but they **do not restore** L3 visibility. Level-3 capacity lives in **Mode S** (OPS as Level-3–consistent **proxy**, not excess3 itself).

---

## 1. Empirical anchors (existing public arms only)

Sources: `docs/ORDINAL_SENSITIVITY_SPECIFICITY_TRADEOFF.md`, `results/ordinal_opc_param_explore_report.md`.

| Arm | All-event sens | FAR / 24 h (NSRDB) | Role |
|-----|----------------|--------------------|------|
| OPC \(L=50\) companion | \(\approx 0.424\) (14/33) | \(\approx 3.733\) | Specificity-leaning collapse |
| SDD (TV) | \(\approx 0.970\) | \(\approx 46.267\) | Sensitivity-leaning divergence |
| abs-z \(\tau_s\) (frozen) | \(\approx 0.909\) | \(\approx 33.734\) | Production primary |

Grid qualitative margins (collapse-only):

- Lower \(\theta_D\), larger \(L\): mean FAR ↓, mean sens ↓ (deeper L3 blindness).  
- Higher \(\theta_D\), smaller \(L\): sens ↑ only with large FAR inflation.  
- Decision: **`keep_baseline`**.

No new multi-cohort FAR/sens numbers are claimed for OPC-G / OPC-BR / OPS below. Directions are conceptual, supported by (i) structural argument, (ii) unit tests on shipped primitives, (iii) a **synthetic stream probe** (not Holter bake-off)—see §8 and scratch log.

---

## 2. Mode C: how OPC-G and OPC-BR can reduce FAR without *deepening* L3 blindness

### 2.1 What non-informative Mode-C false alarms look like

Collapse+persist fires when \(D_t \le \theta_D\) for \(\ge \theta_R\) consecutive windows. False (or low-value) Mode-C alarms tend to be:

1. **Basal-already-collapsed regimes** — paced segments, night-time low variability, sparse alphabets: absolute \(\theta_D\) reads “collapse” that is basal normality.  
2. **Chatter / interrupted noise** — brief low-div flickers that hard \(R_t\) may still accumulate under loose settings, or (with careless \(G\)) gap policy that glues noise into credit.  
3. **Absolute threshold miscalibration on rich controls** — rare sustained low-\(D\) on NSRDB (already only 28 episodes at baseline; further absolute tightening is possible but L3-hostile).

High-synergy Nivel 3 paths that **never** collapse support are **not** Mode-C false alarms—they are Mode-C **structural misses**. Killing them is not a specificity win; it is L3 violation if the *family* has no other arm.

### 2.2 OPC-BR (Proposal C): primary Mode-C FAR lever

**Mechanism (shipped):**  
\[
\mathrm{rel\_cap}=\min(\theta_D,\;\rho\cdot D_{\mathrm{basal}}),\qquad
\ell_t=\mathbf{1}\{D_t\le \mathrm{rel\_cap}\},
\]  
then hard persistence \(\theta_R\). Basal is an early ordinal support reference (same spirit as SDD basal), **not** continuous \(\mu/\sigma\).

**Why FAR can fall without deepening L3 blindness:**

| Effect | Explanation |
|--------|-------------|
| Suppresses absolute false collapses | When \(D_{\mathrm{basal}}\) is already low, \(\rho D_{\mathrm{basal}} < \theta_D\) raises the bar: only *further* collapse alarms. Unit test: low-basal same-alphabet stream → absolute OPC can low-div; BR with \(\rho=0.85\) stays quiet. |
| Preserves true relative collapses | Rich basal \(\to\) deep drop still meets both absolute and relative caps. |
| Does not change positive class ontology | Still collapse. Still silent on high-\(D\) high-synergy. **Does not claim L3 restore.** |
| Better than absolute \(\theta_D\downarrow\) alone | Absolute \(\theta_D=0.30\) on the grid crushed sens (mean sens\(\approx 0.18\)) and zeros many true mild collapses *and* any near-threshold path. BR is **state-dependent**: only tight when basal is already poor. |

**Honesty.** OPC-BR **remains Level-3–blind**. It improves Mode-C *specificity* among collapse-like phenomena. It does **not** open high-synergy paths. Relative to absolute-only tightening, it is **less L3-hostile** because it need not globally lower \(\theta_D\) on rich-basal event Holters.

**Recommended starting configs (Mode C primary candidates):**

| ID | Parameters | Expected FAR direction | Expected sens direction | L3 status |
|----|------------|------------------------|-------------------------|-----------|
| **C0** companion | \(L=50,\theta_D=0.35,\theta_R=5\) | baseline \(\approx 3.73\) | baseline \(\approx 0.42\) | Blind |
| **C-BR1** | companion + \(\rho=0.85\), basal_end = Phase-2 early basal | **↓** on low-basal controls; ~ flat on rich | ~ or slight ↓ if some “collapses” were basal-normal | Blind (not worse than C0) |
| **C-BR2** | companion + \(\rho=0.75\), \(\theta_R=5\) | **↓↓** stricter relative | risk mild sens ↓ on borderline relative drops | Blind |
| **C-BR3** | companion + \(\rho=0.85\), \(\theta_R=6\) | **↓** (relative + longer run) | slight ↓ | Blind |
| **C-BR4** (stack) | \(\theta_D=0.35\), \(\rho=0.80\), \(\theta_R=6\) | **↓** dual gate | slight ↓ | Blind — still prefer over absolute \(\theta_D=0.30\) alone |

**Registration order:** freeze C0 as reference; evaluate C-BR1 first (minimal change); only then C-BR2/3. Do **not** pair BR with free \(\theta_D\uparrow\).

### 2.3 OPC-G (Proposal A): intermittency polish, not a free specificity upgrade

**Mechanism (shipped):** persistence **credit** increments on low-div; up to \(G\) consecutive high-div windows hold credit; more than \(G\) resets. Alarm only when current window is low-div **and** credit \(\ge \theta_R\).

**Why G alone is not a FAR reducer:**

- \(G>0\) **recovers** interrupted true collapses → tends to **↑ sens and ↑ FAR risk**.  
- Unit test: G=0 matches hard OPC; G>0 ≥ hard alarm mass on interrupted streams.  
- Marketing “gaps = more specific” is **false**.

**How G can participate in a net FAR-aware Mode-C design:**

| Co-constraint | Intent |
|---------------|--------|
| Raise \(\theta_R\) (e.g. 5 → 7–8) with \(G=1\) | Require more low-div mass; gaps only glue genuine intermittent locks, not single flickers |
| Slightly lower \(\theta_D\) (e.g. 0.32–0.30) with \(G=1\) | Stricter collapse definition; G recovers intermittent *true* deep locks only |
| Prefer after BR, not instead of BR | BR cuts basal false collapses; G addresses limit (2) intermittency without continuous moments |

**Recommended starting configs:**

| ID | Parameters | Expected FAR | Expected sens | L3 status |
|----|------------|--------------|---------------|-----------|
| **C-G0** | \(G=0\) ≡ hard OPC | = C0 | = C0 | Blind |
| **C-G1** | \(G=1\), same \((L,\theta_D,\theta_R)\) | **↑ or flat** (risk) | **↑** intermittent events | Blind — **not** for FAR-first |
| **C-G2** | \(G=1\), \(\theta_R=8\), \(\theta_D=0.35\) | target **≤ C0** if co-tuned | recover some intermittent vs hard \(\theta_R=8\) alone | Blind |
| **C-G3** | \(G=1\), \(\theta_D=0.30\), \(\theta_R=6\) | **↓** vs C0 likely | **↓** vs C0 (grid warns) | Blind — **deeper** absolute blindness risk |
| **C-GBR** | \(G=1\), \(\rho=0.85\), \(\theta_R=7\) | **↓** if BR dominates false collapses | recover intermittent true relative collapses | Blind |

**Reject:** C-G1 as a “specificity improvement.”  
**Accept:** C-G2 / C-GBR only under pre-registered co-constraints, secondary to BR.

### 2.4 Mode-C summary: FAR without *aggravating* L3

```
Better Mode-C FAR path:     OPC-BR (± longer θ_R)  →  optional G with raised θ_R
Worse Mode-C path for L3:   absolute θ_D↓ / L↑ alone expecting L3 rescue (Proposal E, rejected)
Impossible Mode-C claim:    “BR or G restores Nivel 3”
```

Mode-C refinements that **respect** RECD L3 are those that:

1. Do not force every pre-VF path through collapse; and  
2. Do not globally annihilate mild-diversity event regimes solely to chase FAR, **when** Mode S is absent—or clearly label that Mode C alone remains L3-incomplete.

---

## 3. Mode S: OPS + ordinal specificity so high-synergy is usable

### 3.1 What OPS measures (shipped)

On bivariate factors \((\pi^{(1)},\pi^{(2)})\):

\[
S_t = \mathrm{TV}\big(P_t(\pi_1,\pi_2),\; P_t(\pi_1)\otimes P_t(\pi_2)\big).
\]

Alarm via absolute \(S_t \ge \theta_S\) or basal-relative \(S_t - S_{\mathrm{basal}} \ge \theta_{\Delta S}\), with consecutive persistence \(\theta_R^S\).

**Labeling rule:** OPS is a **Level-3–consistent discrete proxy** (joint surplus beyond independent margins). It is **not** continuous excess3, not nested RECD \(\theta_3\) machinery, and not a drop-in scientific replacement for excess3 readouts.

### 3.2 Why absolute OPS will explode FAR

Empirical joint TV vs product margins fluctuates on ambulatory streams (finite-\(L\) sampling alone yields non-zero \(S_t\) under independence). Absolute \(\theta_S\) near the noise floor produces chatter; too high and true \(\Delta S\) events are missed. **Specificity must be ordinal-only:**

| Lever | Role | Continuous \(z\)? |
|-------|------|-------------------|
| Basal-relative \(\Delta S\) | Alarm only on **increase** over early basal mean \(S\) | No |
| Longer \(\theta_R^S\) (e.g. 8–12) | Sustained surplus, not blips | No |
| Refractory episode counting (protocol 0.5 h) | FAR definition already; may also gate detector re-arm | No |
| Optional gap \(G_S\) | Intermittent high-synergy locks (secondary; same co-tune discipline as OPC-G) | No |
| Optional dual discrete floor | e.g. require TV vs basal *symbol law* above a floor (SDD-like) — labeled fusion-adjacent | No |

### 3.3 Recommended OPS configs (Mode S primary candidates)

| ID | Parameters | Expected FAR | Expected sens (synergy events) | L3 status |
|----|------------|--------------|--------------------------------|-----------|
| **S0** | \(L_S=50\), \(\theta_S=0.15\), \(\theta_R^S=5\) absolute | **↑↑ risk** | high on raw dependence | L3-consistent proxy; **not production-ready** |
| **S1** | \(L_S=50\), \(\theta_{\Delta S}=0.08\), \(\theta_R^S=5\), basal_end = Phase-2 | controlled vs S0 | fires on sustained \(\Delta S\) | Preferred starting Mode S |
| **S2** | \(\theta_{\Delta S}=0.12\), \(\theta_R^S=8\) | **↓** vs S1 | slight ↓; fewer blips | Specificity-leaning Mode S |
| **S3** | \(\theta_{\Delta S}=0.15\), \(\theta_R^S=10\) | **↓↓** | may miss mild surplus rises | Strict Mode S |
| **S-G1** | S2 + \(G_S=1\) (if gap API extended) | co-tune only | intermittent synergy | Secondary |

**Synthetic probe (not Holter bake-off):** on full-support high-synergy streams (joint mass on ridges + floor so \(D_t\) stays moderate/high), \(S\) rose basal→approach while OPC remained a weak/noisy collapse signal; basal-relative OPS with \(\theta_{\Delta S}\in\{0.12,0.15\}\) and longer \(\theta_R\) zeroed independent-control episodes in the toy generator while still locking the high-syn block. **Do not cite as cohort FAR.**

### 3.4 Mode S does not replace Mode C

- Collapse locks (few active joints) can be low-\(S\) or high-\(S\); OPS is not a universal event detector.  
- Diagonal perfect coupling can *look* like collapse in joint alphabet (\(K=36\), only 6 joints)—OPS and OPC may co-fire; that is **not** proof OPC sees L3 in general.  
- True L3 stress case: **high joint support + high surplus** (many coordinated joints). Mode C silent or fragile; Mode S designed for that axis.

---

## 4. Separate-arm evaluation strategy (primary before combination)

Aligned with Contexto.docx (“trabajar Opción 1 y Opción 2 por separado primero”) and proposals doc §2–§3D.

### 4.1 Primary arms (each owns sens + FAR)

| Arm ID | Detector | Positive class | Primary metrics |
|--------|----------|----------------|-----------------|
| **Mode C** | OPC companion ± BR ± G | Sustained ordinal diversity collapse | sens (SDDB/VFDB/all), FAR NSRDB episode×0.5h |
| **Mode S** | OPS (± basal-relative, long persist) | Sustained ordinal synergy surplus | same tables, **separate rows** |
| **Reference** | frozen abs-z \(\tau_s\) | Continuous pilot primary | reference only; **not retuned** |
| **Context** | SDD (TV) | Distributional divergence | optional context row; not fused into Mode C/S primary |

### 4.2 Protocol (reuse existing public protocol)

1. **Encoding:** joint Bandt–Pompe bivariate as today for Mode C (\(K=36\)); explicit factor pair \((\pi^{(1)},\pi^{(2)})\) for Mode S (same embedding factorization).  
2. **Sensitivity:** first alarm in post-basal pre-event window; SDDB + VFDB; report by cohort and pooled.  
3. **FAR:** NSRDB controls; binary episodes; refractory **0.5 h**; Phase-2 basal/search/control cap conventions.  
4. **Basal_end:** identical Phase-2 early basal definition for BR and OPS \(\Delta S\).  
5. **No free retune** of abs-z; no clinical / FDA / S5 claim language.

### 4.3 Reporting tables (required shape)

**Table P — Primary separate arms**

| Arm | Config ID | sens_all | FAR/24h | notes |
|-----|-----------|----------|---------|-------|
| Mode C | C0, C-BR1, … | … | … | collapse-only; L3-blind |
| Mode S | S1, S2, … | … | … | L3-proxy; not excess3 |
| abs-z | frozen | … | … | reference |

**Table S — Secondary combination (labeled only)**

| Rule | Purpose | Must show |
|------|---------|-----------|
| \(A_C \vee A_S\) | ordinal coverage upper bound | separate + union FAR (expect ↑) |
| \(A_S\) confirm \(A_C\) within \(\pm w\) | cascade analogue | gain/loss like SDD→OPC |
| Free AND | high precision exploratory | usually low sens |

**Forbidden as primary score:** union FAR presented as “the” improved detector without Mode C and Mode S rows.

### 4.4 Pre-registration sketch (before any bake-off)

1. Freeze **C0** companion.  
2. Mode C grid (small): {C-BR1, C-BR3, optional C-GBR} — success = FAR ≤ C0 with sens not collapsing below a pre-set floor *or* FAR↓ with only modest sens↓, **and** explicit L3-blind label retained.  
3. Mode S grid (small): {S1, S2, S3} — success = non-saturating control FAR **and** qualitative capture of known high-synergy synthetic + any Holter subset tagged excess3-up / high \(\Delta\)excess3 (diagnostic continuous readout only).  
4. Secondary union/cascade only after P tables exist.  
5. Reject cells that reintroduce \(\mu/\sigma\) or abs-z on excess3.

---

## 5. Concrete parameter proposals (summary sheet)

### 5.1 Mode C — FAR-first stack

| Priority | Config | \(L\) | \(\theta_D\) | \(\theta_R\) | \(G\) | \(\rho\) | Intent |
|----------|--------|------|--------------|--------------|-------|---------|--------|
| 1 | **C-BR1** | 50 | 0.35 | 5 | 0 | 0.85 | Primary FAR cut on low-basal false collapses |
| 2 | **C-BR3** | 50 | 0.35 | 6 | 0 | 0.85 | BR + slightly longer run |
| 3 | **C-GBR** | 50 | 0.35 | 7 | 1 | 0.85 | Intermittency after BR; co-tuned θ_R |
| — | C-G1 alone | 50 | 0.35 | 5 | 1 | — | **Not recommended** for FAR-first |
| — | Absolute \(\theta_D=0.30\) only | 50 | 0.30 | 5 | 0 | — | Grid-like FAR↓ / sens↓; **L3-hostile**; reject as L3 answer |

### 5.2 Mode S — usable L3 proxy

| Priority | Config | \(L_S\) | gate | \(\theta_R^S\) | Intent |
|----------|--------|---------|------|----------------|--------|
| 1 | **S1** | 50 | \(\theta_{\Delta S}=0.08\) | 5 | Starting Mode S |
| 2 | **S2** | 50 | \(\theta_{\Delta S}=0.12\) | 8 | Specificity-leaning Mode S |
| 3 | **S3** | 50 | \(\theta_{\Delta S}=0.15\) | 10 | Strict Mode S |
| — | S0 absolute \(\theta_S=0.15\) | 50 | absolute | 5 | Diagnostic only; FAR risk |

### 5.3 What not to propose as core

- Continuous basal \(\mu/\sigma\) or \(z\)-score on \(D_t\), \(S_t\), or excess3 as OPC/OPS core.  
- Free OR as primary FAR narrative.  
- Another absolute collapse-only grid as Level-3 fix (Proposal E).  
- abs-z(excess3) rebranded ordinal (Proposal F).

---

## 6. Conceptual trade-off assessment

| Design move | Sens (concept) | FAR (concept) | L3 family capacity | Verdict |
|-------------|----------------|---------------|--------------------|---------|
| Stricter absolute OPC only | ↓ | ↓ | **worse** (deeper blindness) | FAR polish only after Mode S exists; not L3 solution |
| OPC-BR alone | ~ / slight ↓ | **↓** on low-basal controls | unchanged (still blind Mode C) | **Best Mode-C FAR lever** |
| OPC-G alone | ↑ intermittent | **↑ risk** | unchanged | Not FAR-first |
| OPC-G + raised \(\theta_R\) / BR | recover intermittent true locks | target ≤ baseline | unchanged | Conditional Mode-C polish |
| OPS absolute | ↑ synergy paths | **↑↑** | opens Mode S | Unusable without gates |
| OPS \(\Delta S\) + long persist | ↑ true \(\Delta S\) | controlled | **opens L3-proxy** | **Main Mode-S path** |
| Dual-mode separate report | per-arm | per-arm | **yes via Mode S** | **Required architecture** |
| Free OR primary | ↑↑ | ↑↑ | yes but muddled FAR | Secondary only |

**Honest summary:**  
You cannot make collapse-only OPC both more specific *and* Level-3–aware. You **can** (i) cut non-informative Mode-C alarms with **OPC-BR** (± disciplined OPC-G), and (ii) place Level-3–consistent detection in **Mode S / OPS** with basal-relative surplus and long persistence, scored separately. That is how the family lowers FAR *where collapse is the wrong positive class* without structurally excluding high-synergy reorganizations.

---

## 7. Ranked recommendation

### Priority 1 — Dual-mode family with separate scoring (Proposal D)

Architecture: Mode C and Mode S as **primary separate arms**; abs-z frozen reference; free OR/AND/cascade only secondary labeled tables.

### Priority 2 — Mode-C FAR: OPC-BR first (Proposal C)

Adopt exploration order **C-BR1 → C-BR3 → optional C-GBR**. This is the most coherent pure-ordinal way to improve Mode-C specificity **without** the L3 hostility of absolute \(\theta_D\downarrow\) grids.

### Priority 3 — Mode-S usability: OPS + \(\Delta S\) + long \(\theta_R^S\) (Proposal B)

Order **S1 → S2 → S3**. This is the only path in-scope that **opens** high-synergy visibility inside ordinal alarm space.

### Priority 4 — OPC-G (Proposal A) only as co-tuned polish

Never alone as FAR upgrade; only with raised \(\theta_R\) and preferably after BR.

### Rejected as answers to this goal

| ID | Why |
|----|-----|
| E absolute collapse-only L3 fix | Grid + ontology fail |
| F abs-z(excess3) as OPC/OPS | Breaks ordinal-native constraint |
| Free OR as primary FAR solution | Multiplicity cost; hides per-arm truth |
| C-G1 alone marketed as “more specific” | False |

### Most coherent tool combination

**OPC-BR (Mode C specificity) + OPS basal-relative (Mode S L3 proxy) + separate-arm evaluation**, with **OPC-G optional** under co-constraints. Gap tolerance alone does not solve FAR or L3; basal-relative collapse alone does not solve L3; OPS without \(\Delta S\)/persistence does not solve FAR.

---

## 8. Level-3 respect criteria (pass / fail)

Use these as **go/no-go tests** for any proposed “improvement.”

### 8.1 Pass (respects RECD Nivel 3)

| ID | Criterion |
|----|-----------|
| **L3-P1** | High-synergy / high-excess3-style paths remain **representable** in the *family* (Mode S or an explicitly labeled continuous excess3 arm)—not required to pass through collapse. |
| **L3-P2** | Mode C changes that lower FAR do **not** claim to restore L3; documentation states Mode C remains collapse-based and L3-blind. |
| **L3-P3** | OPS (if used) is labeled **Level-3–consistent proxy**, not excess3; no identity claim. |
| **L3-P4** | No abs-z / \(\mu/\sigma\) on continuous excess3 (or on \(D_t\), \(S_t\)) as proposed **core** of OPC/OPS. |
| **L3-P5** | Primary evaluation reports Mode C and Mode S **separately**; combination rules are secondary and labeled. |
| **L3-P6** | Stress case acknowledged: high joint support + high surplus may silence Mode C; family still has Mode S. |

### 8.2 Fail (violates or undermines Nivel 3)

| ID | Criterion |
|----|-----------|
| **L3-F1** | Collapse-only retune marketed as solving high-synergy pre-VF paths. |
| **L3-F2** | Absolute \(\theta_D\downarrow\) / \(L\uparrow\) sold as “better RECD alignment” without Mode S. |
| **L3-F3** | Free OR presented as primary detector so Mode-C FAR “inherits” Mode-S hits without separate FAR tables. |
| **L3-F4** | Continuous \(z\)-score on excess3 rebranded as native ordinal OPC/OPS. |
| **L3-F5** | OPS claimed identical to excess3 / nested \(\Phi_3\) science readout. |
| **L3-F6** | OPC-G alone claimed to fix L3 blindness or to improve specificity without co-tuning disclosure. |

### 8.3 Quick self-check for any new config cell

1. Does this cell **only** fire on low \(D_t\)? → Mode C; L3-blind; OK if labeled.  
2. Can a full-support high-\(S\) stream alarm **some** primary arm? → need Mode S (or fail L3-P1).  
3. Did FAR improve only by silencing mild-diversity **and** high-syn paths with no Mode S? → **L3-F2**.  
4. Any \(\mu/\sigma\)? → **L3-F4**.

---

## 9. Primitives map (implementation already in repo)

| Proposal | Function | File |
|----------|----------|------|
| A OPC-G | `opc_detect_gap_tolerant` | `code/ordinal_detectors/opc_refinements.py` |
| C OPC-BR | `opc_detect_basal_relative` | same |
| B OPS surplus | `ordinal_synergy_surplus` | same |
| B OPS detect | `ops_detect` | same |
| Joint encode | `joint_symbols_from_factors` | same |
| Exports | `__init__.py` | `code/ordinal_detectors/__init__.py` |
| Unit tests | `tests/test_opc_refinements.py` | gap vs hard, BR suppresses low-basal, OPS dep vs ind, basal-relative \(\Delta S\) |

Exploratory only — not production ship.

---

## 10. Relation to Contexto.docx

Contexto.docx formalized:

- **Opción 1:** Ordinal Persistence Collapse (OPC) — diversity + run-length.  
- **Opción 2:** Symbolic Distribution Divergence (SDD) — TV/KL vs basal.  
- Method: **work separately first**; fusion only if neither is sufficient.

This note **preserves** separate-first discipline and extends it:

- Mode C = refined Opción 1 (BR, optional G) for **specificity**.  
- Mode S = new ordinal arm for **Nivel-3–consistent surplus** (not SDD; SDD remains distributional shift).  
- SDD stays the sensitivity-leaning **divergence** arm; not merged into Mode C/S primary scores here.  
- Aspirational “collapse → \(\Phi_3\)” language in early formalization is treated as **over-promise**: collapse is Φ₂-adjacent locking, not Level-3 surplus.

---

## 11. Bottom line (deliverable checklist)

1. **Mode-C FAR without deepening L3 blindness:** prefer **OPC-BR** (relative collapse vs basal \(D\)); use **OPC-G only co-tuned** (raised \(\theta_R\) / after BR). Both remain L3-blind—stated explicitly. Absolute collapse-only grids are the wrong L3 answer.  
2. **Mode-S usability:** OPS with **basal-relative \(\Delta S\)**, longer \(\theta_R^S\), refractory episode protocol; absolute OPS is diagnostic-only.  
3. **Evaluation:** Mode C vs Mode S **separate primary tables**; OR/AND/cascade secondary labeled only.  
4. **Configs:** C-BR1/C-BR3/C-GBR and S1/S2/S3 as concrete starting points (§5).  
5. **Recommendation:** **D + C + B**, then optional **A**; reject E/F and free-OR-primary.  
6. **L3 criteria:** §8 pass/fail tests—use before claiming any refinement “respects RECD Nivel 3.”

---

## 12. Synthetic probe note (exploratory, non-cohort)

A small synthetic generator (rich control, low-basal control, true collapse, intermittent collapse, full-support high-synergy) was used to **check directional behavior** of shipped primitives:

- OPC-BR zeroed low-basal false collapses while keeping true collapse.  
- OPC-G alone did not reduce intermittent chatter; needs co-tuning.  
- Full-support high-synergy raised \(S\) basal→approach; basal-relative OPS locked the approach block; independent streams quieted under stricter \(\Delta S\) + \(\theta_R^S\).  
- Perfect diagonal coupling is a **confound** (looks like collapse in joint \(K\)); L3 stress cases must use **rich support + high surplus**.

These are **not** replacements for NSRDB/SDDB/VFDB bake-off numbers. Full empirical FAR/sens for C-BR* / S* remains future work (deferred bake-off).
