# Integrated ordinal detector: Nivel-2 persistence on Nivel-3 synergy surplus

**Status:** Design + shipped exploratory primitive (not a production retune; not a multi-cohort Holter bake-off)  
**Date:** 2026-07-13  
**Root:** Cardiac CCTP pilot  
**Builds on:** `docs/OPC_MODE_C_S_SPECIFICITY_REFINEMENT.md`, `docs/OPC_RECD_LEVEL3_REFINEMENT_PROPOSALS.md` (Proposal B OPS), Contexto.docx (Φ₂ persistence + ordinal space), `code/ordinal_detectors/opc_refinements.py`  
**Companion anchors (existing public arms only):** OPC sens \(\approx 0.424\), FAR \(\approx 3.733\)/24h; SDD \(0.970\) / \(46.267\); abs-z \(0.909\) / \(33.734\)  
**Does not retune:** frozen abs-z; production OPC companion \((L,\theta_D,\theta_R)=(50,0.35,5)\)

---

## 0. Goal and the design pivot

**Objective.** Design an **integrated** ordinal alarm whose **primary** mechanism is **sustained persistence** (Nivel-2 / Φ₂ run-length logic) applied to a **high synergistic-surplus** signal (Nivel-3–consistent ordinal proxy), with collapse/locking allowed only as a **complementary** role—not as the dominant criterion.

**Pivot relative to separate Mode C / Mode S.**  
The Mode-C/S note kept collapse (Mode C) and OPS (Mode S) as **separate primary arms**, scored first, with free OR only secondary. That remains the right **evaluation discipline** for fairness tables. This note is a **labeled next design**: one detector whose ontology is already fused as

\[
\text{alarm primary} = \mathrm{persist}\big(\text{high surplus } S_t\text{ or }\Delta S\big),
\]

with collapse \(D_t \le \theta_D\) only as tag / modulate / confirm—**never** as the sole or dominant predicate.

| Design | Primary predicate | Collapse role | L3 high-support path |
|--------|-------------------|---------------|----------------------|
| Collapse-only OPC | \(\mathrm{persist}(D\le\theta_D)\) | Dominant | **Silent** (structural miss) |
| Separate Mode C ∪ Mode S | Two independent alarms | Equal peers | Open via Mode S |
| Free \(A_C\vee A_S\) | Union of peers | Can dominate FAR narrative | Open, muddled FAR |
| **This integrated detector (OPSP)** | \(\mathrm{persist}(\text{high }S/\Delta S)\) | Subordinate only | **Open** on I0 / tag / modulate |

**Central honesty (non-negotiable):**

1. The surplus signal is OPS / \(S_t=\mathrm{TV}(P_{\mathrm{joint}},P_1\otimes P_2)\) — a **Level-3–consistent proxy**, **not** continuous excess3 / nested \(\Phi_3\).  
2. Collapse is **subordinate**. A pure collapse stream with near-zero surplus **must not** alarm.  
3. This is **not** a refined collapse-only OPC.  
4. Free OR of collapse ∪ surplus is **not** “the” integrated story unless surplus-persistence primacy is explicit.

---

## 1. Signals (pure ordinal)

### 1.1 Nivel-3–consistent surplus (primary signal)

On bivariate ordinal factors \((\pi^{(1)},\pi^{(2)})\), window length \(L\):

\[
S_t = \mathrm{TV}\big(P_t(\pi_1,\pi_2),\; P_t(\pi_1)\otimes P_t(\pi_2)\big).
\]

Basal-relative gate (preferred for FAR control; same spirit as Mode-S S1–S3):

\[
\ell^S_t = \mathbf{1}\{S_t - S_{\mathrm{basal}} \ge \theta_{\Delta S}\}
\quad\text{or absolute}\quad
\mathbf{1}\{S_t \ge \theta_S\}.
\]

### 1.2 Nivel-2 persistence logic (applied to surplus, not to collapse)

Run-length on \(\ell^S\):

\[
R^S_t =
\begin{cases}
R^S_{t-1}+1 & \text{if }\ell^S_t=1,\\
0 & \text{otherwise.}
\end{cases}
\]

**Core (I0) alarm:**

\[
A^{\mathrm{core}}_t = \mathbf{1}\{\ell^S_t=1\}\wedge\mathbf{1}\{R^S_t \ge \theta_R^S\}.
\]

This is Φ₂ **persistence** (Contexto.docx Opción 1 *logic*) on a Φ₃-adjacent **surplus** signal (Proposal B), not on diversity collapse.

### 1.3 Collapse / locking (complementary only)

On joint codes \(\sigma_t = \pi^{(1)}_t\cdot k_2 + \pi^{(2)}_t\), \(K=k_1 k_2\):

\[
D_t = \frac{|\mathrm{supp}(W_t)|}{K},\qquad
\ell^C_t = \mathbf{1}\{D_t \le \theta_D\}.
\]

\(\ell^C\) never defines \(A^{\mathrm{core}}\). It only enters the roles in §2.

---

## 2. Concrete integration rules

Shipped entry point: `opsp_integrated_detect` in `code/ordinal_detectors/opc_refinements.py`.

### 2.1 Role catalogue

| Role ID | `collapse_role` | Alarm formula | Collapse may create alarm alone? | Intent |
|---------|-----------------|---------------|----------------------------------|--------|
| **I0** | `none` | \(A = A^{\mathrm{core}}\) | **No** | Preferred starting integrated core |
| **I-tag** | `tag` | \(A = A^{\mathrm{core}}\); flag co-occurrence | **No** | Diagnostics / subtype labels |
| **I-mod** | `modulate` | \(A^{\mathrm{core}}\) with **relaxed** \(\theta_{\Delta S},\theta_R^S\) when \(\ell^C_t=1\) | **No** | Soft assist when locking co-presents |
| **I-confirm** | `confirm` | \(A = A^{\mathrm{core}} \wedge \exists\,s\in[t\pm w]:\ell^C_s=1\) | **No** | FAR-leaning secondary; may silence pure high-support L3 |

**Universal invariant (all roles):**  
\[
A_t=1 \;\Rightarrow\; \text{core surplus-persistence path was active (or modulate-relaxed surplus path).}
\]
Equivalently: collapse-only \(\Rightarrow A\equiv 0\).

### 2.2 Modulate detail (I-mod)

When \(\ell^C_t=1\):

\[
\theta_{\Delta S}^{\mathrm{eff}} = \max\big(\tfrac12\theta_{\Delta S},\;\theta_{\Delta S}-\delta_S\big)
\quad(\text{if }\theta_{\Delta S}>0),\qquad
\theta_R^{\mathrm{eff}} = \max(1,\;\theta_R^S-\delta_R).
\]

When \(\ell^C_t=0\), full \(\theta_{\Delta S},\theta_R^S\) apply. Surplus mass is still required; collapse only **eases** the surplus gates. The half-nominal floor prevents \(\theta^{\mathrm{eff}}\to 0\) from turning collapse-only (near-zero surplus) streams into false high-syn indicators.

Default exploratory: \(\delta_R=2\), \(\delta_S=0.02\).

### 2.3 Confirm detail (I-confirm)

\[
A_t = A^{\mathrm{core}}_t \wedge \mathbf{1}\Big\{\sum_{s=t-w}^{t+w}\ell^C_s \ge 1\Big\}.
\]

Collapse **filters** core alarms; it does not invent them.  
**L3 cost:** high-support high-\(S\) streams with \(\ell^C\equiv 0\) lose alarms. Label I-confirm as **secondary FAR-leaning**, not the RECD-default integrated arm.

### 2.4 Recommended named configs

| ID | Role | \(L\) | Surplus gate | \(\theta_R^S\) | Collapse params | Use |
|----|------|-------|--------------|----------------|-----------------|-----|
| **I0** | none | 50 | \(\theta_{\Delta S}=0.08\) | 5 | — | **Primary recommended** |
| **I0-strict** | none | 50 | \(\theta_{\Delta S}=0.12\) | 8 | — | Specificity-leaning core |
| **I-tag** | tag | 50 | \(\theta_{\Delta S}=0.08\) | 5 | \(w=5\) (samples) | Same alarms as I0 + labels |
| **I-mod** | modulate | 50 | \(\theta_{\Delta S}=0.08\) | 5 | \(\theta_D=0.35\), \(\delta_R=2\), \(\delta_S=0.02\) | Soft N2 assist |
| **I-confirm** | confirm | 50 | \(\theta_{\Delta S}=0.12\) | 8 | \(\theta_D=0.35\), \(w=5\) | Secondary FAR filter |
| **Reject** | — | — | — | — | collapse as primary | Refined OPC; out of scope |

Absolute \(S_t\ge\theta_S\) remains **diagnostic only** (same FAR risk as Mode-S S0).

### 2.5 Explicit non-rules (rejected as “integration”)

| Rejected | Why |
|----------|-----|
| \(A = \mathrm{persist}(D\le\theta_D)\) with surplus as optional tag | Collapse-dominant; reverts to OPC |
| \(A = A_C \vee A_S\) marketed as integrated primary | Free OR; no surplus primacy; muddled FAR |
| \(A = A_C\) then “confirm with surplus” as primary story | Inverts hierarchy (cascade SDD→OPC style on wrong axis) |
| abs-z / \(\mu/\sigma\) on excess3 as core | Breaks ordinal-native constraint |
| OPS ≡ excess3 claim | Proxy only |

---

## 3. Ontological implications (N2 + N3 vs separate arms)

### 3.1 What the integrated detector *is*

- **Nivel 2 contribution:** the *form* of evidence — sustained run-length / locking-in-time (Φ₂ language from Contexto Opción 1).  
- **Nivel 3 contribution:** the *content* of evidence — joint surplus beyond independent margins (Proposal B / OPS).  
- **Nested RECD spirit:** higher-order reorganization is read as **persistent surplus**, not as “few symbols.”  
- **Collapse (Φ₂-adjacent locking of support):** optional co-factor (co-locking), not the definition of the positive class.

Positive class (I0): *sustained increase (or high absolute) of ordinal synergistic surplus*.

### 3.2 What separate arms did better / worse

| Aspect | Separate Mode C / Mode S | Integrated OPSP |
|--------|--------------------------|-----------------|
| Fair FAR tables | **Best** — each arm owns FAR | Need still to report I0 vs I-confirm separately |
| L3 high-support visibility | Open via Mode S | Open via I0 / I-tag / I-mod |
| Collapse-only events | Visible on Mode C | **May miss** if surplus low (honest gap) |
| Single narrative detector | Weaker (two stories) | **Stronger** one primary story |
| Risk of collapse dominance | Low if scored separate | Controlled by role catalogue + invariant |
| Free-OR temptation | Explicit secondary | Explicitly rejected as primary |

### 3.3 Honest gap vs dual-arm family

Integrated OPSP **prioritizes** high-synergy reorganizations (with or without collapse).  
Events that are **pure collapse with little surplus** (some Φ₂ locks, sparse alphabets without joint structure) may be **silent** under I0—exactly because collapse is not dominant.  

If those events matter clinically/scientifically, keep **Mode C as a separate report row** (companion OPC / OPC-BR). Integration does **not** delete Mode C from the family; it defines a **surplus-primary** detector for the N2×N3 story.

### 3.4 Relation to Contexto.docx

Contexto formalized Opción 1 as collapse+persist and asked for separate-first work before fusion.  
This design:

- **Reuses** Opción 1’s *persistence* operator.  
- **Replaces** the Opción 1 *diversity* operand with OPS surplus.  
- Treats aspirational “collapse → \(\Phi_3\)” language as **over-promise**: collapse is not Level-3 detection; persistent surplus is the Level-3–consistent ordinal path.  
- Does **not** violate separate-first *evaluation* discipline: §4 still requires comparison tables against Mode C, Mode S, and free OR as **contrasts**.

---

## 4. Evaluation (conceptual + synthetic on shipped code)

### 4.1 Conceptual sens / FAR directions

Using only **existing public anchors** and structural argument—no invented Holter bake-off numbers for I*.

| Comparator | Existing anchor | Conceptual I0 (surplus-persist) | Conceptual I-confirm | Conceptual free \(A_C\vee A_S\) |
|------------|-----------------|----------------------------------|----------------------|--------------------------------|
| Collapse OPC companion | sens \(\approx 0.42\), FAR \(\approx 3.73\)/24h | Higher sens on **high-syn** events OPC misses; FAR depends on \(\theta_{\Delta S},\theta_R^S\) (risk \(>\) OPC if absolute OPS; controlled if S1-like gates) | FAR **↓** vs I0; sens **↓** on pure high-support L3 | FAR **↑↑** (union) |
| Mode S OPS alone | no Holter FAR published | I0 **is** Mode-S core with optional collapse assist | Stricter than Mode S | Union with Mode C |
| Free OR | not primary | Better ontology (one primary) | Not free OR | — |
| abs-z frozen | sens \(0.909\), FAR \(33.7\) | Not a production replacement claim | — | — |

**Qualitative trade-off sheet:**

| Design move | Sens (concept) | FAR (concept) | L3 high-support | Collapse-only events |
|-------------|----------------|---------------|-----------------|----------------------|
| Collapse OPC only | low on L3 paths | low | **blind** | detected |
| Mode S I0-like alone | up on \(\Delta S\) paths | medium if gated | **open** | often silent |
| **I0 integrated core** | same primary as Mode S | same as Mode S gates | **open** | silent if no surplus |
| **I-mod** | slight ↑ when co-lock | slight ↑ risk | open | still needs surplus |
| **I-confirm** | ↓ vs I0 if no collapse | **↓** vs I0 | **risk of silence** | only if surplus also |
| Free \(A_C\vee A_S\) | high | high | open | open |
| Collapse-dominant “integration” | — | — | **fails RECD** | — |

### 4.2 Synthetic probe (unit-backed, not Holter)

Shipped tests in `tests/test_opc_refinements.py` drive `opsp_integrated_detect` on synthetic streams:

| Probe | Expected (verified by tests) |
|-------|------------------------------|
| Indep basal → sustained dependence | I0 / core alarms; matches `ops_detect` high_syn / alarm |
| Constant single joint (collapse, \(S\approx 0\)) | **No** alarm under none/tag/modulate/confirm |
| Ridge+floor high-support high-surplus | I0 alarms; mean \(D>0.35\); confirm ⊆ core |
| Coupled low-alphabet (surplus + collapse) | modulate ≥ core alarm mass possible |
| Unknown `collapse_role` | `ValueError` |
| Named I0 / I-mod / I-confirm | callable; I0 fires on \(\Delta S\) rise |

**Do not cite synthetic alarm counts as NSRDB FAR or SDDB/VFDB sens.**

### 4.3 What a future Holter bake-off should report (not run here)

1. **Primary rows:** I0, I0-strict (each: sens_all, FAR/24h episode×0.5h).  
2. **Complementary rows:** I-mod, I-confirm (labeled secondary).  
3. **Contrast rows:** Mode C companion, Mode S S1/S2, free OR (secondary), abs-z frozen reference.  
4. **L3 diagnostic subset:** records/windows with high continuous excess3 / high \(\Delta\)excess3 (readout only) — does I0 fire when OPC is silent?  
5. Success sketch: I0 FAR not saturating controls under S1-like gates **and** I0 captures a non-empty share of high-excess3 approaches that OPC misses.

---

## 5. RECD / Nivel-3 respect criteria

### 5.1 Pass

| ID | Criterion |
|----|-----------|
| **IN-P1** | Primary alarm is persistence on surplus (\(A^{\mathrm{core}}\)), not on collapse. |
| **IN-P2** | Collapse-only streams do not alarm under any shipped role. |
| **IN-P3** | High-support high-surplus streams can alarm on I0 / I-tag / I-mod without requiring \(\ell^C=1\). |
| **IN-P4** | OPS / \(S_t\) labeled Level-3–consistent **proxy**, not excess3. |
| **IN-P5** | No \(\mu/\sigma\) or abs-z as core of the integrated detector. |
| **IN-P6** | I-confirm’s L3 cost is disclosed; I-confirm is not sold as the default RECD-aligned arm. |
| **IN-P7** | Free OR and collapse-dominant fusions are rejected as the integrated primary story. |

### 5.2 Fail

| ID | Criterion |
|----|-----------|
| **IN-F1** | Collapse persistence is the main predicate; surplus is optional decoration. |
| **IN-F2** | Free OR presented as “the integrated detector” without surplus primacy. |
| **IN-F3** | I-confirm (or any collapse-required gate) claimed to fully respect all high-support L3 paths. |
| **IN-F4** | Continuous \(z\) on excess3 rebranded as OPSP. |
| **IN-F5** | Fabricated multi-cohort FAR/sens for I* without a real run. |

### 5.3 Self-check for any new cell

1. Can a collapse-only series alarm? → If yes, **IN-F1**.  
2. Can a high-\(D\), high-\(S\) series alarm on the **default** config? → Need I0-like; if only I-confirm, mis-labeled.  
3. Is \(A_t=1\) possible with \(\ell^S\) never active (even under modulate)? → Bug / fail.  
4. Any \(\mu/\sigma\)? → **IN-F4**.

---

## 6. Primitives map

| Piece | Symbol | Implementation |
|-------|--------|----------------|
| Surplus \(S_t\) | OPS | `ordinal_synergy_surplus` |
| Mode-S-only persist | — | `ops_detect` |
| **Integrated OPSP** | I0 / I-tag / I-mod / I-confirm | **`opsp_integrated_detect`** |
| Joint codes | \(\sigma\) | `joint_symbols_from_factors` |
| Diversity \(D_t\) | complementary | internal `_windowed_diversity` |
| Collapse-only OPC | contrast | `opc_detect` / BR / G |
| Exports | — | `code/ordinal_detectors/__init__.py` |
| Tests | — | `tests/test_opc_refinements.py` |

Exploratory only — not production ship; does not replace frozen abs-z.

---

## 7. Ranked recommendation

1. **Adopt I0 as the integrated core** — basal-relative \(\Delta S\) + \(\theta_R^S\) persistence; collapse ignored for alarm. This is the cleanest N2-logic × N3-signal design.  
2. **Use I-tag in parallel** for subtype reporting (surplus±collapse) without changing decisions.  
3. **Explore I-mod** only if co-locking events are common and I0 misses short surplus runs that co-occur with collapse—pre-register \(\delta_R,\delta_S\).  
4. **Treat I-confirm as secondary FAR filter**, never as the sole “RECD integrated” claim.  
5. **Keep Mode C companion (and Mode-S-alone) as separate contrast rows** for bake-offs; do not delete separate-arm discipline.  
6. **Reject** collapse-dominant integration, free-OR-as-primary, and absolute-OPS-only production claims.

---

## 8. Bottom line

| Deliverable item | Where |
|------------------|--------|
| Integrated proposal | §0–§1: persist surplus primary |
| Concrete integration rules | §2: I0 / I-tag / I-mod / I-confirm |
| Evaluation | §4 conceptual + synthetic unit probes |
| Ontology vs separate arms | §3 |
| RECD / L3 checks | §5 |
| Shipped code | `opsp_integrated_detect` |

**One sentence:** The integrated detector applies Nivel-2 **persistence** to a Nivel-3–consistent **synergy surplus** signal, and allows collapse only as tag, soft modulate, or secondary confirm—so high-synergy reorganizations remain first-class even without diversity collapse, while pure collapse never masquerades as the primary alarm.
