# OPC refinements under RECD: specificity without blinding Level 3

**Status:** Conceptual / design exploration (not a bake-off; not a production retune)  
**Date:** 2026-07-13  
**Companion code (exploratory primitives only):** `code/ordinal_detectors/opc_refinements.py`  
**Does not retune:** frozen abs-z, production OPC companion \((L,\theta_D,\theta_R)=(50,0.35,5)\)

---

## 0. Framing and an honest tension

### 0.1 Empirical starting point (existing public arms only)

| Arm | All-event sens | FAR / 24 h (NSRDB) | Role |
|-----|----------------|--------------------|------|
| OPC \(L=50\) | \(\approx 0.424\) (14/33) | \(\approx 3.733\) | Specificity-leaning collapse |
| SDD (TV) | \(\approx 0.970\) | \(\approx 46.267\) | Sensitivity-leaning divergence |
| abs-z \(\tau_s\) (frozen) | \(\approx 0.909\) | \(\approx 33.734\) | Production primary |

Sources: `docs/ORDINAL_SENSITIVITY_SPECIFICITY_TRADEOFF.md`, `results/ordinal_opc_param_explore_report.md` (`keep_baseline`).  
No new empirical bake-off numbers are invented here.

### 0.2 What “improve specificity without blinding Level 3” can mean

Current OPC already has **low exploratory FAR** relative to abs-z / SDD. Further pure tightening of collapse (\(\theta_D \downarrow\), \(L \uparrow\), \(\theta_R \uparrow\)) can shave FAR further but **deepens** Level-3 blindness and already failed the modest grid’s “higher sens, FAR in slack” test (`keep_baseline`).

Therefore this note **separates** two design problems:

1. **Specificity levers** that reduce *non-informative* alarms (spurious collapse, basal already-collapsed regimes, chatter from hard resets) without requiring continuous \(\mu,\sigma\).
2. **Multi-mode ordinal visibility** so that high-synergy (Level-3–consistent) reorganizations are not structurally excluded from the *family* of native detectors—even if they are not collapse events.

Tightening collapse-only knobs addresses (1) at the expense of (2). Viable directions treat OPC as a **Φ₂-collapse arm inside a RECD-aligned family**, not as a single predicate that must do everything.

### 0.3 Hard constraints (non-negotiable here)

- No abs-z / continuous basal \(\mu,\sigma\) as proposed OPC core.
- No claim that OPC is superior overall; no clinical / FDA / S5 readiness.
- Level 3 (excess3 / irreducible surplus) must remain **representable** or at least not discarded by design.
- excess3 in the pilot is a continuous *readout* of an ordinal hierarchy. Integrating Level 3 **must not** smuggle “abs-z on excess3” and rebrand it as OPC. Any continuous readout used only as diagnostic must be labeled as such; alarm proposals below prefer **discrete** Level-3 proxies on the symbol stream.

---

## 1. How current OPC conflicts with Level-3 / high-synergy reorganization

### 1.1 What RECD Level 3 is (manuscript operational sense)

Level 3 isolates joint ordinal organization **not explained by Levels 1–2 alone**. Operationally the pilot tracks continuous **excess3**—a surplus contribution of the highest nested layer beyond lower-order structure—and \(\Delta\mathrm{excess3}\) as the approach-vs-basal Level-3 contrast. Empirically, \(\Delta\tau_s\) and \(\Delta\mathrm{excess3}\) are often **sign-concordant** (8/10 discovery records) even when classical variance stories reverse. That is: pre-VF reorganization is frequently a **relational depth** story, not a univariate CSD story.

Critically, Level 3 is **not** defined as “few active symbols.” High synergistic surplus can coexist with:

- moderate or high support diversity \(D_t\) (many joint symbols visited);
- rich AF-like alphabets;
- intermittent rather than contiguous locking.

### 1.2 What OPC computes

\[
D_t=\frac{|\mathrm{supp}(W_t)|}{K},\qquad
A_{\mathrm{OPC}}(t)=\mathbf{1}\{D_t\le\theta_D\}\wedge\mathbf{1}\{R_t\ge\theta_R\}
\]

with \(R_t\) a **hard consecutive** run of low-diversity windows. Positive class = *sustained collapse of ordinal repertoire*.

### 1.3 Why this structurally misses Level-3 paths

| Pre-VF path (RECD-consistent) | OPC response | Why |
|-------------------------------|--------------|-----|
| excess3 rises while many joint symbols remain active | Silent | \(D_t\) stays above \(\theta_D\); no collapse |
| Directed \(\tau_s\) reorganization (sign of \(\Delta\tau_s\)) without support collapse | Silent | Support cardinality carries no direction |
| Intermittent near-collapse clusters | Often silent | Single high-\(D\) window resets \(R_t\) |
| High-entropy AF basal with subtle synergy shift | Silent / unreachable | Collapse threshold hard to meet without FAR blow-up if loosened |
| True sustained collapse (Φ₂ locking of few states) | Can alarm | Overlap with OPC’s positive class |

The manuscript’s claim that OPC is a “path toward \(\Phi_3\)” was **aspirational**: prolonged low diversity *may* concentrate mass on few joints (sometimes necessary for deeper conjunction), but it is **neither necessary nor sufficient** for Level-3 surplus change. Collapse is a Φ₂-adjacent locking signature, not a Level-3 detector.

### 1.4 Four known limits (context for proposals)

1. **Level-3 / excess3 blindness** — ontological for collapse-only OPC.  
2. **Intermittency / persistence reset** — operational fragility of hard \(R_t\).  
3. **High-entropy substrates** — collapse rare when many symbols remain active.  
4. **Loss of directional relational information** — support ignores \(\mathrm{sign}(\Delta\tau_s)\).

Any refinement that only tightens \((L,\theta_D,\theta_R)\) addresses specificity in a narrow band while **worsening** (1) and (3). Grid evidence already shows sensitivity gains require FAR inflation (`keep_baseline`).

---

## 2. Design principle: multi-mode ordinal family, mode-specific specificity

```
                    RECD symbol stream σ_t
                           │
           ┌───────────────┼───────────────┐
           ▼               ▼               ▼
     Mode C: Collapse   Mode S: Ordinal   Mode R: Relation
     (classic OPC /     synergy surplus   persistence
      OPC-G)            (Level-3 proxy)   (Φ₂ codes)
           │               │               │
     specificity       specificity       specificity
     gates (persist,   gates (basal-     gates (run length,
     basal-relative D) relative S,       refractory)
                       refractory)
           │               │               │
           └───────────────┴───────────────┘
                           │
              report arms SEPARATELY first;
              free OR/AND only as labeled secondary
```

**Specificity is improved arm-by-arm** (refractory, basal-relative gates, gap policy, longer mode-specific persistence), not by forcing every transition through collapse.

**Level-3 visibility** lives in Mode S (and partly Mode R), not by pretending Mode C sees synergy.

---

## 3. Concrete proposals

### Proposal A — OPC-G: gap-tolerant persistence (intermittency without continuous moments)

**Mechanism.**  
Replace hard reset of \(R_t\) with a **gap budget** \(G \ge 0\):

- If \(D_t \le \theta_D\): increment persistence credit; reset gap counter.  
- If \(D_t > \theta_D\) and gap_counter \(< G\): allow a gap (do not zero credit; increment gap_counter).  
- If more than \(G\) consecutive high-diversity windows: reset credit to 0.

Alarm when persistence credit \(\ge \theta_R\) and current window is in an allowed active state (low-div, or within an open gap run after enough low-div mass—implementation default: require low-div at alarm time *or* credit \(\ge \theta_R\) with last low-div within \(G\) steps; see code).

**Intended FAR / specificity effect.**  
Alone, \(G>0\) tends to **raise** sensitivity and **risk FAR** (recovers interrupted collapses). To *improve* specificity while using gaps:

- pair with **stricter** \(\theta_D\) or larger \(\theta_R\), **or**
- require a minimum number of low-div windows inside the credit window (e.g. \(\ge \theta_R\) low-div hits, not merely credit length), **or**
- add a control **refractory** after each episode (already used in FAR protocol; can be part of detector policy).

The specificity win is **anti-chatter**: intermittent noise that never accumulates enough low-div mass still fails; contiguous false collapses remain gated by \(\theta_R\).

**Level-3 / excess3 visibility.**  
**Does not restore Level 3.** Still collapse-based. Honest status: fixes limit (2) only.

**RECD philosophy.**  
Compatible: discrete, no \(\mu/\sigma\), Φ₂ persistence language preserved (softened memory of locking).

**Viability:** **Promising as a Mode-C refinement** (conditional on joint retune of \(G\) with \(\theta_D/\theta_R\), pre-registered). Not a Level-3 solution.

**Reject if:** \(G\) is used alone as a sensitivity rescue and then marketed as “more specific.”

---

### Proposal B — OPS: Ordinal Persistence of Synergy surplus (discrete Level-3 proxy)

**Mechanism.**  
On bivariate ordinal factors \((\pi^{(1)}_t,\pi^{(2)}_t)\) (or joint code with known factorization into margins), for window \(W_t\):

1. Estimate joint empirical law \(P_t(\pi_1,\pi_2)\).  
2. Estimate product of margins \(P_t(\pi_1)P_t(\pi_2)\).  
3. Define a **synergy / dependence surplus** (purely discrete), e.g.

\[
S_t = \mathrm{TV}\big(P_t,\; P_t^{(1)}\otimes P_t^{(2)}\big)
= \tfrac12\sum_{\pi_1,\pi_2}\big|P_t(\pi_1,\pi_2)-P_t(\pi_1)P_t(\pi_2)\big|
\]

or discrete mutual information \(I_t(\pi_1;\pi_2)\) (same qualitative role).  
4. Low-synergy / high-synergy indicators relative to a **fixed early basal ordinal law** (not continuous \(z\)):

\[
\ell^S_t=\mathbf{1}\{S_t \ge \theta_S\}\quad\text{or}\quad
\mathbf{1}\{S_t - S_{\mathrm{basal}} \ge \theta_{\Delta S}\}.
\]

5. Persist \(\ell^S\) with hard run or OPC-G-style gaps → alarm \(A_{\mathrm{OPS}}\).

**Why this is Level-3–aligned.**  
Level 3 is surplus joint organization beyond lower-order structure. \(S_t\) is exactly a windowed, discrete measure of joint structure beyond independent margins—**the same conceptual axis as synergistic surplus**, without requiring abs-z on continuous excess3. High \(S_t\) can occur with **high** support diversity (many coordinated joints), so OPS is **not** collapse-blind in the Level-3 sense: it can fire when OPC is silent.

**Intended FAR / specificity effect.**  
Naive OPS (absolute \(\theta_S\), short persistence) will likely **raise FAR** on ambulatory Holters (dependence fluctuates). Specificity levers that stay ordinal:

- basal-relative \(\Delta S\) (only *change* in surplus, not raw dependence);  
- longer \(\theta_R^S\);  
- refractory episode counting;  
- optional dual gate: OPS alarms only if \(\Delta S\) high **and** support is *not* pure noise (e.g. TV vs basal symbol law above a floor)—still discrete.

**RECD philosophy.**  
Strong alignment with Level 3 **if labeled as a separate arm** (Mode S), not as a silent redefinition of OPC. Continuous excess3 remains the pilot’s hierarchical readout; OPS is a **native alarm-layer proxy** on the same symbol factors.

**Risks / honesty.**

- OPS is **not** identical to excess3; it is a dependence surplus on factorized bivariate symbols, while excess3 encodes nested RECD hierarchy with its own \(\theta_3\) machinery. Treat OPS as **Level-3–consistent**, not as a drop-in replacement for excess3 science.  
- Do **not** reintroduce abs-z on excess3 and call it OPS.  
- Free \(A_{\mathrm{OPC}}\vee A_{\mathrm{OPS}}\) without separate FAR accounting is a fusion experiment, not a free upgrade.

**Viability:** **Most promising for Level-3 visibility** inside ordinal space. Specificity is **conditional** on basal-relative gates + persistence; not free.

---

### Proposal C — OPC-BR: basal-relative collapse gate (specificity on Mode C)

**Mechanism.**  
Classic OPC uses absolute \(\theta_D\). Instead (or additionally):

\[
\ell^{\mathrm{rel}}_t=\mathbf{1}\Big\{D_t \le \min(\theta_D,\; \rho\cdot D_{\mathrm{basal}})\Big\}
\quad\text{or}\quad
\mathbf{1}\{D_{\mathrm{basal}}-D_t \ge \delta_D\},
\]

with \(D_{\mathrm{basal}}\) the support diversity on a fixed early basal segment (ordinal only). Persist with hard or gap-tolerant runs.

**Intended FAR / specificity effect.**  
**Directly targets false collapses** when the basal repertoire is already small (paced segments, night-time low variability, artifact-sparse alphabets): absolute \(\theta_D\) may fire “collapse” that is just basal normality. Requiring *relative* collapse should **lower FAR** on such controls. On event Holters that truly collapse from a rich basal, sensitivity is preserved; if basal is already collapsed, Mode C stays appropriately quiet (Mode S/R may still fire).

**Level-3 visibility.**  
**No.** Still collapse mode. Complements OPS rather than replacing it.

**RECD philosophy.**  
Compatible: basal is a reference ordinal law (as SDD already uses), not continuous moments of \(\tau_s\).

**Viability:** **Promising, low-complexity Mode-C specificity lever.** Should be explored before further absolute \(\theta_D\) tightening.

---

### Proposal D — Dual-mode ordinal locking (collapse **or** high-synergy lock), reported as a family

**Mechanism.**  
Define two predicates on the same stream:

- \(A_C\): OPC or OPC-G + OPC-BR (collapse lock).  
- \(A_S\): OPS (synergy surplus lock).  
- Optional \(A_R\): persistence of a **discrete relation code** (e.g. sustained EQ/GT/LT or fixed joint transition)—Φ₂ proper, not support cardinality—if relation codes are already available from RECD machinery.

**Primary reporting rule (recommended):** score \(A_C\), \(A_S\) (and \(A_R\)) **separately** with their own sens/FAR.  
**Secondary only:** labeled rules such as

- \(A_C \vee A_S\) (union: maximize ordinal coverage; expect FAR ↑),  
- \(A_S\) with confirm \(A_C\) within \(\pm w\) (or vice versa)—cascade costs already documented for SDD→OPC.

**Intended FAR / specificity effect.**  
- Mode C alone: push FAR **down** via OPC-BR + refractory + (optional) stricter absolute thresholds.  
- Mode S alone: control FAR via basal-relative \(\Delta S\) + long persistence—not by forcing collapse.  
- Union: **not** a specificity improvement; useful only as a sensitivity upper bound for “any ordinal lock.”

**Level-3 visibility.**  
**Yes**, via Mode S (and partially Mode R). Collapse no longer monopolizes the definition of “native ordinal alarm.”

**RECD philosophy.**  
Best match to nested RECD: different levels / phenomena get different alarm geometries. Avoids the category error “one diversity threshold = whole RECD.”

**Viability:** **Recommended architectural direction.** Empirical work should be multi-arm, not a single re-tuned OPC cell.

---

### Proposal E — (Rejected / not viable as stated) Tighter collapse-only grid as Level-3 fix

Further absolute tightening of \((L,\theta_D,\theta_R)\) without multi-mode structure:

- May lower FAR slightly;  
- Cannot create Level-3 sensitivity;  
- Already explored modestly → `keep_baseline`;  
- Contradicts the goal of remaining able to see high-synergy paths.

**Viability:** **Not viable** as the answer to this goal. Useful only as a Mode-C FAR polish *after* Mode S exists.

---

### Proposal F — (Rejected) abs-z on continuous excess3 labeled as “ordinal OPC”

Would reintroduce continuous basal \(\mu/\sigma\) on a scalar series. May detect Level-3 amplitude excursions but **violates** the ordinal-native constraint and the project’s frozen-abs-z separation.  
**Viability:** **Not viable** under stated constraints. Continuous excess3 remains a scientific readout; alarm design for it—if ever revisited—must be a **named continuous arm**, not OPC.

---

## 4. Per-proposal summary matrix

| Proposal | Mechanism (short) | FAR / specificity intent | Level-3 visibility | RECD stance | Viability |
|----------|-------------------|---------------------------|--------------------|-------------|-----------|
| **A OPC-G** | Gap-tolerant \(R_t\) | Anti-chatter; pair with stricter collapse for net FAR↓ | No (still collapse) | Φ₂ persistence OK | Promising (Mode C) |
| **B OPS** | TV/MI joint vs product margins + persist | FAR control via basal-relative \(\Delta S\), long run | **Yes** (discrete surplus) | Level-3–consistent proxy | **Most promising for L3** |
| **C OPC-BR** | Collapse relative to basal \(D\) | FAR↓ when basal already low-div | No | Ordinal basal OK | Promising (Mode C) |
| **D Dual-mode family** | Separate \(A_C\), \(A_S\) (+ optional \(A_R\)) | Specificity per arm; no free OR | **Yes** (via \(A_S\)) | Nested RECD best match | **Recommended architecture** |
| E Tighter collapse-only | \(\theta_D\downarrow\) etc. | FAR↓ possible | **Worse** | Misaligned for L3 goal | **Not viable** as solution |
| F abs-z on excess3 | Continuous \(z\) | Spec. of another continuous arm | Sees L3 amplitude | **Breaks** ordinal-native OPC | **Rejected** |

---

## 5. Conceptual trade-off (no new bake-off)

| Design move | Expected sens (conceptual) | Expected FAR (conceptual) | Notes |
|-------------|----------------------------|---------------------------|--------|
| Stricter absolute OPC only | ↓ | ↓ | Deepens L3 blindness; grid already discourages free sens gains |
| OPC-G alone | ↑ | ↑ or flat | Needs co-tuning |
| OPC-BR alone | ~ or slight ↓ | ↓ on low-basal-div controls | Cheap specificity |
| OPS alone (absolute \(S\)) | ↑ on synergy events | ↑ risk | Needs basal-relative gate |
| OPS basal-relative + long persist | ↑ on true \(\Delta S\) events | controlled | Main L3 arm candidate |
| Dual-mode, report separate | per-arm | per-arm | Correct scientific reporting |
| Dual-mode free OR | ↑↑ | ↑↑ | Secondary only; multiplicity cost |
| SDD→OPC cascade (existing) | ~OPC | ~OPC | Does not open L3 |

Honest summary: **you cannot make collapse-only OPC both more specific and Level-3–aware.** You can (i) make Mode C more specific, and (ii) add Mode S so the *family* is Level-3–aware, with mode-specific FAR accounting.

---

## 6. Ranked recommendation (what to explore next)

### Priority 1 — Dual-mode architecture (Proposal D) with OPS (B) + OPC-BR (C)

1. Implement / freeze formal OPS on bivariate factor symbols (TV vs product margins + basal-relative option)—primitives started in `opc_refinements.py`.  
2. Add OPC-BR to Mode C (basal-relative collapse).  
3. Score **three primary exploratory arms separately** on the existing public protocol:  
   - Mode C: OPC companion ± BR ± G  
   - Mode S: OPS  
   - Frozen abs-z (reference only; not retuned)  
4. Do **not** lead with free OR. Optional secondary: cascade Mode S → Mode C confirm (symmetric to existing SDD→OPC lesson).

### Priority 2 — OPC-G (A) as Mode-C intermittency polish

Only after BR, and only with pre-registered \((G,\theta_R,\theta_D)\) co-constraints so “gaps” do not become a silent sensitivity cheat.

### Explicitly deferred

- Full multi-cohort bake-off and production ship.  
- Retuning frozen abs-z.  
- Clinical / S5 / FDA claims.  
- Claiming OPS ≡ excess3.  
- Manuscript promotion of dual-mode as superior without data.  
- Free fusion with SDD as primary.

### What not to do next

- Another absolute collapse grid expecting Level-3 rescue.  
- Rebranding abs-z(excess3) as “native ordinal.”  
- Union rules without separate FAR tables.

---

## 7. Relation to Contexto.docx and prior formalization

`Contexto.docx` introduced OPC as collapse+persistence and SDD as distributional divergence, keeping options **separate first**. That separation remains correct. This note extends the program by:

- admitting that OPC’s “path to \(\Phi_3\)” language over-promised;  
- adding a **third ordinal arm class** (OPS / Mode S) for Level-3–consistent surplus;  
- using basal-relative and gap machinery as **specificity tools**, not as continuous \(z\)-scores.

---

## 8. Minimal formal sketch (OPS + OPC-G) for later bake-off registration

**OPC-G.** Parameters: \(L,\theta_D,\theta_R,G\). Credit updates as in §3A.  

**OPS.** Parameters: \(L_S\), \(\theta_{\Delta S}\) or \(\theta_S\), \(\theta_R^S\), optional \(G_S\); basal length fixed early.  

**OPC-BR.** Parameters: \(\rho\) or \(\delta_D\) plus standard OPC.  

**Registration rule (suggested):** freeze Mode C companion first; add BR; only then open a small OPS grid on held-out logic—never retune abs-z.

---

## 9. Bottom line

- Current OPC is a **good specificity-leaning Φ₂-collapse detector** and a **bad Level-3 detector**. That is structural.  
- Further FAR reduction should come from **basal-relative collapse, refractory discipline, and careful gap policy**—not from deeper absolute collapse alone.  
- Level-3 respect requires a **separate ordinal synergy/surplus arm (OPS)** inside a dual-mode family, scored honestly, without abs-z on continuous excess3.  
- Most promising exploration path: **D + B + C**, then **A**; reject **E** and **F** as answers to this goal.

---

## 10. Follow-on specificity note

For Mode-C FAR levers (OPC-BR / OPC-G co-tuning), Mode-S OPS gates, separate-arm evaluation design, concrete parameter sheets (C-BR1…, S1…), and Level-3 pass/fail criteria, see:

**`docs/OPC_MODE_C_S_SPECIFICITY_REFINEMENT.md`**
