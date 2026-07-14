# Delivery note — §3.11 Toward native ordinal detectors (rigor upgrade)

**Date:** 2026-07-14 (updated: I0 surplus + structural + final ranking)  
**Status:** Integrated into `manuscript/CCTP_SDBB_manuscript.md` (replaced prior §3.11 in place; extended with surplus track + ranking).  
**Companion full-text extract:** same prose lives in-manuscript; no parallel competing section.

---

## (a) Full proposed section text

See **`manuscript/CCTP_SDBB_manuscript.md`**, heading `## 3.11 Toward native ordinal detectors (exploratory)` (immediately after §3.10, before `# 4. Discussion`).

Structural blocks present:

| Block | Location in §3.11 |
|-------|-------------------|
| Original OPC motivation | *Motivation for native ordinal alternatives* + *OPC* definition |
| Trade-off / specificity strengths / non-superiority | *Exploratory sensitivity–FAR trade-off* + Table `tab:ordinal_tradeoff` |
| Four structural limitations | *Structural limitations of OPC (four named weaknesses)* — paragraphs (1)–(4) |
| Process self-critique | *Process self-critique: reactive discovery of structural limits* |
| Keep baseline + cascade secondary | *Modest OPC parameter exploration*; *Light cascade…* |
| **I0 surplus + structural arms** | *Surplus-primary I0 and structural arms* (0/16 grid; 0/6 structural; stop) |
| **Final ranking for pre-VF hit rate** | *Practical ranking for pre-VF / ventricular-arrhythmia event detection* |
| Frozen primary abs-z | Stated in mismatch, ranking, and *Honest reading*; echoed in §4.1 + §6 |

---

## (b) Placement recommendation

**Keep as §3.11**, after Phase 2 public FAR (§3.10) and before Discussion (§4).

Rationale:

1. The frozen abs-z specificity bottleneck must be fully stated before ordinal alternatives appear; otherwise OPC’s lower FAR can be misread as the main result of the paper.
2. Discussion §4.1 already frames “relational signal real / abs-z operating point not yet specific.” The upgraded §3.11 feeds that discussion with a sharper second clause: *even an ordinal FAR-leaning alternative is not a free upgrade*, because OPC discards Level-3 and directional τ_s information.
3. Do **not** promote §3.11 into Methods as a primary detector section, and do **not** duplicate it under Discussion. One epilogue is enough.
4. Optional one-line pointer in Discussion §4.1 (if a later edit pass touches Discussion): mention that ordinal alternatives expose a trade-off surface and that OPC’s structural limits were under-stressed at design time. Not required for this delivery.

---

## (c) Narrative-balance guidance (achievements vs limitations)

| Weight | Content |
|--------|---------|
| **Keep strong** | Discovery relational message (§3.1–3.8); Phase 1/2 external honesty on abs-z FAR; existence of a clean three-arm trade-off surface. |
| **Keep explicit but not triumphant** | OPC FAR ≈ 3.73/24 h on NSRDB under fixed params; sens/FAR efficiency if someone prioritizes specificity. |
| **Elevate (this upgrade)** | Four structural OPC limits; keep_baseline as evidence that the miss is not a single bad cell; process self-critique (reactive vs anticipatory). |
| **Never imply** | OPC “better than abs-z overall”; I0 “better than abs-z”; clinical/S5/FDA readiness; cascade as rescue; that Level-3 alignment of OPC was operational rather than aspirational. |
| **Final ranking (locked)** | For pre-VF event hit rate on public Holter: **abs-z τ_s first**; SDD sensitivity ceiling; I0 not promoted; OPC specificity companion only. |

**Balance rule of thumb for future prose:** for every sentence that praises OPC’s specificity, require a nearby sentence that states its sensitivity cost *and* one named structural reason (synergy blindness, high-entropy substrate, intermittency reset, or loss of directional τ_s). Do not let FAR reduction stand alone as the takeaway of §3.11. For every mention of I0/surplus, restate non-promotion vs abs-z.

**What was achieved (honest):** rigorous discrete predicates; independent evaluation; transparent numbers; keep_baseline discipline; correction of an over-elegant design story.

**What was sacrificed (honest):** event hit rate; Level-3 visibility; directed relational information; design-time adversarial stress-testing.

---

## Numbers locked to artifacts

| Claim | Source |
|-------|--------|
| OPC sens all 0.424 (14/33), FAR 3.733 | `docs/ORDINAL_SENSITIVITY_SPECIFICITY_TRADEOFF.md` |
| SDD 0.970 / 46.267 | same |
| abs-z 0.909 / 33.734 | same |
| keep_baseline L=50, θ_D=0.35, θ_R=5 | `results/ordinal_opc_param_explore_report.md` |
| Cascade sens 0.394, FAR ≈3.87 | prior §3.11 / cascade doc |
| I0 default ~0.79 / FAR ~40; grid 0/16 | `docs/I0_SURPLUS_PARAM_GRID.md` |
| Structural R0–R5 0/6; stop | `docs/I0_STRUCTURAL_ARMS.md` |
| abs-z preferred primary for event hit rate | §3.11 ranking + §4.1 + §6 Conclusions |

---

## Delta vs prior §3.11 (summary)

- Prior: trade-off surface + keep baseline + soft “honest reading.”
- Mid: dedicated **four structural limitations**; **process self-critique**; OPC **not** a net advance over abs-z.
- **2026-07-14:** I0 surplus grid + structural arms closed; **final ranking** for ventricular pre-event detection written into §3.11, Discussion §4.1, and Conclusions §6 (abs-z primary; OPC companion; I0 not promoted).
