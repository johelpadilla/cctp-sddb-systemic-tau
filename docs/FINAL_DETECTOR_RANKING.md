# Final detector ranking (public Holter bake-off)

**Status:** Locked into the final report narrative  
**Date:** 2026-07-14  
**Authoritative prose:** `manuscript/CCTP_SDBB_manuscript.md` §3.11 (ranking subsection), §4.1, §6 Conclusions  
**No clinical / FDA / deployability claim.**

---

## Question answered

Among the detectors measured on SDDB+VFDB (events) and NSRDB (controls), which is least wrong when the Holter **actually** harbors pre-ventricular-event / SCD-class dynamics?

---

## Ranking (event hit rate priority)

| Rank | Detector | Role in final report | Sens (all events) | FAR/24h (NSRDB) |
|------|----------|----------------------|-------------------|-----------------|
| **1** | **abs-z on continuous \(\tau_s\)** | **Preferred primary** operating point | ≈ **0.91** | ≈ 33.7 |
| **2** | SDD | Sensitivity ceiling only (not practical primary) | ≈ 0.97 | ≈ 46.3 |
| **3** | I0 / OPS surplus-primary | Mode-S narrative / L3 proxy; **not promoted** | ≈ 0.76–0.79 | ≈ 40 |
| **4** | OPC \(L=50\) | **Specificity-leaning ordinal companion** | ≈ 0.42 | ≈ **3.73** |
| **5** | I-confirm | Secondary FAR filter; L3 cost | ≈ 0.27 | ≈ 1.73 |

Structural surplus best dual lever (R3 MAD): sens ≈ 0.61, FAR ≈ 23 — still not competitive with abs-z for hit rate nor OPC for FAR.

---

## One-sentence takeaway

**Parametric abs-z on Systemic Tau \(\tau_s\) remains the preferred primary for not missing public pre-ventricular events; OPC is the clean low-FAR ordinal reference; surplus I0 is closed as a FAR-competitive primary.**

---

## Supporting closed experiments

| Experiment | Outcome |
|------------|---------|
| OPC 36-cell grid | keep_baseline |
| I0 \(4\times4\) \(\theta\)-grid | **0/16** clear advances vs OPC |
| Structural R0–R5 | **0/6** structural wins; **stop surplus-primary** |

Details: `docs/I0_SURPLUS_PARAM_GRID.md`, `docs/I0_STRUCTURAL_ARMS.md`, `docs/ORDINAL_SENSITIVITY_SPECIFICITY_TRADEOFF.md`.
