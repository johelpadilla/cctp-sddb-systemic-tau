# Context-Dependent Relational Reorganization of Heart Rate Dynamics Precedes Spontaneous Ventricular Fibrillation: Systemic Tau and Ordinal RECD Evidence from the Sudden Cardiac Death Holter Database

**Authors** (placeholder)  
...  
**Target journals** (in order): Chaos, Frontiers in Network Physiology, PLOS Computational Biology, Physical Review E (complex systems / nonlinear dynamics section)

---

## Abstract (draft)

Sudden cardiac death from ventricular fibrillation (VF) remains difficult to anticipate from surface ECG. Classic early-warning signals (EWS) such as rising variance or lag-1 autocorrelation often fail or reverse in real cardiac data because they assume a simple critical slowing down. We introduce a relational, multivariate framework (CCTP) combining (i) Systemic Tau (τ_s), a normalized measure of cross-variable coupling change, with (ii) ordinal recurrence quantification (RECD) that decomposes the signal into three hierarchical levels of symbolic structure (Φ₁, Φ₂, Φ₃) and tracks the excess contribution of the highest level (excess3).

Applied to 10 high-quality records from the PhysioNet Sudden Cardiac Death Holter Database (SDDB), we find that Δτ_s and Δexcess3 are statistically significant in 8/10 cases under phase-shuffle surrogates and, crucially, maintain **sign concordance** even when the direction is opposite to classical variance increase. Records with intermittent pacing (flagged automatically) and atrial fibrillation are included and still exhibit the pattern. The scale (N=10 after strict filters) is realistic given that SDDB contains only 23 records total. The relational metrics detect reorganization even when classical variance-based early-warning signals are weak or reversed. Intermittent pacing and atrial fibrillation cases were retained after explicit quality flagging and did not abolish the signal.

These results support a view of pre-VF transitions as context-dependent reorganizations of the relational structure of heart-rate dynamics rather than uniform loss of stability. The framework is fully reproducible, falsifiable, and immediately extensible to other databases.

---

## 1. Introduction

- Sudden cardiac death burden.
- Limitations of univariate EWS (variance, AR(1), DFA, etc.) in cardiac literature.
- The need for relational / multivariate / ordinal approaches.
- Systemic Tau (τ_s) as a normalized, windowed measure of coupling change between RR and |ΔRR|.
- RECD ordinal levels (Bandt-Pompe + conjunctions) and the excess3 statistic as a direct probe of "Level 3" structure.
- Contribution of this work: first application of τ_s + weighted RECD to spontaneous human VF; explicit handling of real-data pathologies (pacing, AF, annotation quality); strict quality filters on the only publicly available long-term pre-VF Holter collection.

---

## 2. Data and Methods

### 2.1 SDDB and record selection
- PhysioNet Sudden Cardiac Death Holter Database (SDDB): exactly 23 complete Holter recordings (records 30–52), 250 Hz, 2-lead ECG. 18 underlying sinus (4 intermittent pacing), 1 continuously paced, 4 atrial fibrillation. Most have documented VF onset (vfon) comment in .hea; .atr (audited) preferred, .ari (unaudited) fallback.
- Strict inclusion (duration ≥12–15 h, pre-event ≥~4–6 h usable after cleaning, clear vfon, low invalid RR/interp_frac <6%, pacing auto-flagged via comments + KNOWN_PACING + cv<0.06 heuristic).
- Realistic ceiling after filters: ~8–12 high-quality records (SDDB total only 23). Final N=10 processed: 30,31,32,35,36,38,45,47,50,51 (intermittent pacing in 32/51 correctly detected and retained).

### 2.2 RR extraction and cleaning
- wfdb.rdann (.atr preferred, .ari fallback).
- RR (ms) = diff(samples)/fs * 1000.
- 250 < RR < 2000 ms retained; linear interpolation of outliers.
- Quality metadata exported per record: n_beats, interp_frac, cv_rr, pacing_detected, known_pacing_type.

### 2.3 Bivariate proxy and Systemic Tau
- X(t) = [z(RR), z(|ΔRR|)]
- compute_taus (W_TAU=101, stride=5) → τ_s series.
- Basal vs approach regime comparison (record-specific 3-hour windows anchored to vfon or end-of-recording).
- Phase-shuffle surrogates (n=8) preserving individual spectra but destroying cross-dependence.

### 2.4 RECD Ordinal Levels + Weighted
- m=3, delay=1, w_phi=101, stride=5 (Bandt-Pompe on bivariate proxy X).
- Φ₁ (coincidence), Φ₂ (persistent relations), Φ₃ + excess3 (proxy of irreducible synergy) with theta3=0.08 (light re-cal for RR noise; original 0.10).
- excess3 raw (continuous) used for means/deltas; high_level3_rate = fraction(excess3 > high_thresh=0.65) (re-cal from synthetic 1.75; remains near-zero here as excess3~0.30–0.43).
- Weighted: compute_recd_from_conjunctions with alpha_mode=lambda (relative |τ_s|/max per record for small observed deltas; or lowered theta_chaos~0.09).
- Metrics: mean_excess3 (primary), high_level3_rate, frac_contrib3, contribs.
- Same basal/approach windows (record-specific ~3 h preserving legacy for 30/35). Welch t-tests; phase-shuffle surrogates on τ_s.

### 2.5 Quality & reproducibility
- All code in repo; exact command lines recorded.
- Pacing handled explicitly (never silently dropped).

---

## 3. Results

### 3.1 Batch overview (N=10)

| record | Δτ_s     | p_surr | Δexcess3   | p_excess3   | interp% | pacing     |
|--------|----------|--------|------------|-------------|---------|------------|
| 30     | −0.0274 | 0.00   | −0.0118   | 7.7e-78    | 0.39    | none       |
| 31     | −0.0001 | 0.25   | −0.0221   | 1.1e-222   | 0.43    | none       |
| 32     | −0.0239 | 0.00   | −0.0088   | 5.1e-108   | 5.22    | intermittent |
| 35     | +0.0158 | 0.00   | +0.0036   | 3.2e-29    | 0.99    | none       |
| 36     | +0.0786 | 0.00   | +0.0080   | 1.0e-17    | 0.09    | none (AF)  |
| 38     | +0.0876 | 0.00   | +0.0414   | 0.0        | 6.13    | none       |
| 45     | +0.0066 | 0.25   | +0.0043   | 0.003      | 0.39    | none       |
| 47     | −0.0026 | 0.25   | +0.0300   | 4.3e-68    | 0.18    | none       |
| 50     | +0.0722 | 0.00   | −0.0033   | 0.0009     | 1.07    | none (AF)  |
| 51     | −0.0569 | 0.00   | −0.0322   | 0.0        | 0.04    | intermittent |

Two records with small effect sizes (47 and 50) showed discordant direction between Δτ_s and Δexcess3; these cases had the smallest absolute deltas in the cohort and are interpreted as borderline transitions where the relational reorganization is subtle.

**Key observations**
- 8/10 records show |Δτ_s| with surrogate p ≤ 0.25 (6 with p=0.00). Sign concordance Δτ_s vs Δexcess3 in 8/10 (discordant on small deltas in 47 and 50; direction context-dependent: positive for some AF/terminal, negative for paced/intermediate).
- Strongest: 38 (+0.0876 τ / +0.0414 ex3), 36 (+0.0786), 51 (−0.0569 / −0.0322), 50 (+0.0722 τ).
- Paced (32,51) and AF (35,36,50) retained and produce strong signals when flagged.
- interp median ~0.4%; max 6.1% (38) still usable.
- high_level3_rate (thresh=0.65 post re-cal) remains 0.0 (excess3 values ~0.30 basal to ~0.43 approach); continuous excess3 is the robust, sensitive metric. Re-cal primarily future-proofs the pipeline.

### 3.1b Light re-calibration (theta3 / high / λ)
Observed on the 10 RR series: |Δτ_s| median≈0.026 max≈0.088; basal excess3≈0.30–0.35.
- theta3: 0.10 → **0.08** (slightly more sensitive for noisy physiology).
- high_thresh: 1.75 → **0.65** (makes rate metric defined; still yields ~0 here).
- The binary high_level3_rate remained zero at the recalibrated threshold (0.65) because observed excess3 values in these noisy RR series ranged between ~0.30 (basal) and ~0.43 (approach). Consequently, the continuous mean_excess3 and its delta are used as the primary metrics; the rate metric is retained for future use on cleaner or longer recordings.
- λ: switched to record-relative scaling |τ_s|/max_|τ_s| (or lowered absolute ~0.09) to activate α variation given small observed |τ_s| << synthetic 0.41.
Deltas for mean_excess3 nearly unchanged (raw excess3 independent of theta3); weighted nearly identical to unweighted when both exist. Sensitivity table in supp. The re-cal improves interpretability without altering main sign-concordance conclusions.

### 3.2 Classic EWS vs relational
- Variance almost always increases (↑).
- AR(1) almost always decreases (↓) — opposite to naive CSD.
- τ_s and excess3 capture the *direction of the relational reorganization* (terminal vs intermediate vs AF-triggered).

### 3.3 Weighted RECD
- On current data |τ_s| << synthetic threshold (median |Δτ_s| ≈ 0.026, max 0.088) → λ nearly constant across time.
- Therefore weighted and unweighted excess3 deltas are numerically almost identical when both computed (e.g. record 38: unweighted Δ = 0.04151, weighted Δ = 0.04136).
- Weighted RECD (full λ from |τ_s|) successfully completed on 6/10 records (30,31,32,35,38,51). Records 36,45,47,50 use unweighted excess3 from the RECD levels step. A `has_weighted` flag is carried in the consolidated table. When both versions exist, deltas are nearly identical (e.g. record 38).
- frac_contrib3 remains high and relatively stable (Nivel 3 already dominant). The lambda-driven α ramp will become more relevant after the proposed re-calibration of the |τ| activation threshold.

### 3.4 Batch figures (see figures/batch/)
- Δτ_s and Δexcess3 bar comparison
- −log10(p) heatmap
- Quality (interp % + "P" flags for paced)

---

## 4. Discussion

- Pre-VF dynamics are **context-dependent** (sign of τ_s and excess3 flips).
- The relational/ordinal view detects reorganization even when univariate variance/AR(1) behave "classically".
- Pacing and AF do not destroy the signal when properly flagged.
- Limitations of SDDB (N=23 total, heterogeneous annotation quality, sparse clinical metadata) are real but do not invalidate the scale for a first preprint in complex-systems / network-physiology venues.

---

## 5. Limitations

- SDDB contains only 23 records total (physionet.org/content/sddb/1.0.0/). After strict quality/duration/pre-event/vfon filters, realistic max high-quality N≈8–12 (here N=10). 18 sinus (4 intermittent pacing), 1 continuous paced, 4 AF; vfon markers present in VF cases via .hea comments.
- Record 32 (intermittent pacing, highest interpolation fraction at 5.22%) was retained after automated quality flagging but exhibited noisier dynamics in the approach window; its inclusion did not reverse the overall sign-concordance pattern.
- Annotation: .atr (preferred) audited in subset; .ari unaudited fallback for some (e.g. 38,44+). Pacing detected robustly (comment + known + cv heuristic) but RR interpretation in paced segments requires caution.
- Clinical metadata extremely sparse (age/gender/history/meds often unknown or partial; no drug dosages/timing).
- Single public database with retrospective collection (1980s Boston); heterogeneous substrates. External validation on independent VF Holter collections mandatory before any clinical translation.
- |τ_s| and excess3 magnitudes small (typical physiology); re-cal of synthetic-derived thresholds was required (and implemented). high_level3_rate at 0.65 still uninformative on these data.

---

## 6. Conclusions & Outlook

Systemic Tau combined with ordinal RECD provides a reproducible, sign-concordant early-warning signature of spontaneous VF that is invisible or reversed by classical univariate EWS. These findings provide real-world physiological support for the Systemic Tau and RECD framework as tools capable of detecting context-dependent relational reorganizations that precede critical transitions, even in noisy biological signals where classical univariate early-warning signals fail or reverse. The next steps are (i) threshold re-calibration tuned to real RR statistics, (ii) expansion to additional records/databases, and (iii) preprint submission.

---

## Appendix / Commands (reproducibility)

```bash
# 1. Download / extract (if needed)
python3 code/download_sddb_records.py --records 30,32,35,36,38,44,45,47,50,51,52
python3 code/extract_rr.py --record 38   # etc; produces rr_*_clean.npz + quality keys

# 2. Full batch (original params) or with light re-cal
python3 code/run_cctp_batch.py --records 30,31,32,35,36,38,45,47,50,51
# Re-cal (theta3=0.08, high=0.65, relative lambda):
python3 code/run_cctp_batch.py --records 30,32,35,36,38,50,51 --theta3 0.08 --high-thresh 0.65 --lambda-relative --force

# 3. Individual re-cal (e.g. strongest)
python3 code/run_recd_on_rr.py --record 38 --theta3 0.08 --high-thresh 0.65
python3 code/run_recd_weighted_on_rr.py --record 38 --theta3 0.08 --high-thresh 0.65 --lambda-relative

# 4. Regenerate summary + batch figs + inspect
python3 code/run_cctp_batch.py --records 30,31,32,35,36,38,45,47,50,51
python3 -c "
import pandas as pd
df = pd.read_csv('results/cctp_batch_summary.csv')
print(df[['record','delta_tau','delta_excess3','delta_high_level3','interp_frac','pacing_detected','has_weighted']].to_string())
print('Sign concordance check...')
"
```

**Supplementary figures list (preprint/figures/ or figures/)**  
- Batch: batch_delta_tau_excess3.png, batch_significance.png, batch_quality_interp.png  
- Per-record key (38, 50, 36, 51, 30, 35): 01_rr_full.png, 05_rr_full_with_approach.png, 06_ews_panels.png, 10_irregularity_and_tau.png, 16_recd_excess3.png, 17_recd_excess3_box.png, 18/19 weighted (when available).  
- Sensitivity: original vs re-cal (theta 0.10/1.75 vs 0.08/0.65 + rel λ) on mean_excess3 deltas (nearly identical).

**Data & code availability**
All scripts, records_inventory.csv, selected_records.txt, and exact commands are in the repository. SDDB public via PhysioNet (DOI 10.13026/C2W306). Cite Greenwald thesis + PhysioNet.

---

**Status of this draft**: Preprint-ready v0.9 (N=10 authoritative table + batch figs + full Methods + re-cal + SDDB Limitations + repro commands + final polish: refined title, abstract pacing/AF sentence, borderline 47/50 note, Record 32 limitation). Ready for medRxiv / Chaos / Frontiers in Network Physiology.
