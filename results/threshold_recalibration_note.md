# Light Threshold Re-calibration (CCTP v0.1 on SDDB RR)

Date: 2026-07-08
Based on: 10 records (30,31,32,35,36,38,45,47,50,51)

## Observed ranges on real Holter RR (bivariate proxy)

### |Δτ_s|
- median |Δτ_s| ≈ 0.026
- max |Δτ_s| ≈ 0.088 (records 38, 36, 50 strongest)
- Many small deltas (31, 45, 47) near 0.00–0.006

### excess3 (Nivel 3 proxy)
- Typical basal ~0.30–0.35
- Approach deltas: −0.032 to +0.042
- Strongest positive: 38 (+0.0415), 47 (+0.030), 36 (+0.008)
- high_level3_rate (current high_thresh=1.75) → always ~0.0 across records

### lambda (from |τ_s|)
- |τ_s(t)| values observed << synthetic 0.41 threshold (median |Δτ_s| ≈ 0.026, max ≈ 0.088 across 10 records)
- Current lambda ≈ 0 or constant on real data
- α3 weighting effectively flat → weighted and unweighted excess3 deltas very similar when both computed (e.g. rec 38: 0.04151 vs 0.04136)
- Weighted RECD completed on 6/10 records; the other 4 use unweighted excess3 (still fully valid and sign-concordant with τ_s)

## Proposed light re-cal (for RR physiology, not synthetic)

1. **theta3** (in phi3 / ordinal level detection)
   - Current: 0.10
   - Proposed: **0.08**
   - Rationale: RR series noisier than synthetic; slightly lower threshold captures more Level-3 structure without exploding false positives. Can be sensitivity-tested.

2. **high_thresh** (for high_level3_rate and "high excess3" classification)
   - Current: 1.75 (completely useless on real data)
   - Proposed: **0.65**
   - Rationale: Observed excess3 ~0.30 basal, up to ~0.38–0.42 in approach for strong transitions. 0.65 gives headroom above mean+2sd in good cases while producing non-zero rates in approach for 36/38/50 etc. Will make high_level3_rate informative.

3. **lambda activation / scaling**
   - Current: hard threshold ~0.41 on |τ_s|
   - Proposed options (document both):
     a) **Relative scaling** (preferred for robustness): λ(t) = |τ_s(t)| / max_|τ|_(record)   or use a soft ramp.
     b) **Lower absolute threshold**: activate at 0.08–0.10 (≈ median |Δτ|).
   - This will allow α3 ramp to actually vary on real pre-VF data.

## Suggested test after change
Re-run RECD (levels + weighted) on the 4–5 strongest records (36,38,50,32,51) with new params and compare:
- Δ high_level3_rate now >0 in approach for strong cases?
- frac_contrib3 dynamics more visible?
- Concordance with Δτ_s preserved?

## Code locations
- recd_ordinal_levels.py : theta3, high_thresh defaults + compute_lambda
- run_recd_on_rr.py and run_recd_weighted_on_rr.py : pass the params
- For preprint: report original + sensitivity table with the two sets.

## Impact on current results
Current conclusions (sign concordance of τ_s and excess3) are robust because they use continuous excess3, not the binary high_rate. Re-cal mainly improves interpretability and high_rate metric for future work.

Next step: implement a --theta3 / --high-thresh CLI or a recal_params.json and re-process strongest records.
