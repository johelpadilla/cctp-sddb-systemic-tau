# Exploratory light cascade: SDD → OPC L=50 (±5 min confirm)

**Generated:** 2026-07-13T18:11:19.799998+00:00

**Status:** Exploratory only. No clinical claim, no S5 claim, no superiority claim, no production threshold change.

## Cascade rule

1. **Candidate filter:** SDD (L_c=50, θ_TV=0.35, θ_S=1).
2. **Confirmation:** An SDD alarm at time *t* becomes a cascade alarm only if **OPC L=50** (θ_D=0.35, θ_R=5) also alarms within a **closed window [t−W, t+W]** with **W = ±5 minutes** (0.083333 h).
3. Cascade marks only confirmed **SDD-on** samples (OPC-only never starts a cascade alarm).
4. **Causality (events):** only OPC with `t_opc < event_hr` may confirm; `decision_time = max(t_SDD, t_OPC_confirm) < event_hr` (no post-event look-ahead).

### Window justification

±5 min is a local co-occurrence confirm: short VFDB horizons and multi-hour SDDB leads both only need same short-term instability epoch; not a second independent long-horizon detector. Causality forbids post-event OPC from confirming pre-event SDD.

## Head-to-head comparison

| Detector | SDDB sens (n=11) | VFDB sens (n=22) | All sens (n=33) | NSRDB FAR (/24h) | Episodes |
|----------|------------------|------------------|-----------------|------------------|----------|
| **Cascade SDD→OPC** | 0.545 (6/11) | 0.318 (7/22) | **0.394** (13/33) | **3.867** | 29 |
| OPC L=50 alone | 0.545 (6/11) | 0.364 (8/22) | 0.424 (14/33) | 3.733 | 28 |
| SDD alone | 1.000 (11/11) | 0.955 (21/22) | 0.970 (32/33) | 46.267 | 347 |
| abs-z τ_s frozen | 1.000 (11/11) | 0.864 (19/22) | 0.909 (30/33) | 33.734 | 253 |

## Events gained / lost (cascade vs singletons)

| Comparison | Gained | Lost | Both on | Both off |
|------------|--------|------|---------|----------|
| vs OPC L=50 | 0 | 1 | 13 | 19 |
| vs SDD | 0 | 19 | 13 | 1 |
| vs abs-z | 0 | 17 | 13 | 3 |

## Observations

- Cascade all-event sensitivity = 0.394 (OPC L=50 0.424, SDD 0.970, abs-z 0.909).
- Cascade pooled NSRDB FAR = 3.867/24h (OPC 3.733, SDD 46.267, abs-z 33.734).
- Cascade FAR ≤ SDD FAR on pooled NSRDB, as expected under OPC confirmation (confirmation can only drop or keep SDD episodes, never invent SDD-less alarms).
- Events vs SDD: gained 0, lost 19 (lost expected when SDD fires without local OPC confirmation).
- Events vs OPC L=50: gained 0, lost 1 (gains require SDD candidate co-located with OPC; losses when OPC alone was enough).
- Events vs abs-z: gained 0, lost 17.
- Cascade is SDD-first: it cannot detect events that SDD misses; sensitivity ceiling is SDD's sensitivity (minus unconfirmed SDD hits).
- NSRDB is healthy Holter, not device-matched to VFDB/SDDB telemetry — FAR is an interim public upper-bound estimate only.
- Exploratory only: no clinical claim, no S5 claim, no superiority claim, no production threshold change.

## Limitations

- Cascade cannot recover events that SDD misses (SDD is hard filter).
- Confirmation may drop true SDD positives when OPC does not co-fire locally pre-event.
- Causal rule: post-event OPC cannot confirm (no look-ahead); sensitivity can drop vs non-causal co-occurrence.
- VFDB short records: ±5 min window can miss staggered OPC/SDD peaks.
- NSRDB device mismatch; FAR not ICU/telemetry-matched.
- Timebases: abs-z on strided τ_s vs symbol endpoints for ordinal/cascade.
- Fixed params only; no threshold or window grid search (avoids fishing).
- Exploratory fusion arm — singleton experiments remain fusion=False.

## Executive summary

Light cascade SDD→OPC (±5 min): all-event sens=0.394 (SDDB 0.545, VFDB 0.318); NSRDB FAR=3.867/24h vs OPC 3.733, SDD 46.267, abs-z 33.734. Recommendation: low_priority. Exploratory only — no clinical claims.

## Recommendation

**low_priority:** Cascade does not improve on OPC L=50 in both sens and FAR under these fixed params; further cascade work is low priority unless a new confirm rule is justified a priori (avoid post-hoc window fishing).

## Reproduce

```bash
cd Investigaciones/Cardiac_CCTP_Pilot
python3 code/run_ordinal_cascade_fusion.py --write-doc
python3 -m pytest tests/test_ordinal_cascade_fusion.py -q
```

## Artifacts

- `results/ordinal_cascade_comparison.csv`
- `results/ordinal_cascade_per_record.csv`
- `results/ordinal_cascade_nsrdb_far_per_record.csv`
- `results/ordinal_cascade_gain_loss.csv`
- `results/ordinal_cascade_summary.json`
- Pure merger: `code/ordinal_detectors/cascade_fusion.py`
- Entry: `code/run_ordinal_cascade_fusion.py`

