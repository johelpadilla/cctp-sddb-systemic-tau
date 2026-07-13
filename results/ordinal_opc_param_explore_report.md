# OPC parameter exploration (exploratory)

Generated: 2026-07-13T18:42:06.173999+00:00

## Framing

- Exploratory only — not clinical optimization, not S5, not FDA, not superiority.
- Abs-z remains frozen (z=2, min_run=3); not retuned here.
- No fusion with SDD.

## Parameter grid (justified modest)

- **L** ∈ [40, 50, 60, 70] — Around L=50 companion used for K=36 so min(L,K)/K can exceed θ_D (need L > θ_D·K ≈ 12.6). Explore ±10–20: {40,50,60,70}.
- **θ_D** ∈ [0.3, 0.35, 0.4] — Around 0.35 baseline diversity collapse threshold. Lower (0.30) is stricter collapse → fewer alarms; higher (0.40) is looser → more sensitive.
- **θ_R** ∈ [4, 5, 6] — Around 5 consecutive low-diversity windows. Lower (4) eases alarm (more sens / more FAR); higher (6) requires longer collapse.
- Modest 4×3×3=36 cell product — not exhaustive grid search or nested CV. Exploratory only; abs-z remains frozen; no fusion with SDD.

- Cells evaluated: **36**

## Methodology

- Sensitivity: first OPC alarm in post-basal pre-event window (SDDB n=11, VFDB n=22).
- FAR: binary alarm episodes, refractory 0.5 h, Phase-2 basal/search, control cap 12 h, NSRDB n=18.
- FAR formula: total_episodes / total_search_hours × 24.
- Encoding: joint bivariate Bandt–Pompe m=3, K=36 (unchanged).

## Baseline (L=50, θ_D=0.35, θ_R=5)

- sens_sddb = 0.5455 (6/11)
- sens_vfdb = 0.3636 (8/22)
- sens_all = 0.4242 (14/33)
- FAR = 3.7334 /24h (28 ep / 180.00 h search)

## Recommendation

- **Decision:** `keep_baseline`
- **Reason:** No cell increased all-event sensitivity while keeping FAR within exploratory slack (≤ max(1.5× baseline, baseline+2.0)). Some cells raised sensitivity only with substantially higher FAR.
- **Adopt params (if any):** `{'L': 50, 'theta_D': 0.35, 'theta_R': 5}`

## Best candidates (higher sens, FAR within slack)

_None._

## Qualitative parameter effects

- L=40: mean sens_all=0.545, mean FAR=14.00 | L=50: mean sens_all=0.418, mean FAR=6.27 | L=60: mean sens_all=0.313, mean FAR=2.92 | L=70: mean sens_all=0.222, mean FAR=1.44
- θ_D=0.30: mean sens_all=0.179, mean FAR=0.40 | θ_D=0.35: mean sens_all=0.374, mean FAR=3.94 | θ_D=0.40: mean sens_all=0.571, mean FAR=14.12
- θ_R=4: mean sens_all=0.379, mean FAR=6.52 | θ_R=5: mean sens_all=0.376, mean FAR=6.16 | θ_R=6: mean sens_all=0.369, mean FAR=5.79
- L effect (high vs low): Δmean_sens=-0.323, Δmean_FAR=-12.56 (larger L averages diversity over more symbols — tends to dampen brief collapses).
- θ_D effect (high vs low): Δmean_sens=+0.391, Δmean_FAR=+13.72 (higher θ_D is a looser diversity threshold → easier low-div declaration).
- θ_R effect (high vs low): Δmean_sens=-0.010, Δmean_FAR=-0.73 (higher θ_R requires longer consecutive collapse before alarm).
- Framing: descriptive grid margins only — not causal clinical effects; n is small.

## Full grid (sorted by sens_all desc, then FAR asc)

| is_baseline | L | θ_D | θ_R | sens_sddb | sens_vfdb | sens_all | FAR | class |
|---:|---:|---:|---:|---:|---:|---:|---:|---|
| 0 | 40 | 0.40 | 6 | 1.000 | 0.682 | 0.788 | 28.134 | higher_sens_far_worse |
| 0 | 40 | 0.40 | 5 | 1.000 | 0.682 | 0.788 | 29.734 | higher_sens_far_worse |
| 0 | 40 | 0.40 | 4 | 1.000 | 0.682 | 0.788 | 30.800 | higher_sens_far_worse |
| 0 | 50 | 0.40 | 6 | 0.909 | 0.500 | 0.636 | 14.534 | higher_sens_far_worse |
| 0 | 50 | 0.40 | 5 | 0.909 | 0.500 | 0.636 | 14.934 | higher_sens_far_worse |
| 0 | 50 | 0.40 | 4 | 0.909 | 0.500 | 0.636 | 15.867 | higher_sens_far_worse |
| 0 | 40 | 0.35 | 6 | 0.727 | 0.455 | 0.545 | 9.733 | higher_sens_far_worse |
| 0 | 40 | 0.35 | 5 | 0.727 | 0.455 | 0.545 | 10.800 | higher_sens_far_worse |
| 0 | 40 | 0.35 | 4 | 0.727 | 0.455 | 0.545 | 12.133 | higher_sens_far_worse |
| 0 | 60 | 0.40 | 4 | 0.727 | 0.409 | 0.515 | 8.000 | higher_sens_far_worse |
| 0 | 60 | 0.40 | 6 | 0.636 | 0.409 | 0.485 | 7.600 | higher_sens_far_worse |
| 0 | 60 | 0.40 | 5 | 0.636 | 0.409 | 0.485 | 7.867 | higher_sens_far_worse |
| 1 | 50 | 0.35 | 5 | 0.545 | 0.364 | 0.424 | 3.733 | baseline |
| 0 | 50 | 0.35 | 4 | 0.545 | 0.364 | 0.424 | 3.867 | same_sens_far_ok |
| 0 | 50 | 0.35 | 6 | 0.545 | 0.318 | 0.394 | 3.333 | lower_sens_lower_far |
| 0 | 70 | 0.40 | 6 | 0.545 | 0.273 | 0.364 | 3.733 | other |
| 0 | 70 | 0.40 | 5 | 0.545 | 0.273 | 0.364 | 4.000 | lower_sens_higher_far |
| 0 | 70 | 0.40 | 4 | 0.545 | 0.273 | 0.364 | 4.267 | lower_sens_higher_far |
| 0 | 60 | 0.35 | 5 | 0.455 | 0.273 | 0.333 | 0.933 | lower_sens_lower_far |
| 0 | 60 | 0.35 | 4 | 0.455 | 0.273 | 0.333 | 1.067 | lower_sens_lower_far |
| 0 | 60 | 0.35 | 6 | 0.364 | 0.273 | 0.303 | 0.800 | lower_sens_lower_far |
| 0 | 40 | 0.30 | 6 | 0.455 | 0.227 | 0.303 | 1.333 | lower_sens_lower_far |
| 0 | 40 | 0.30 | 5 | 0.455 | 0.227 | 0.303 | 1.600 | lower_sens_lower_far |
| 0 | 40 | 0.30 | 4 | 0.455 | 0.227 | 0.303 | 1.733 | lower_sens_lower_far |
| 0 | 50 | 0.30 | 5 | 0.364 | 0.136 | 0.212 | 0.000 | lower_sens_lower_far |
| 0 | 50 | 0.30 | 4 | 0.364 | 0.136 | 0.212 | 0.133 | lower_sens_lower_far |
| 0 | 70 | 0.35 | 5 | 0.364 | 0.136 | 0.212 | 0.267 | lower_sens_lower_far |
| 0 | 70 | 0.35 | 6 | 0.364 | 0.136 | 0.212 | 0.267 | lower_sens_lower_far |
| 0 | 70 | 0.35 | 4 | 0.364 | 0.136 | 0.212 | 0.400 | lower_sens_lower_far |
| 0 | 50 | 0.30 | 6 | 0.273 | 0.136 | 0.182 | 0.000 | lower_sens_lower_far |
| 0 | 60 | 0.30 | 4 | 0.273 | 0.045 | 0.121 | 0.000 | lower_sens_lower_far |
| 0 | 60 | 0.30 | 5 | 0.273 | 0.045 | 0.121 | 0.000 | lower_sens_lower_far |
| 0 | 60 | 0.30 | 6 | 0.273 | 0.045 | 0.121 | 0.000 | lower_sens_lower_far |
| 0 | 70 | 0.30 | 4 | 0.273 | 0.000 | 0.091 | 0.000 | lower_sens_lower_far |
| 0 | 70 | 0.30 | 5 | 0.273 | 0.000 | 0.091 | 0.000 | lower_sens_lower_far |
| 0 | 70 | 0.30 | 6 | 0.273 | 0.000 | 0.091 | 0.000 | lower_sens_lower_far |

