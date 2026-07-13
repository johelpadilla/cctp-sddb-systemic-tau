# Exploratory bake-off: OPC vs SDD vs frozen abs-z (SDDB + VFDB)

**Status:** Exploratory observational study only  
**Date:** 2026-07-13  
**Clinical / FDA / superiority claim:** **NONE**  
**Fusion of OPC and SDD:** **NONE** (options run completely separately)  
**Production abs-z retune:** **NONE** (frozen z≥2 sustained ≥3 on τ_s)

---

## 1. What was tested

Two formalized ordinal detectors from `docs/ORDINAL_ALARM_DETECTORS.md` and
`code/ordinal_detectors/`, applied independently to real cleaned RR:

| Option | Detector | Entry point | Default params used |
|--------|----------|-------------|---------------------|
| **1 — OPC** | Ordinal Persistence Collapse | `opc_detect` | \(L=8\), \(\theta_D=0.35\), \(\theta_R=5\) |
| **1 companion** | OPC with K-aware window | `opc_detect` | \(L=50\), same \(\theta_D,\theta_R\) (not abs-z retune; see §3) |
| **2 — SDD** | Symbolic Distribution Divergence (TV) | `sdd_detect` | \(L_c=50\), \(\theta_{\mathrm{TV}}=0.35\), \(\theta_S=1\), fixed basal, `mask_basal=True` |
| **Baseline** | Continuous abs-z on \(\tau_s\) | `detect_lead_time` | \(z\ge 2\), min consecutive 3 (frozen) |

Runner: `code/run_ordinal_exploratory_bakeoff.py`

### 1.1 Symbol stream (no pre-baked ordinal NPZ)

The repository **does not** store dedicated ordinal-symbol NPZ under `data/sddb/` or
`data/vfdb/` (those hold PhysioNet raw annotations). Symbols are **derived at runtime**
from cleaned RR:

1. Bivariate proxy \(X = [z(\mathrm{RR}),\, z(|\Delta\mathrm{RR}|)]\) via `build_bivariate_proxy`
2. Bandt–Pompe per channel: \(m=3\), delay \(=1\) via `generate_multivariate_symbols`
3. **Joint code** \(\sigma_t = s_{t,0}\cdot 6 + s_{t,1} \in \{0,\ldots,35\}\), **\(K=36\)**

Inventory: `results/ordinal_data_inventory.txt`

### 1.2 Records

| Cohort | Records | Windows |
|--------|---------|---------|
| **SDDB analytic** | 30, 31, 32, 35, 36, 38, 45, 47, 50, 51 | Holter basal/approach via `get_event_and_windows` |
| **SDDB extension** | 44 | Same Holter windows |
| **VFDB** | 22 cleaned RR (`data/rr_external/vfdb_*_clean.npz`) | `short_db_windows` (thirds of pre-event) |
| **NSRDB** (optional control flavor) | 8 cleaned RR (first ~6 h) | Post-basal alarm fraction only — **not** FAR/24h |

---

## 2. Per-option results (exploratory)

### 2.1 Option 1 — OPC (default \(L=8\))

| Cohort | n | Detected | Exploratory “sensitivity” | Median lead (h) |
|--------|---|----------|---------------------------|-----------------|
| SDDB (11) | 11 | 11 | **1.00** | ~7.0 |
| VFDB (22) | 22 | 22 | **1.00** | ~0.17 |
| All events | 33 | 33 | **1.00** | — |

**Basal / control behavior (critical):**

- SDDB mean basal alarm fraction ≈ **1.0** (alarms essentially everywhere after warm-up).
- NSRDB (8 controls): post-basal alarm fraction = **1.0** on every record.

**Qualitative:** With joint \(K=36\) and \(L=8\), maximum support diversity is
\(D_{\max}=L/K=8/36\approx 0.222 \le \theta_D=0.35\). Therefore \(\ell_t=1\) almost
always and after \(\theta_R=5\) consecutive windows **OPC always alarms**. This is a
**parameter–alphabet interaction**, not evidence of pre-VF specificity. Default
suggested ranges in the formalization were written with small alphabets (e.g. \(K=6\))
in mind.

### 2.2 Option 1 companion — OPC with \(L=50\) (K-aware scale)

Same \(\theta_D,\theta_R\); \(L\) raised so that \(D_t=\mathrm{n_{unique}}/K\) can
**exceed** \(\theta_D\) when many joint patterns appear. **Not** a production threshold
search and **not** a retune of abs-z.

| Cohort | n | Detected | Exploratory “sensitivity” | Median lead among detections (h) |
|--------|---|----------|---------------------------|----------------------------------|
| SDDB (11) | 11 | **6** | **0.55** | ~4.8 |
| VFDB (22) | 22 | **8** | **0.36** | ~0.05 |
| All events | 33 | 14 | 0.42 | — |

**Basal:** SDDB mean basal alarm fraction ≈ **0.005** (orders of magnitude lower than \(L=8\)).

**SDDB detections (L=50):** 31, 32, 36, 38, 44, 47 (misses 30, 35, 45, 50, 51).

**Qualitative:** When the window is long enough for diversity to be informative, OPC
behaves as a **true collapse detector**: it fires on a subset of events with low basal
activity, consistent with “few ordinal patterns sustained,” not continuous amplitude.

### 2.3 Option 2 — SDD (TV, independent of OPC)

| Cohort | n | Detected | Exploratory “sensitivity” | Median lead among detections (h) |
|--------|---|----------|---------------------------|----------------------------------|
| SDDB (11) | 11 | **11** | **1.00** | ~7.0 |
| VFDB (22) | 22 | **21** | **0.95** | ~0.14 |
| All events | 33 | 32 | 0.97 | — |

VFDB miss: record **419** (pre-event ≈ 0.011 h, stratum `too_short`).

**Basal / control:**

- SDDB mean diagnostic basal high-TV fraction ≈ **0.07** (primary alarms use
  `mask_basal=True`, so basal is not used as a search region).
- Mean max TV in approach (SDDB) ≈ **0.58** (well above \(\theta_{\mathrm{TV}}=0.35\)).
- NSRDB mean post-basal alarm fraction ≈ **0.11** (non-zero; structural drift vs fixed
  early basal on long Holters — exploratory false-alarm flavor, not FAR/24h).

**Qualitative:** SDD fires when the **current empirical law on \(\Sigma\)** diverges
from a fixed basal law (TV). It detects **distributional reorganization**, not
necessarily support collapse. On this data it is highly sensitive on SDDB and on
processable VFDB, with moderate control-side activity under a fixed early basal.

### 2.4 Frozen abs-z baseline (\(\tau_s\))

| Cohort | n | Detected | Exploratory “sensitivity” | Median lead among detections (h) |
|--------|---|----------|---------------------------|----------------------------------|
| SDDB (11) | 11 | **11** | **1.00** | ~5.9 |
| VFDB (22) | 22 | **19** | **0.86** | ~0.10 |
| All events | 33 | 30 | 0.91 | — |

VFDB misses under abs-z: **419**, **607**, **615** (all `too_short` pre-event strata).
Parameters unchanged from production discovery freeze.

---

## 3. Side-by-side comparison (observational)

### 3.1 SDDB analytic + 44

| Record | OPC \(L=8\) | OPC \(L=50\) | SDD (TV) | abs-z \(\tau_s\) |
|--------|:-----------:|:------------:|:--------:|:----------------:|
| 30 | yes | no | yes | yes |
| 31 | yes | yes | yes | yes |
| 32 | yes | yes | yes | yes |
| 35 | yes | no | yes | yes |
| 36 | yes | yes | yes | yes |
| 38 | yes | yes | yes | yes |
| 44 | yes | yes | yes | yes |
| 45 | yes | no | yes | yes |
| 47 | yes | yes | yes | yes |
| 50 | yes | no | yes | yes |
| 51 | yes | no | yes | yes |

Full numeric table: `results/ordinal_opc_sdd_absz_comparison.csv`

### 3.2 Notable differences (not superiority)

1. **OPC \(L=8\) + \(K=36\) is not a usable specificity probe** — it saturates by
   construction (\(D\le L/K\le\theta_D\)). Any “100% sensitivity” here is **not**
   comparable to abs-z performance.
2. **OPC \(L=50\)** is more selective than abs-z on SDDB (6/11 vs 11/11) with much
   cleaner basal behavior; misses several Holters that abs-z and SDD catch. It is a
   **different predicate** (collapse + persistence), not a drop-in replacement.
3. **SDD** matches abs-z sensitivity on SDDB (11/11) and is slightly higher on all
   VFDB (21/22 vs 19/22), but NSRDB control-side alarm fractions (~0.11 of endpoints)
   show it is **not** free of false structure under fixed basal on long Holters.
4. **Lead times** for SDD and OPC \(L=8\) on SDDB often start near the end of basal /
   start of search (very long leads ~7 h) — early first-crossing of a loose ordinal
   predicate, not necessarily a tight pre-VF localization. abs-z median leads on SDDB
   are somewhat shorter (~5.9 h) but still early relative to event.
5. Options remain **independent**: no AND/OR fusion was defined or scored.

### 3.3 VFDB short_15_60min subset (n=11, Phase-1 style)

| Detector | Detected |
|----------|----------|
| OPC \(L=8\) | 11/11 (saturated) |
| OPC \(L=50\) | 5/11 |
| SDD | 11/11 |
| abs-z \(\tau_s\) | 11/11 |

Short pre-event records remain difficult to interpret for any detector; many VFDB
files are stratum `too_short` (&lt;15 min pre-event).

---

## 4. Artifacts

| File | Content |
|------|---------|
| `results/ordinal_data_inventory.txt` | RR inventory; note on derived symbols |
| `results/ordinal_opc_per_record.csv` | OPC \(L=8\) per record |
| `results/ordinal_opc_L50_per_record.csv` | OPC \(L=50\) companion |
| `results/ordinal_sdd_per_record.csv` | SDD per record |
| `results/ordinal_absz_per_record.csv` | Frozen abs-z on \(\tau_s\) |
| `results/ordinal_opc_sdd_absz_comparison.csv` | Side-by-side |
| `results/ordinal_opc_nsrdb_controls.csv` | Control flavor OPC |
| `results/ordinal_sdd_nsrdb_controls.csv` | Control flavor SDD |
| `results/ordinal_exploratory_summary.json` | Aggregate JSON |
| `code/run_ordinal_exploratory_bakeoff.py` | Reproducible runner |

Reproduce:

```bash
python3 code/run_ordinal_exploratory_bakeoff.py
```

---

## 5. Honest conclusions (exploratory)

1. **Native ordinal alarms can be run end-to-end on real SDDB/VFDB RR** without
   moments of the continuous series inside OPC/SDD.
2. **Default OPC params require alphabet-aware \(L\)** (or \(\theta_D\) re-expressed
   as absolute support size). With joint \(K=36\), \(L=8\) is not informative.
3. With a longer observation window, **OPC looks like a collapse/persistence rule**:
   lower event hit rate than abs-z, low basal activity — a different operating point.
4. **SDD (TV)** is highly sensitive on these event sets and tracks distributional
   shift relative to fixed basal; control Holters show non-zero post-basal alarms —
   needs careful basal policy before any claim of specificity.
5. **abs-z on \(\tau_s\) remains the frozen empirical baseline** for the pilot; OPC/SDD
   are **methodological alternatives**, not replacements.
6. **No fusion, no clinical claim, no production retune** in this bake-off.

---

## 6. What this does *not* claim

- Superiority of OPC or SDD over abs-z
- Clinical utility, deployability, or FDA readiness
- Optimized ordinal thresholds (no grid search / validation-set retune)
- A fused OPC∧SDD detector
- Phase-2 institutional specificity (public NSRDB control flavor only)
