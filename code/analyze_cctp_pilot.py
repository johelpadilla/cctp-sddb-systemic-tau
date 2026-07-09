#!/usr/bin/env python3
"""
Cardiac Critical Transitions Protocol (CCTP) v0.1 - Multi-record pilot on PhysioNet SDDB

Surgical implementation:
- RR from SDDB records (default 35 best long pre-VF; 30 for intermediate event)
- Classic EWS: rolling variance + lag-1 AR(1)
- Systemic Tau via systemictau.compute_taus on minimal proxy (z(RR), z(|dRR|))
- Many diagnostic + comparative figures (per record)
- Basal vs Approach statistical comparison + JSON summary
- Reuses systemictau exactly; supports --record for easy multi-record runs

Usage:
  python3 code/analyze_cctp_pilot.py --record 35
  python3 code/analyze_cctp_pilot.py --record 30
"""

import os
import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import warnings
warnings.filterwarnings("ignore")

parser = argparse.ArgumentParser()
parser.add_argument("--record", default="35", help="SDDB record id (e.g. 35, 30)")
args = parser.parse_args()
RECORD = args.record
print(f"Running CCTP pilot for record {RECORD}")

# Make systemictau importable (pip package, or local SYSTEMICTAU_SRC / monorepo)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _bootstrap import import_systemictau_core
compute_taus, systemic_tau, HAS_SYSTEMICTAU = import_systemictau_core()
if not HAS_SYSTEMICTAU:
    print("Warning: could not import systemictau.core. Falling back to manual tau.")
    print("Install with: pip install systemictau  (or set SYSTEMICTAU_SRC)")

# Paths (per-record outputs for clean multi-record work)
BASE = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
DATA_DIR = os.path.join(BASE, "data")
FIG_DIR = os.path.join(BASE, "figures", RECORD)
RES_DIR = os.path.join(BASE, "results")
os.makedirs(FIG_DIR, exist_ok=True)
os.makedirs(RES_DIR, exist_ok=True)

# Load cleaned RR (rr_{RECORD}_clean.npz must exist)
rr_npz = os.path.join(DATA_DIR, f"rr_{RECORD}_clean.npz")
if not os.path.exists(rr_npz):
    raise FileNotFoundError(f"Missing {rr_npz}. Generate it first (see extraction logic).")
d = np.load(rr_npz)
rr = d["rr_ms"].astype(float)
t_sec = d["t_sec"].astype(float)
vfon_sec = float(d["vfon_sec"])
total_h = float(d["total_hours"])

print("="*60)
print(f"CCTP Pilot - Record {RECORD}")
print(f"RR length: {len(rr):,} beats | total span: {total_h:.2f} h")
print(f"VF onset marker (comment): {vfon_sec/3600:.2f} h")
print("="*60)

# === Event alignment ===
# For record 35 (event near end of file) we historically use t.max() to avoid
# marker placed after actual data. For records with clear intermediate event
# (e.g. 30) we use the vfon comment time.
event_idx = np.argmin(np.abs(t_sec - vfon_sec))
t_hr = t_sec / 3600.0
vfon_hr = vfon_sec / 3600.0
total_h = float(t_hr.max())
if RECORD == "35" or (vfon_hr > total_h - 0.8):
    event_hr = total_h
else:
    event_hr = vfon_hr
print(f"Event aligned at beat index ~{event_idx} / {len(rr)} (t={t_hr[event_idx]:.2f}h)")
print(f"Using event_hr={event_hr:.2f}h for windows/markers (RECORD={RECORD})")

# === Define windows (in #beats) ===
# For cardiac, use clinically meaningful windows: ~5-15 min
# avg RR ~850ms => ~70 beats/min => W=350 ~5min, W=700~10min
W_EWS = 501   # ~7 min local window (odd for centering)
W_TAU = 101   # smaller for tau to keep resolution (~1.4 min)
STRIDE = 5    # denser for much better visual in plots (n=100k is small)

def rolling_var(x, w, stride=1):
    """Rolling variance (population)"""
    n = len(x)
    out = np.full(n, np.nan)
    for i in range(w-1, n, stride):
        out[i] = np.var(x[i-w+1:i+1])
    return out

def rolling_ar1(x, w, stride=1):
    """Rolling lag-1 autocorrelation (Pearson)"""
    n = len(x)
    out = np.full(n, np.nan)
    for i in range(w-1, n, stride):
        seg = x[i-w+1:i+1]
        if len(seg) > 2 and np.std(seg) > 1e-9:
            r = np.corrcoef(seg[:-1], seg[1:])[0, 1]
            out[i] = r if not np.isnan(r) else np.nan
    return out

def rolling_mean(x, w, stride=1):
    n = len(x)
    out = np.full(n, np.nan)
    for i in range(w-1, n, stride):
        out[i] = np.mean(x[i-w+1:i+1])
    return out

# Compute classic EWS
print("Computing rolling EWS metrics (var, ar1, mean)...")
var = rolling_var(rr, W_EWS, STRIDE)
ar1 = rolling_ar1(rr, W_EWS, STRIDE)
rmean = rolling_mean(rr, W_EWS, STRIDE)

# === Systemic Tau (reusing framework) ===
print("Computing Systemic Tau (reusing compute_taus)...")
# Surgical bivariate proxy for univariate RR:
#   Module 0: standardized RR
#   Module 1: standardized |first difference| (captures local irregularity)
drr = np.abs(np.diff(rr, prepend=rr[0]))
rr_z = (rr - np.mean(rr)) / (np.std(rr) + 1e-12)
drr_z = (drr - np.mean(drr)) / (np.std(drr) + 1e-12)
X = np.column_stack([rr_z, drr_z])   # (T, 2)

if HAS_SYSTEMICTAU:
    taus_global, taus_per = compute_taus(X, window_size=W_TAU, stride=STRIDE)
else:
    # Fallback manual (simple mean pairwise)
    taus_global = np.full(len(rr), np.nan)
    for i in range(W_TAU-1, len(rr), STRIDE):
        win = X[i-W_TAU+1:i+1]
        # rank corr approx using pearson on ranks for speed
        r = np.corrcoef(stats.rankdata(win[:,0]), stats.rankdata(win[:,1]))[0,1]
        taus_global[i] = r if not np.isnan(r) else 0.0
    taus_per = np.zeros((len(rr), 2))

tau_s = taus_global.copy()

# Also compute a "pure" univariate proxy: rolling kendall of RR vs lagged RR inside window? But ar1 is close.
# For comparison, compute simple rolling spearman-like on consecutive inside window already covered by ar1.

# t_hr and event_hr already computed above (parameterized per RECORD).
# For metrics, the value at position i represents window ending at i (causal)

# Downsample for clean plots
def downsample_for_plot(y, t, stride):
    idx = np.arange(0, len(y), stride)
    return t[idx], y[idx]

def get_valid_plot_data(t, y):
    """Return only non-NaN points (critical for strided rolling outputs)."""
    valid = ~np.isnan(y)
    return t[valid], y[valid]

# === Basal / Approach windows (parameterized) ===
# Always define relative to event_hr for the chosen record.
approach_start = max(0.0, event_hr - 3.0)
approach_end = event_hr

# Robust basal selection (preserve exact 30/35 windows; sensible default for others)
if RECORD == "35":
    basal_start, basal_end = 6.0, 16.0
elif RECORD == "30":
    basal_start, basal_end = 0.5, 3.5
else:
    b_end = max(3.5, approach_start - 4.0)
    basal_start = max(0.5, b_end - 3.0)
    basal_end = b_end

print(f"Windows -> basal: [{basal_start}, {basal_end}]h | approach: [{approach_start:.2f}, {approach_end:.2f}]h")

# === Figure batch: Raw + EWS indicators (many panels) ===
print("Generating figures...")

# 1. Full RR with event + shaded approach
fig, ax = plt.subplots(1,1, figsize=(14, 4.2))
ax.plot(t_hr, rr, lw=0.35, color="#1f77b4", alpha=0.65, label="RR (cleaned)")
ax.axvline(event_hr, color="#d62728", lw=2.5, ls="--", label=f"VF onset / event ({event_hr:.2f}h)")
ax.set_xlabel("Time from recording start (hours)")
ax.set_ylabel("RR interval (ms)")
ax.set_title(f"Record {RECORD} (SDDB) — Full cleaned R-R series + event ground truth")
ax.legend(loc="upper right")
ax.grid(True, alpha=0.25)
ax.axvspan(approach_start, event_hr, alpha=0.15, color="red", label="Approach window (3h)")
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "05_rr_full_with_approach.png"), dpi=160)
plt.close()

# 2. Multi-panel EWS
fig, axes = plt.subplots(4, 1, figsize=(13, 10.8), sharex=True)

# RR (downsampled for viz)
t_ds, rr_ds = downsample_for_plot(rr, t_hr, 5)
axes[0].plot(t_ds, rr_ds, lw=0.4, color="#1f77b4")
axes[0].axvline(event_hr, color="#d62728", lw=2, ls="--")
axes[0].axvspan(approach_start, event_hr, alpha=0.12, color="red")
axes[0].set_ylabel("RR (ms)")
axes[0].set_title("R-R interval (subsampled)")

# Variance
t_v, v_ds = get_valid_plot_data(t_hr, var)
axes[1].plot(t_v, v_ds, lw=1.0, color="#ff7f0e")
axes[1].axvline(event_hr, color="#d62728", lw=2, ls="--")
axes[1].axvspan(approach_start, event_hr, alpha=0.12, color="red")
axes[1].set_ylabel("Rolling Var\n(W=501 beats ~7min)")

# AR(1)
t_a, a_ds = get_valid_plot_data(t_hr, ar1)
axes[2].plot(t_a, a_ds, lw=1.0, color="#2ca02c")
axes[2].axvline(event_hr, color="#d62728", lw=2, ls="--")
axes[2].axvspan(approach_start, event_hr, alpha=0.12, color="red")
axes[2].set_ylabel("Lag-1 AR(1)\n(W=501)")
axes[2].set_ylim(-0.1, 1.05)

# Systemic Tau
t_tau, tau_ds = get_valid_plot_data(t_hr, tau_s)
axes[3].plot(t_tau, tau_ds, lw=1.1, color="#9467bd")
axes[3].axvline(event_hr, color="#d62728", lw=2, ls="--")
axes[3].axvspan(approach_start, event_hr, alpha=0.12, color="red")
axes[3].set_ylabel("Systemic Tau τ_s\n(W=101, proxy RR+|ΔRR|)")
axes[3].set_xlabel("Time (hours)")
axes[3].set_ylim(-0.1, 1.05)

for ax in axes:
    ax.grid(True, alpha=0.25)

fig.suptitle(f"Early Warning Signals (EWS) + Systemic Tau — Record {RECORD}\nVertical red = event marker | Red band = approach window",
             fontsize=11, y=0.975)
fig.tight_layout(rect=[0, 0, 1, 0.92])
fig.savefig(os.path.join(FIG_DIR, "06_ews_panels.png"), dpi=160, bbox_inches="tight", pad_inches=0.25)
plt.close()

# 3. Zoomed on approach (last 6h): all indicators
pre6 = max(0, event_hr - 6)
mask = (t_hr >= pre6)

fig, axes = plt.subplots(4, 1, figsize=(11, 9), sharex=True)
axes[0].plot(t_hr[mask], rr[mask], lw=0.5, color="#1f77b4")
axes[0].axvline(event_hr, color="#d62728", lw=2.5, ls="--")
axes[0].set_ylabel("RR (ms)")
axes[0].set_title(f"Zoom: last 6h before event (Record {RECORD}) + indicators")

for i, (y, lab, col) in enumerate([
    (var, "Var (W≈7min)", "#ff7f0e"),
    (ar1, "AR(1)", "#2ca02c"),
    (tau_s, "τ_s (Systemic)", "#9467bd")
]):
    t_y, y_valid = get_valid_plot_data(t_hr, y)
    m = (t_y >= pre6)
    axes[i+1].plot(t_y[m], y_valid[m], lw=1.2, color=col)
    axes[i+1].axvline(event_hr, color="#d62728", lw=2, ls="--")
    axes[i+1].set_ylabel(lab)
axes[3].set_xlabel("Time (h)")
for ax in axes: ax.grid(alpha=0.25)
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "07_ews_zoom_6h.png"), dpi=160)
plt.close()

# 4. Basal vs Approach comparison (boxplots + stats)
# Use the parameterized windows defined earlier
def get_valid(y, m):
    v = y[m]
    return v[~np.isnan(v)]

basal_mask = (t_hr >= basal_start) & (t_hr <= basal_end)
approach_mask = (t_hr >= approach_start) & (t_hr < approach_end)

basal_var = get_valid(var, basal_mask)
app_var = get_valid(var, approach_mask)
basal_ar1 = get_valid(ar1, basal_mask)
app_ar1 = get_valid(ar1, approach_mask)
basal_tau = get_valid(tau_s, basal_mask)
app_tau = get_valid(tau_s, approach_mask)

print(f"\n=== Basal ({basal_start:.1f}-{basal_end:.1f}h) vs Approach (last ~3h pre-event) ===")
for name, b, a in [("Var", basal_var, app_var), ("AR1", basal_ar1, app_ar1), ("Tau_s", basal_tau, app_tau)]:
    mb, ma = np.nanmean(b), np.nanmean(a)
    print(f"{name:6s} | basal mean={mb:.5f}  approach mean={ma:.5f}  delta={ma-mb:+.5f}  (n_b={len(b)}, n_a={len(a)})")

# Boxplot figure
fig, axes = plt.subplots(1, 3, figsize=(12, 4))
basal_label = f"Basal\n({basal_start:.1f}-{basal_end:.1f}h)"
app_label = f"Approach\n(last ~3h pre)"
labels = [basal_label, app_label]
for ax, b, a, title in zip(axes, 
                           [basal_var, basal_ar1, basal_tau],
                           [app_var, app_ar1, app_tau],
                           ["Rolling Variance", "AR(1) autocorr", "Systemic Tau τ_s"]):
    ax.boxplot([b, a], tick_labels=labels, patch_artist=True,
               boxprops=dict(facecolor="#cce5ff"), medianprops=dict(color="red", lw=2))
    ax.set_title(title)
    ax.grid(True, alpha=0.3, axis="y")
fig.suptitle(f"Distribution shift: Basal vs Approach (Record {RECORD})", fontsize=11)
fig.tight_layout(rect=[0, 0, 1, 0.94])
fig.savefig(os.path.join(FIG_DIR, "08_basal_vs_approach_boxplots.png"), dpi=150, bbox_inches="tight", pad_inches=0.2)
plt.close()

# 5. Scatter / trajectory of (AR1, Var) colored by time-to-event
fig, ax = plt.subplots(figsize=(7,6))
# Sample only valid (non-NaN) positions within each regime
n_samp = 400
valid_b = np.where(basal_mask & ~np.isnan(var))[0]
valid_a = np.where(approach_mask & ~np.isnan(var))[0]
idx_b = valid_b[::max(1, len(valid_b)//n_samp)]
idx_a = valid_a[::max(1, len(valid_a)//n_samp)]
ax.scatter(ar1[idx_b], var[idx_b], c="#1f77b4", s=8, alpha=0.5, label="Basal")
ax.scatter(ar1[idx_a], var[idx_a], c="#d62728", s=12, alpha=0.7, label="Approach (pre-event)")
ax.set_xlabel("AR(1)")
ax.set_ylabel("Variance")
ax.set_title(f"Phase plane: AR(1) vs Variance (Record {RECORD})\n(color = regime relative to event)")
ax.legend()
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "09_phase_plane_ar1_var.png"), dpi=150)
plt.close()

# 6. Delta RR (irregularity) over time + tau relation
fig, ax = plt.subplots(2,1, figsize=(12,6), sharex=True)
ax[0].plot(t_hr[::5], drr[::5], lw=0.3, color="#e377c2", alpha=0.6)
ax[0].axvline(event_hr, color="r", lw=2, ls="--")
ax[0].axvspan(approach_start, event_hr, alpha=0.1, color="red")
ax[0].set_ylabel("|ΔRR| (ms)")
ax[0].set_title(f"Local beat-to-beat irregularity |ΔRR| (Record {RECORD})")

t_tau_valid, tau_valid = get_valid_plot_data(t_hr, tau_s)
ax[1].plot(t_tau_valid, tau_valid, lw=1.1, color="#9467bd")
ax[1].axvline(event_hr, color="r", lw=2, ls="--")
ax[1].set_ylabel("τ_s")
ax[1].set_xlabel("Time (h)")
for a in ax: a.grid(alpha=0.25)
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "10_irregularity_and_tau.png"), dpi=150)
plt.close()

# 7. Detailed last 90 min before event (high res)
pre90m = max(0, event_hr - 1.5)
m90 = t_hr >= pre90m
fig, axes = plt.subplots(3,1, figsize=(10,7), sharex=True)
axes[0].plot(t_hr[m90], rr[m90], lw=0.7, color="#1f77b4")
axes[0].axvline(event_hr, color="r", lw=2, ls="--")
axes[0].set_ylabel("RR (ms)")

t_var, var_valid = get_valid_plot_data(t_hr, var)
m90v = t_var >= pre90m
axes[1].plot(t_var[m90v], var_valid[m90v], lw=1.2, color="#ff7f0e")
axes[1].axvline(event_hr, color="r", lw=2, ls="--")
axes[1].set_ylabel("Var")

t_ar, ar_valid = get_valid_plot_data(t_hr, ar1)
m90a = t_ar >= pre90m
axes[2].plot(t_ar[m90a], ar_valid[m90a], lw=1.2, color="#2ca02c")
axes[2].axvline(event_hr, color="r", lw=2, ls="--")
axes[2].set_ylabel("AR1")
axes[2].set_xlabel("Time (h)")
fig.suptitle(f"Last ~90 min before event — high-resolution EWS (Record {RECORD})")
for a in axes: a.grid(alpha=0.25)
fig.tight_layout(rect=[0, 0, 1, 0.94])
fig.savefig(os.path.join(FIG_DIR, "11_last90min_detail.png"), dpi=160, bbox_inches="tight", pad_inches=0.2)
plt.close()

# === Additional diagnostic figures (rich visuals as per protocol emphasis) ===

# 12. RR histogram: basal vs approach
fig, ax = plt.subplots(figsize=(8,4.5))
basal_rr = rr[basal_mask]
app_rr = rr[approach_mask]
ax.hist(basal_rr, bins=40, alpha=0.55, label=f"Basal ({basal_start:.1f}-{basal_end:.1f}h)", color="#1f77b4", density=True)
ax.hist(app_rr, bins=40, alpha=0.55, label=f"Approach (pre-event)", color="#d62728", density=True)
ax.set_xlabel("RR (ms)")
ax.set_ylabel("Density")
ax.set_title(f"RR interval distribution (Record {RECORD})")
ax.legend()
ax.grid(True, alpha=0.2, axis="y")
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "12_rr_hist.png"), dpi=150)
plt.close()

# 13. Normalized overlay of EWS (z-score) in approach window (last ~8h pre)
pre8 = max(0, event_hr - 8)
m8 = (t_hr >= pre8) & (t_hr <= event_hr)
t8 = t_hr[m8]
def z(x):
    x = np.asarray(x)
    return (x - np.nanmean(x)) / (np.nanstd(x) + 1e-12)
fig, ax = plt.subplots(figsize=(11, 4.5))
# sample the strided metrics for the window
t_v8, v8 = get_valid_plot_data(t_hr[m8], var[m8])
t_a8, a8 = get_valid_plot_data(t_hr[m8], ar1[m8])
t_t8, tt8 = get_valid_plot_data(t_hr[m8], tau_s[m8])
ax.plot(t_v8, z(v8), label="Var (z)", color="#ff7f0e", lw=1.0)
ax.plot(t_a8, z(a8), label="AR(1) (z)", color="#2ca02c", lw=1.0)
ax.plot(t_t8, z(tt8), label="τ_s (z)", color="#9467bd", lw=1.2)
ax.axvline(event_hr, color="#d62728", lw=2, ls="--", label="event")
ax.set_xlabel("Time (h)")
ax.set_ylabel("z-scored (std)")
ax.set_title(f"Normalized EWS overlay (last 8h pre-event, Record {RECORD})")
ax.legend(loc="upper left", fontsize=8)
ax.grid(True, alpha=0.25)
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "13_ews_normalized_overlay.png"), dpi=150)
plt.close()

# 14. Heart rate (bpm) derived + event
bpm = 60000.0 / (rr + 1e-9)   # instantaneous
fig, ax = plt.subplots(figsize=(14, 3.8))
ax.plot(t_hr[::3], bpm[::3], lw=0.4, color="#17becf", alpha=0.7)
ax.axvline(event_hr, color="#d62728", lw=2.5, ls="--")
ax.axvspan(approach_start, event_hr, alpha=0.12, color="red")
ax.set_xlabel("Time (h)")
ax.set_ylabel("Heart rate (bpm)")
ax.set_title(f"Instantaneous heart rate (Record {RECORD})")
ax.grid(True, alpha=0.2)
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "14_hr_bpm.png"), dpi=150)
plt.close()

# 15. Full trace + markers for context (especially useful for intermediate event record 30)
fig, ax = plt.subplots(figsize=(14, 3.5))
ax.plot(t_hr[::2], rr[::2], lw=0.3, color="#1f77b4", alpha=0.6)
ax.axvline(event_hr, color="#d62728", lw=2, ls="--", label=f"event @ {event_hr:.2f}h")
ax.axvspan(approach_start, event_hr, alpha=0.12, color="red")
ax.set_xlabel("Time (h)")
ax.set_ylabel("RR (ms)")
ax.set_title(f"Full R-R trace with event marker (Record {RECORD})")
ax.legend(loc="upper right")
ax.grid(alpha=0.2)
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "01_rr_full.png"), dpi=140)
plt.close()

print("Additional diagnostic figures (12-15 / 01) generated.")

# === Save summary metrics JSON ===
import json
from datetime import datetime

summary = {
    "protocol": "CCTP v0.1",
    "record": RECORD,
    "source": "Sudden Cardiac Death Holter Database (PhysioNet sddb)",
    "vfon_sec": vfon_sec,
    "total_hours": total_h,
    "n_beats": int(len(rr)),
    "windows": {"W_EWS": W_EWS, "W_TAU": W_TAU, "stride": STRIDE},
    "basal_window_h": [float(basal_start), float(basal_end)],
    "approach_window_h": [float(approach_start), float(approach_end)],
    "metrics": {
        "var": {
            "basal_mean": float(np.nanmean(basal_var)),
            "approach_mean": float(np.nanmean(app_var)),
            "delta": float(np.nanmean(app_var) - np.nanmean(basal_var)),
            "p_welch": float(stats.ttest_ind(basal_var, app_var, equal_var=False, nan_policy="omit").pvalue)
        },
        "ar1": {
            "basal_mean": float(np.nanmean(basal_ar1)),
            "approach_mean": float(np.nanmean(app_ar1)),
            "delta": float(np.nanmean(app_ar1) - np.nanmean(basal_ar1)),
            "p_welch": float(stats.ttest_ind(basal_ar1, app_ar1, equal_var=False, nan_policy="omit").pvalue)
        },
        "tau_s": {
            "basal_mean": float(np.nanmean(basal_tau)),
            "approach_mean": float(np.nanmean(app_tau)),
            "delta": float(np.nanmean(app_tau) - np.nanmean(basal_tau)),
            "p_welch": float(stats.ttest_ind(basal_tau, app_tau, equal_var=False, nan_policy="omit").pvalue)
        }
    },
    "has_systemictau": HAS_SYSTEMICTAU,
    "timestamp": datetime.now().isoformat()
}

out_json = os.path.join(RES_DIR, f"cctp_pilot_summary_{RECORD}.json")
with open(out_json, "w") as f:
    json.dump(summary, f, indent=2)

print(f"\nSaved {out_json}")
print("\nAll figures and metrics generated for pilot.")
print("Key figures in:", FIG_DIR)
print("Next: inspect plots + decide whether to extend (TDA, RECD ordinal, more records).")
