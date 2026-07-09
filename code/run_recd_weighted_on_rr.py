#!/usr/bin/env python3
"""
run_recd_weighted_on_rr.py
Full weighted RECD (RECD ponderado) on real RR series using λ(t) derived from |τ_s|.

Surgical design (per spec):
- Exact same rr_{rec}_clean.npz inputs, same proxy X=[z(RR), z(|ΔRR|)], same windows.
- m=3, delay=1 ; w_phi=101 (for phi3), stride=5 for consistency.
- Compute τ_s series (W_TAU=101, stride=1 for smooth λ) → λ(t) via compute_lambda.
- alpha_mode="lambda": pass lam_override derived from empirical |τ_s|.
- Call compute_recd_from_conjunctions with beta1=2.0, gamma2=1.5, gamma3=6.0 (stronger Nivel 3), delta3=2.0.
- Metrics per regime (basal vs approach):
    - mean_excess3, high_level3_rate (ref)
    - contrib1, contrib2, contrib3 (α_k · Φ_k)
    - frac_contrib3 = mean(contrib3) / mean(total)
    - mean_delta_recd
- Welch tests on mean_excess3, contrib3, and pointwise frac_contrib3 series.
- Outputs: recd_weighted_rr_*.json + diagnostic plots (excess3(t), contrib stacked, frac_contrib3(t), boxes).
- Reuses recd_ordinal_levels.py + systemictau.core.compute_taus exactly (no reimplementation).

Usage:
  python3 code/run_recd_weighted_on_rr.py --record 35
  python3 code/run_recd_weighted_on_rr.py --record 30
  python3 code/run_recd_weighted_on_rr.py --record all

This tests whether the weighted contribution of Nivel 3 (α3·Φ3) increases in approach
in a manner consistent with the synthetic post-Feigenbaum behavior, using λ(t) from observed τ_s.
"""

import os
import sys
import argparse
import json
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from datetime import datetime
import warnings
warnings.filterwarnings("ignore")

# === Imports: vendored RECD levels + systemictau for λ(t) from τ_s ===
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
try:
    from _bootstrap import import_recd_ordinal_levels, import_systemictau_core
    _recd = import_recd_ordinal_levels()
    generate_multivariate_symbols = _recd["generate_multivariate_symbols"]
    compute_phi1 = _recd["compute_phi1"]
    compute_phi2 = _recd["compute_phi2"]
    compute_phi3 = _recd["compute_phi3"]
    compute_recd_from_conjunctions = _recd["compute_recd_from_conjunctions"]
    compute_weighted_contributions = _recd["compute_weighted_contributions"]
    compute_lambda = _recd["compute_lambda"]
    high_level3_rate = _recd["high_level3_rate"]
    HAS_RECD = True
except Exception as e:
    print(f"ERROR: could not import recd_ordinal_levels: {e}")
    HAS_RECD = False
    sys.exit(1)

compute_taus, _systemic_tau, HAS_SYSTEMICTAU = import_systemictau_core()
if not HAS_SYSTEMICTAU:
    print("ERROR: could not import systemictau.core. Required for λ(t) from τ_s.")
    print("Install with: pip install systemictau  (or set SYSTEMICTAU_SRC)")
    sys.exit(1)

# === Paths (match pilot) ===
BASE = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
DATA_DIR = os.path.join(BASE, "data")
FIG_ROOT = os.path.join(BASE, "figures")
RES_DIR = os.path.join(BASE, "results")
os.makedirs(RES_DIR, exist_ok=True)

# === Parameters (surgical + consistent with prior RECD + τ_s pilot) ===
M = 3
DELAY = 1
W_PHI = 101          # window for phi3 / excess3 (matches prior RECD port)
STRIDE = 5           # for reporting density (phi3 stride inside)
W_TAU = 101          # for τ_s (same as pilot)
THETA3 = 0.10
HIGH_THRESH = 1.75
LAMBDA_THETA = 0.41  # default chaos thresh for compute_lambda; re-cal ~0.08-0.10 or use relative

# Alpha defaults per user spec for this surgical step (stronger Nivel 3 ramp)
BETA1 = 2.0
GAMMA2 = 1.5
GAMMA3 = 6.0   # elevated vs module default (3.0) to emphasize Level 3
DELTA3 = 2.0

# Re-cal note (2026-07-08): theta3=0.08, high=0.65, lambda_theta~0.09 or relative for real RR physiology.


def get_event_and_windows(record, t_hr, vfon_hr):
    total_h = float(t_hr.max())
    if record == "35" or (vfon_hr > total_h - 0.8):
        event_hr = total_h
    else:
        event_hr = vfon_hr
    approach_start = max(0.0, event_hr - 3.0)
    approach_end = event_hr

    if record == "35":
        basal_start, basal_end = 6.0, 16.0
    elif record == "30":
        basal_start, basal_end = 0.5, 3.5
    else:
        b_end = max(3.5, approach_start - 4.0)
        basal_start = max(0.5, b_end - 3.0)
        basal_end = b_end
    return event_hr, (basal_start, basal_end), (approach_start, approach_end)


def get_valid(y, m):
    """Non-NaN values under boolean mask."""
    v = y[m]
    return v[~np.isnan(v)]


def build_proxy(rr):
    """Exact same minimal bivariate proxy used for τ_s and prior RECD levels."""
    drr = np.abs(np.diff(rr, prepend=rr[0]))
    rr_z = (rr - np.mean(rr)) / (np.std(rr) + 1e-12)
    drr_z = (drr - np.mean(drr)) / (np.std(drr) + 1e-12)
    return np.column_stack([rr_z, drr_z])


def compute_tau_series(X, W_TAU=W_TAU, stride=1):
    """High-resolution τ_s (stride=1) for smooth λ(t). Stride=5 version available for plots if needed."""
    taus, _ = compute_taus(X, window_size=W_TAU, stride=stride)
    return taus


def compute_recd_weighted_on_X(X, t_hr_for_X, lam_override=None, theta3=None):
    """
    Run full pipeline:
      S → Φ1/Φ2/Φ3/excess3 → λ (from |τ_s| or override) → α(λ) → delta_recd, T_recd, contribs.
    Returns the dict from compute_recd_from_conjunctions + contrib arrays.
    """
    t3 = theta3 if theta3 is not None else THETA3
    recd = compute_recd_from_conjunctions(
        X,
        lam_override=lam_override,
        m=M,
        d=4,
        theta3=t3,
        window_tau=W_PHI,
        beta1=BETA1,
        gamma2=GAMMA2,
        gamma3=GAMMA3,
        delta3=DELTA3,
    )

    # Derived weighted contributions (pointwise)
    contribs = compute_weighted_contributions(recd)
    c1 = recd["alpha1"] * np.nan_to_num(recd["phi1"], nan=0.0)
    c2 = recd["alpha2"] * np.nan_to_num(recd["phi2"], nan=0.0)
    c3 = recd["alpha3"] * np.nan_to_num(recd["phi3"], nan=0.0)
    total_c = c1 + c2 + c3
    frac_c3 = np.divide(c3, total_c + 1e-12, where=(total_c > 1e-12))

    recd["contrib1"] = c1
    recd["contrib2"] = c2
    recd["contrib3"] = c3
    recd["frac_contrib3"] = frac_c3
    recd["total_contrib"] = total_c

    return recd


def extract_weighted_regime_stats(recd, t_recd, basal_start, basal_end, approach_start, approach_end, high_thresh=HIGH_THRESH):
    """Per-regime means + Welch on excess3, contrib3, frac_contrib3, delta_recd."""
    basal_mask = (t_recd >= basal_start) & (t_recd <= basal_end)
    app_mask = (t_recd >= approach_start) & (t_recd < approach_end)

    def gv(arr, mask):
        return get_valid(arr, mask)

    ex_b = gv(recd["excess3"], basal_mask)
    ex_a = gv(recd["excess3"], app_mask)

    if len(ex_b) < 5 or len(ex_a) < 5:
        return {"valid": False, "n_basal": len(ex_b), "n_app": len(ex_a)}

    # excess3
    mean_ex_b = float(np.nanmean(ex_b))
    mean_ex_a = float(np.nanmean(ex_a))
    delta_ex = mean_ex_a - mean_ex_b
    p_ex = float(stats.ttest_ind(ex_b, ex_a, equal_var=False, nan_policy="omit").pvalue)

    # high_level3_rate (reference)
    rate_b = float(np.nanmean(ex_b > high_thresh))
    rate_a = float(np.nanmean(ex_a > high_thresh))
    p_rate = None
    if (np.std(ex_b > high_thresh) > 1e-12) or (np.std(ex_a > high_thresh) > 1e-12):
        p_rate = float(stats.ttest_ind((ex_b > high_thresh).astype(float),
                                       (ex_a > high_thresh).astype(float),
                                       equal_var=False).pvalue)

    # contribs
    c3_b = gv(recd["contrib3"], basal_mask)
    c3_a = gv(recd["contrib3"], app_mask)
    c1_b = gv(recd["contrib1"], basal_mask)
    c1_a = gv(recd["contrib1"], app_mask)
    c2_b = gv(recd["contrib2"], basal_mask)
    c2_a = gv(recd["contrib2"], app_mask)

    mc1_b, mc1_a = float(np.nanmean(c1_b)), float(np.nanmean(c1_a))
    mc2_b, mc2_a = float(np.nanmean(c2_b)), float(np.nanmean(c2_a))
    mc3_b, mc3_a = float(np.nanmean(c3_b)), float(np.nanmean(c3_a))
    delta_c3 = mc3_a - mc3_b
    p_c3 = None
    if (np.std(c3_b) > 1e-12 or np.std(c3_a) > 1e-12):
        p_c3 = float(stats.ttest_ind(c3_b, c3_a, equal_var=False, nan_policy="omit").pvalue)

    # frac_contrib3 pointwise then mean + Welch on the frac series
    f3_b = gv(recd["frac_contrib3"], basal_mask)
    f3_a = gv(recd["frac_contrib3"], app_mask)
    mf3_b = float(np.nanmean(f3_b))
    mf3_a = float(np.nanmean(f3_a))
    delta_f3 = mf3_a - mf3_b
    p_f3 = float(stats.ttest_ind(f3_b, f3_a, equal_var=False, nan_policy="omit").pvalue) if (np.std(f3_b) > 1e-12 or np.std(f3_a) > 1e-12) else None

    # delta_recd
    dr_b = gv(recd["delta_recd"], basal_mask)
    dr_a = gv(recd["delta_recd"], app_mask)
    mdr_b = float(np.nanmean(dr_b))
    mdr_a = float(np.nanmean(dr_a))
    delta_dr = mdr_a - mdr_b
    p_dr = None
    if (np.std(dr_b) > 1e-12 or np.std(dr_a) > 1e-12):
        p_dr = float(stats.ttest_ind(dr_b, dr_a, equal_var=False, nan_policy="omit").pvalue)

    # lambda (diagnostic)
    lam_b = gv(recd["lambda"], basal_mask)
    lam_a = gv(recd["lambda"], app_mask)
    mlam_b = float(np.nanmean(lam_b))
    mlam_a = float(np.nanmean(lam_a))

    return {
        "valid": True,
        "n_basal": len(ex_b),
        "n_app": len(ex_a),
        "mean_excess3": {"basal": mean_ex_b, "approach": mean_ex_a, "delta": delta_ex, "p_welch": p_ex},
        "high_level3_rate": {"basal": rate_b, "approach": rate_a, "delta": rate_a - rate_b, "p_welch": p_rate},
        "contrib": {
            "contrib1": {"basal": mc1_b, "approach": mc1_a, "delta": mc1_a - mc1_b, "p_welch": None},
            "contrib2": {"basal": mc2_b, "approach": mc2_a, "delta": mc2_a - mc2_b, "p_welch": None},
            "contrib3": {"basal": mc3_b, "approach": mc3_a, "delta": delta_c3, "p_welch": p_c3},
        },
        "frac_contrib3": {"basal": mf3_b, "approach": mf3_a, "delta": delta_f3, "p_welch": p_f3},
        "mean_delta_recd": {"basal": mdr_b, "approach": mdr_a, "delta": delta_dr, "p_welch": p_dr},
        "mean_lambda": {"basal": mlam_b, "approach": mlam_a},
        "high_thresh": high_thresh,
    }


def plot_weighted_time_series(t_recd, recd, event_hr, basal_start, basal_end, approach_start, approach_end,
                              rec, out_dir, stats_dict):
    """Rich diagnostic: excess3 + stacked/layered contribs + frac_contrib3 with regime shading."""
    fig, axes = plt.subplots(3, 1, figsize=(13, 9), sharex=True)

    # 1. excess3
    ax = axes[0]
    valid = ~np.isnan(recd["excess3"])
    ax.plot(t_recd[valid], recd["excess3"][valid], lw=0.85, color="#8e44ad", label="excess3 (Nivel 3 proxy)")
    ax.axhline(high_t, color="#e74c3c", lw=1.2, ls=":", label=f"high thresh={high_t} (ref)")
    ax.axvline(event_hr, color="#d62728", lw=1.8, ls="--", label=f"VF event ({event_hr:.2f}h)")
    ax.axvspan(basal_start, basal_end, alpha=0.10, color="#1f77b4", label=f"basal [{basal_start:.1f}-{basal_end:.1f}h]")
    ax.axvspan(approach_start, approach_end, alpha=0.13, color="red", label="approach (last ~3h)")
    ax.set_ylabel("excess3")
    ax.set_title(f"Record {rec} — Weighted RECD (λ from |τ_s| , α3 ramp γ3={GAMMA3})\n"
                 f"mean_excess3 basal={stats_dict['mean_excess3']['basal']:.4f} → approach={stats_dict['mean_excess3']['approach']:.4f} (Δ={stats_dict['mean_excess3']['delta']:+.5f})")
    ax.legend(loc="upper right", fontsize=7)
    ax.grid(True, alpha=0.25)
    ax.set_ylim(bottom=max(-0.05, np.nanmin(recd["excess3"][valid]) - 0.05))

    # 2. contribs (overlaid or stacked)
    ax = axes[1]
    c1v = recd["contrib1"]
    c2v = recd["contrib2"]
    c3v = recd["contrib3"]
    valid_c = ~np.isnan(c3v)
    ax.plot(t_recd[valid_c], c1v[valid_c], lw=0.7, color="#3498db", label="contrib1 (α1·Φ1)")
    ax.plot(t_recd[valid_c], c2v[valid_c], lw=0.7, color="#2ecc71", label="contrib2 (α2·Φ2)")
    ax.plot(t_recd[valid_c], c3v[valid_c], lw=0.9, color="#e74c3c", label="contrib3 (α3·Φ3) ← key")
    ax.axvline(event_hr, color="#d62728", lw=1.5, ls="--")
    ax.axvspan(basal_start, basal_end, alpha=0.10, color="#1f77b4")
    ax.axvspan(approach_start, approach_end, alpha=0.13, color="red")
    ax.set_ylabel("α_k · Φ_k  (contrib)")
    ax.legend(loc="upper right", fontsize=7)
    ax.grid(True, alpha=0.25)

    # 3. frac_contrib3
    ax = axes[2]
    fv = recd["frac_contrib3"]
    valid_f = ~np.isnan(fv)
    ax.plot(t_recd[valid_f], fv[valid_f], lw=0.9, color="#9b59b6", label="frac_contrib3 = c3 / (c1+c2+c3)")
    ax.axvline(event_hr, color="#d62728", lw=1.5, ls="--")
    ax.axvspan(basal_start, basal_end, alpha=0.10, color="#1f77b4")
    ax.axvspan(approach_start, approach_end, alpha=0.13, color="red")
    ax.axhline(1.0/3, color="#7f8c8d", lw=0.8, ls=":", label="1/3 uniform ref")
    ax.set_ylabel("frac contrib3")
    ax.set_xlabel("Time (hours)")
    ax.legend(loc="upper right", fontsize=7)
    ax.grid(True, alpha=0.25)
    ax.set_ylim(0, max(0.6, np.nanmax(fv[valid_f]) * 1.05))

    plt.tight_layout()
    fpath = os.path.join(out_dir, "18_recd_weighted_contribs.png")
    plt.savefig(fpath, dpi=150)
    plt.close()
    print(f"Saved {fpath}")
    return fpath


def plot_weighted_boxplots(ex_b, ex_a, c3_b, c3_a, f3_b, f3_a, rec, basal_label, app_label, out_dir):
    """Comparative boxplots for key weighted metrics."""
    fig, axes = plt.subplots(1, 3, figsize=(12, 4.2))

    # excess3
    ax = axes[0]
    ax.boxplot([ex_b, ex_a], tick_labels=[basal_label, app_label], patch_artist=True,
               boxprops=dict(facecolor="#e8daef"), medianprops=dict(color="#8e44ad", lw=2))
    ax.set_ylabel("excess3")
    ax.set_title("excess3 (Nivel 3 proxy)")
    ax.grid(True, alpha=0.3, axis="y")

    # contrib3
    ax = axes[1]
    ax.boxplot([c3_b, c3_a], tick_labels=[basal_label, app_label], patch_artist=True,
               boxprops=dict(facecolor="#fadbd8"), medianprops=dict(color="#c0392b", lw=2))
    ax.set_ylabel("contrib3 (α3·Φ3)")
    ax.set_title("Weighted Nivel 3 contrib")
    ax.grid(True, alpha=0.3, axis="y")

    # frac
    ax = axes[2]
    ax.boxplot([f3_b, f3_a], tick_labels=[basal_label, app_label], patch_artist=True,
               boxprops=dict(facecolor="#d5f5e3"), medianprops=dict(color="#27ae60", lw=2))
    ax.set_ylabel("frac_contrib3")
    ax.set_title("Fraction of total ΔRECD from Nivel 3")
    ax.grid(True, alpha=0.3, axis="y")

    fig.suptitle(f"Record {rec} — Weighted RECD: Basal vs Approach", fontsize=11)
    plt.tight_layout()
    fpath = os.path.join(out_dir, "19_recd_weighted_box.png")
    plt.savefig(fpath, dpi=150)
    plt.close()
    print(f"Saved {fpath}")
    return fpath


def run_weighted_for_record(record, theta3=None, high_thresh=None, lambda_theta=None, lambda_relative=False):
    theta = theta3 if theta3 is not None else THETA3
    high_t = high_thresh if high_thresh is not None else HIGH_THRESH
    lam_t = lambda_theta if lambda_theta is not None else LAMBDA_THETA
    use_rel = lambda_relative

    rr_npz = os.path.join(DATA_DIR, f"rr_{record}_clean.npz")
    if not os.path.exists(rr_npz):
        raise FileNotFoundError(f"Missing {rr_npz}")

    d = np.load(rr_npz)
    rr = d["rr_ms"].astype(float)
    t_sec = d["t_sec"].astype(float)
    vfon_sec = float(d["vfon_sec"])

    t_hr = t_sec / 3600.0
    vfon_hr = vfon_sec / 3600.0

    event_hr, (basal_start, basal_end), (approach_start, approach_end) = get_event_and_windows(record, t_hr, vfon_hr)

    print(f"\n=== Record {record} (Weighted RECD) ===")
    print(f"event_hr={event_hr:.2f}h | basal=[{basal_start},{basal_end}] | approach=[{approach_start:.2f},{approach_end:.2f}]")
    print(f"Params: m={M}, delay={DELAY}, W_TAU={W_TAU}, w_phi={W_PHI}, theta3={theta}, stride={STRIDE}  high={high_t}")
    print(f"alpha: beta1={BETA1}, gamma2={GAMMA2}, gamma3={GAMMA3}, delta3={DELTA3}  (alpha_mode=lambda)  lam_theta={lam_t} rel={use_rel}")

    X = build_proxy(rr)

    # τ_s high-res for λ(t)
    taus_hr = compute_tau_series(X, W_TAU=W_TAU, stride=1)
    if use_rel:
        max_abs = np.nanmax(np.abs(taus_hr)) or 1.0
        lam_hr = np.abs(taus_hr) / max_abs
    else:
        lam_hr = compute_lambda(taus_hr, theta_chaos=lam_t)
    # NaN-safe
    lam_hr = np.nan_to_num(lam_hr, nan=0.0)

    # Align lam to symbol time base (offset for m=3, delay=1)
    offset = (M - 1) * DELAY
    T_full = X.shape[0]
    T_eff = T_full - offset
    lam_for_recd = lam_hr[offset : offset + T_eff] if len(lam_hr) > offset else np.zeros(T_eff)

    # Full weighted RECD
    recd = compute_recd_weighted_on_X(X, t_hr, lam_override=lam_for_recd, theta3=theta)

    t_recd = t_hr[offset : offset + len(recd["phi1"])]

    # Regime stats
    stats_dict = extract_weighted_regime_stats(
        recd, t_recd, basal_start, basal_end, approach_start, approach_end, high_thresh=high_t
    )

    pex = stats_dict['mean_excess3']['p_welch']
    pc3 = stats_dict['contrib']['contrib3']['p_welch']
    pf3 = stats_dict['frac_contrib3']['p_welch']
    pex_str = f"{pex:.2e}" if pex is not None else "null"
    pc3_str = f"{pc3:.2e}" if pc3 is not None else "null"
    pf3_str = f"{pf3:.2e}" if pf3 is not None else "null"
    print(f"mean_excess3  basal={stats_dict['mean_excess3']['basal']:.5f}  app={stats_dict['mean_excess3']['approach']:.5f}  Δ={stats_dict['mean_excess3']['delta']:+.5f}  p={pex_str}")
    print(f"contrib3      basal={stats_dict['contrib']['contrib3']['basal']:.6f}  app={stats_dict['contrib']['contrib3']['approach']:.6f}  Δ={stats_dict['contrib']['contrib3']['delta']:+.6f}  p={pc3_str}")
    print(f"frac_contrib3 basal={stats_dict['frac_contrib3']['basal']:.5f}  app={stats_dict['frac_contrib3']['approach']:.5f}  Δ={stats_dict['frac_contrib3']['delta']:+.5f}  p={pf3_str}")
    print(f"mean_lambda   basal={stats_dict['mean_lambda']['basal']:.5f}  app={stats_dict['mean_lambda']['approach']:.5f}")

    # Plots
    fig_dir = os.path.join(FIG_ROOT, record)
    os.makedirs(fig_dir, exist_ok=True)

    ex_b = get_valid(recd["excess3"], (t_recd >= basal_start) & (t_recd <= basal_end))
    ex_a = get_valid(recd["excess3"], (t_recd >= approach_start) & (t_recd < approach_end))
    c3_b = get_valid(recd["contrib3"], (t_recd >= basal_start) & (t_recd <= basal_end))
    c3_a = get_valid(recd["contrib3"], (t_recd >= approach_start) & (t_recd < approach_end))
    f3_b = get_valid(recd["frac_contrib3"], (t_recd >= basal_start) & (t_recd <= basal_end))
    f3_a = get_valid(recd["frac_contrib3"], (t_recd >= approach_start) & (t_recd < approach_end))

    fig18 = plot_weighted_time_series(
        t_recd, recd, event_hr, basal_start, basal_end, approach_start, approach_end,
        record, fig_dir, stats_dict
    )

    fig19 = plot_weighted_boxplots(
        ex_b, ex_a, c3_b, c3_a, f3_b, f3_a,
        record,
        f"Basal\n({basal_start:.1f}-{basal_end:.1f}h)",
        "Approach\n(last ~3h)",
        fig_dir
    )

    # JSON
    out = {
        "protocol": "CCTP v0.1 RECD weighted (alpha_mode=lambda from |τ_s|)",
        "record": record,
        "source": "SDDB Holter (same rr_*_clean.npz)",
        "params": {
            "m": M, "delay": DELAY,
            "W_TAU": W_TAU, "w_phi": W_PHI, "stride": STRIDE,
            "theta3": theta, "high_thresh": high_t, "lambda_theta": lam_t, "lambda_relative": use_rel,
            "alpha_mode": "lambda",
            "beta1": BETA1, "gamma2": GAMMA2, "gamma3": GAMMA3, "delta3": DELTA3
        },
        "windows": {
            "basal": [basal_start, basal_end],
            "approach": [approach_start, approach_end]
        },
        "n_recd": int(len(t_recd)),
        "mean_excess3": stats_dict["mean_excess3"],
        "high_level3_rate": stats_dict["high_level3_rate"],
        "contrib": stats_dict["contrib"],
        "frac_contrib3": stats_dict["frac_contrib3"],
        "mean_delta_recd": stats_dict["mean_delta_recd"],
        "mean_lambda": stats_dict["mean_lambda"],
        "figures": [fig18, fig19],
        "timestamp": datetime.now().isoformat()
    }
    out_path = os.path.join(RES_DIR, f"recd_weighted_rr_{record}.json")
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2)
    print(f"Saved {out_path}")

    return {
        "record": record,
        "json": out_path,
        "delta_ex": stats_dict["mean_excess3"]["delta"],
        "p_ex": stats_dict["mean_excess3"]["p_welch"],
        "delta_c3": stats_dict["contrib"]["contrib3"]["delta"],
        "p_c3": stats_dict["contrib"]["contrib3"]["p_welch"],
        "delta_f3": stats_dict["frac_contrib3"]["delta"],
        "p_f3": stats_dict["frac_contrib3"]["p_welch"],
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--record", default="all", help="Record id: 30 | 35 | all")
    parser.add_argument("--theta3", type=float, default=None, help="Override theta3 (re-cal 0.08)")
    parser.add_argument("--high-thresh", type=float, default=None, help="Override high_thresh (re-cal 0.65)")
    parser.add_argument("--lambda-theta", type=float, default=None, help="Override theta_chaos for lambda (re-cal ~0.09); or use --lambda-relative")
    parser.add_argument("--lambda-relative", action="store_true", help="Compute lambda as |tau_s| / max_|tau_s| (record-relative; good for small observed deltas ~0.026-0.088)")
    args = parser.parse_args()

    theta = args.theta3 if args.theta3 is not None else THETA3
    high_t = args.high_thresh if args.high_thresh is not None else HIGH_THRESH
    lam_theta = args.lambda_theta if args.lambda_theta is not None else LAMBDA_THETA
    use_rel = args.lambda_relative

    if args.record == "all":
        records = ["30", "35"]
    else:
        records = [args.record]

    print("=" * 72)
    print("RECD Weighted (Φ + α(λ) from |τ_s|) on real RR series")
    print("alpha_mode=lambda | same proxy/windows as prior τ_s + RECD levels steps")
    print(f"Effective: theta3={theta} high_thresh={high_t} lambda_theta={lam_theta} relative={use_rel}")
    print("=" * 72)

    results = []
    for rec in records:
        res = run_weighted_for_record(rec, theta3=theta, high_thresh=high_t, lambda_theta=lam_theta, lambda_relative=use_rel)
        results.append(res)

    if len(results) >= 2:
        print("\n" + "=" * 72)
        print("Comparative summary (Weighted RECD on RR proxy)")
        print(f"{'Record':<8} {'Δ excess3':>12} {'p_ex':>10} {'Δ contrib3':>12} {'p_c3':>10} {'Δ frac3':>10} {'p_f3':>10}")
        print("-" * 72)
        for r in results:
            pex = r['p_ex']
            pc3 = r['p_c3']
            pf3 = r['p_f3']
            print(f"{r['record']:<8} {r['delta_ex']:>+12.5f} {pex if pex is None else f'{pex:>10.2e}'} "
                  f"{r['delta_c3']:>+12.6f} {pc3 if pc3 is None else f'{pc3:>10.2e}'} "
                  f"{r['delta_f3']:>+10.5f} {pf3 if pf3 is None else f'{pf3:>10.2e}'}")
        print("=" * 72)
        print("Key test: Does frac_contrib3 or contrib3 rise in approach (esp. record 35 terminal)?")
        print("Direction of Δ should be read against prior Δτ_s and Δexcess3 (context-dependent).")

    print("\nDone. JSONs: results/recd_weighted_rr_*.json")
    print("Figures: figures/{30,35}/18_recd_weighted_contribs.png + 19_recd_weighted_box.png")
    print("Next: inspect whether Nivel 3 weighted weight increases coherently with τ_s sign; update docs.")


if __name__ == "__main__":
    main()
