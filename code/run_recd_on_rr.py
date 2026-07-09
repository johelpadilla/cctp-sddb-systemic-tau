#!/usr/bin/env python3
"""
run_recd_on_rr.py
Port of RECD ordinal levels (Φ₁, Φ₂, Φ₃ / excess3) to real RR series from CCTP pilot.

Surgical, falsifiable design (per spec):
- Input: exactly the same rr_{rec}_clean.npz (no re-cleaning).
- Same bivariate proxy as τ_s: X = [z(RR), z(|ΔRR|)]  (N=2 "variables" for conjunctions).
- Ordinal embedding: m=3, delay=1 (adjacent beats; delay=2 optional).
- Windows for Φ levels: W_PHI ~101 beats (matches W_TAU scale), stride=5.
- Metrics (reusing recd_ordinal_levels.py verbatim):
    Φ₁ (coincidence), Φ₂ (persistent relations), Φ₃ + excess3 (proxy of irreducible synergy, theta3=0.10)
    high_level3_rate = fraction of (excess3 > 1.75)
    mean_excess3 basal vs approach (primary for stats)
- Same basal/approach windows per record (reused logic from prior).
- Stats: Welch t-test on mean_excess3 (and on high-indicator series) between regimes.
- Between-record comparison (does Level 3 behave differently for terminal vs intermediate?).
- Outputs: per-record JSON + diagnostic plots (excess3(t) + context + boxplots).
- No lambda/RECD accumulation yet (pure levels test).

Usage (surgical):
  python3 code/run_recd_on_rr.py --record 35
  python3 code/run_recd_on_rr.py --record 30
  python3 code/run_recd_on_rr.py --record all

This tests directly whether the increase in nested ordinal Level 3 (excess3) explains / tracks
the relational reorganization previously captured by Δτ_s (context-dependent sign).
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

# === Import nested RECD levels (vendored in code/recd_ordinal_levels.py) ===
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
try:
    from _bootstrap import import_recd_ordinal_levels
    _recd = import_recd_ordinal_levels()
    generate_multivariate_symbols = _recd["generate_multivariate_symbols"]
    compute_phi1 = _recd["compute_phi1"]
    compute_phi2 = _recd["compute_phi2"]
    compute_phi3 = _recd["compute_phi3"]
    high_level3_rate = _recd["high_level3_rate"]
    HAS_RECD = True
except Exception as e:
    print(f"ERROR: could not import recd_ordinal_levels: {e}")
    HAS_RECD = False
    sys.exit(1)

# === Paths (match pilot/surrogates) ===
BASE = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
DATA_DIR = os.path.join(BASE, "data")
FIG_ROOT = os.path.join(BASE, "figures")
RES_DIR = os.path.join(BASE, "results")
os.makedirs(RES_DIR, exist_ok=True)

# === Parameters (surgical, per user design) ===
M = 3
DELAY = 1          # 1 = adjacent beats; try 2 for sensitivity if needed
W_PHI = 101        # ~ matches W_TAU scale for RR (beats)
STRIDE = 5
THETA3 = 0.10
HIGH_THRESH = 1.75  # for high_level3_rate

# Re-cal note (2026-07-08): for real RR use theta3=0.08, high=0.65 to make high_level3_rate informative. 
# These can be overridden via CLI.

# === Window logic (exact reuse) ===
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
    """Return non-NaN slice under mask (for strided/NaN-padded outputs)."""
    v = y[m]
    return v[~np.isnan(v)]

def build_proxy(rr):
    """Exact same minimal bivariate proxy used for τ_s."""
    drr = np.abs(np.diff(rr, prepend=rr[0]))
    rr_z = (rr - np.mean(rr)) / (np.std(rr) + 1e-12)
    drr_z = (drr - np.mean(drr)) / (np.std(drr) + 1e-12)
    return np.column_stack([rr_z, drr_z])

def compute_recd_levels_on_X(X, t_hr_for_X, m=M, delay=DELAY, w_phi=W_PHI, theta=THETA3, stride=STRIDE):
    """
    Compute Φ levels + excess3 on the proxy.
    Returns dict with arrays aligned to a symbol time base t_sym (length of S).
    """
    S = generate_multivariate_symbols(X, m=m, delay=delay)
    T_sym = S.shape[0]
    offset = (m - 1) * delay
    t_sym = t_hr_for_X[offset : offset + T_sym]

    phi1 = compute_phi1(S)
    phi2 = compute_phi2(S, d=4)

    # phi3/excess3 respect stride internally; still length T_sym with NaNs
    phi3, excess3 = compute_phi3(S, window=w_phi, theta=theta, stride=stride)

    return {
        "t_sym": t_sym,
        "S": S,
        "phi1": phi1,
        "phi2": phi2,
        "phi3": phi3,
        "excess3": excess3,
        "params": {"m": m, "delay": delay, "w_phi": w_phi, "theta3": theta, "stride": stride}
    }

def extract_regime_stats(excess3, t_sym, basal_start, basal_end, approach_start, approach_end, high_thresh=HIGH_THRESH):
    """Return means, deltas, rates, and Welch ps for the two regimes."""
    basal_mask = (t_sym >= basal_start) & (t_sym <= basal_end)
    app_mask = (t_sym >= approach_start) & (t_sym < approach_end)

    ex_b = get_valid(excess3, basal_mask)
    ex_a = get_valid(excess3, app_mask)

    if len(ex_b) < 5 or len(ex_a) < 5:
        return {"n_basal": len(ex_b), "n_app": len(ex_a), "valid": False}

    mean_b = float(np.nanmean(ex_b))
    mean_a = float(np.nanmean(ex_a))
    delta = mean_a - mean_b

    # Welch on the excess3 values themselves
    p_ex = float(stats.ttest_ind(ex_b, ex_a, equal_var=False, nan_policy="omit").pvalue)

    # high_level3_rate (scalar per regime) + t-test on the indicator series for a p-value
    ind_b = (ex_b > high_thresh).astype(float)
    ind_a = (ex_a > high_thresh).astype(float)
    rate_b = float(np.mean(ind_b))
    rate_a = float(np.mean(ind_a))
    delta_rate = rate_a - rate_b
    p_rate = float(stats.ttest_ind(ind_b, ind_a, equal_var=False).pvalue) if (np.std(ind_b) > 1e-12 or np.std(ind_a) > 1e-12) else None

    # Also basic phi means (for context)
    # Note: caller passes full arrays; we slice here with masks (dense for phi1/2)
    # We will compute outside for simplicity in main loop.

    return {
        "valid": True,
        "n_basal": len(ex_b),
        "n_app": len(ex_a),
        "mean_excess3": {"basal": mean_b, "approach": mean_a, "delta": delta, "p_welch": p_ex},
        "high_level3_rate": {"basal": rate_b, "approach": rate_a, "delta": delta_rate, "p_welch": p_rate},
        "high_thresh": high_thresh
    }

def plot_excess3(t_sym, excess3, event_hr, basal_start, basal_end, approach_start, approach_end,
                 rec, out_dir, rate_b, rate_a, high_thresh=HIGH_THRESH):
    """Main diagnostic: excess3(t) with regime shading, event, threshold."""
    fig, ax = plt.subplots(figsize=(13, 5))

    t_v, ex_v = t_sym, excess3  # use get_valid for clean line
    valid = ~np.isnan(ex_v)
    ax.plot(t_v[valid], ex_v[valid], lw=0.9, color="#8e44ad", label="excess3 (Level 3 proxy)")

    ax.axvline(event_hr, color="#d62728", lw=2.0, ls="--", label=f"event ({event_hr:.2f}h)")
    ax.axvspan(basal_start, basal_end, alpha=0.12, color="#1f77b4", label=f"basal [{basal_start:.1f}-{basal_end:.1f}h]")
    ax.axvspan(approach_start, approach_end, alpha=0.15, color="red", label=f"approach (last ~3h)")

    ax.axhline(high_thresh, color="#e74c3c", lw=1.5, ls=":", label=f"high3 thresh = {high_thresh}")

    ax.set_xlabel("Time (hours)")
    ax.set_ylabel("excess3 (combined synergy proxy)")
    ax.set_title(f"Record {rec} — RECD Level 3 (excess3) on RR proxy\n"
                 f"mean basal={rate_b:.3f} vs approach={rate_a:.3f} | high3_rate basal={rate_b:.3f} app={rate_a:.3f}")
    ax.legend(loc="upper right", fontsize=8)
    ax.grid(True, alpha=0.25)
    ax.set_ylim(bottom=max(-0.05, np.nanmin(ex_v[valid]) - 0.05))

    plt.tight_layout()
    fpath = os.path.join(out_dir, "16_recd_excess3.png")
    plt.savefig(fpath, dpi=150)
    plt.close()
    print(f"Saved {fpath}")
    return fpath

def plot_excess3_boxplots(ex_basal, ex_app, rec, basal_label, app_label, out_dir):
    """Boxplot comparison for excess3 (primary Level 3 metric)."""
    fig, ax = plt.subplots(figsize=(6, 4.5))
    ax.boxplot([ex_basal, ex_app], tick_labels=[basal_label, app_label],
               patch_artist=True,
               boxprops=dict(facecolor="#e8daef"), medianprops=dict(color="#8e44ad", lw=2))
    ax.set_ylabel("excess3")
    ax.set_title(f"excess3: Basal vs Approach (Record {rec})\nLevel 3 proxy")
    ax.grid(True, alpha=0.3, axis="y")
    plt.tight_layout()
    fpath = os.path.join(out_dir, "17_recd_excess3_box.png")
    plt.savefig(fpath, dpi=150)
    plt.close()
    print(f"Saved {fpath}")
    return fpath

def run_recd_for_record(record, theta3=None, high_thresh=None):
    theta = theta3 if theta3 is not None else THETA3
    high_t = high_thresh if high_thresh is not None else HIGH_THRESH
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

    print(f"\n=== Record {record} ===")
    print(f"event_hr={event_hr:.2f}h | basal=[{basal_start},{basal_end}] | approach=[{approach_start:.2f},{approach_end:.2f}]")
    print(f"RECD params: m={M}, delay={DELAY}, w_phi={W_PHI}, theta3={theta}, stride={STRIDE}  high_thresh={high_t}")

    # Proxy (identical to τ_s)
    X = build_proxy(rr)

    # RECD levels
    recd = compute_recd_levels_on_X(X, t_hr, m=M, delay=DELAY, w_phi=W_PHI, theta=theta, stride=STRIDE)
    t_sym = recd["t_sym"]
    excess3 = recd["excess3"]
    phi1 = recd["phi1"]
    phi2 = recd["phi2"]
    phi3 = recd["phi3"]

    # Regime extraction + stats
    stats_dict = extract_regime_stats(excess3, t_sym, basal_start, basal_end, approach_start, approach_end, high_thresh=high_t)

    # Slice for phis (dense) and report means
    basal_mask = (t_sym >= basal_start) & (t_sym <= basal_end)
    app_mask = (t_sym >= approach_start) & (t_sym < approach_end)

    phi1_b = get_valid(phi1, basal_mask)
    phi1_a = get_valid(phi1, app_mask)
    phi2_b = get_valid(phi2, basal_mask)
    phi2_a = get_valid(phi2, app_mask)
    phi3_b = get_valid(phi3, basal_mask)
    phi3_a = get_valid(phi3, app_mask)

    def safe_mean(a): return float(np.nanmean(a)) if len(a) > 0 else np.nan

    phi_means = {
        "phi1": {"basal": safe_mean(phi1_b), "approach": safe_mean(phi1_a)},
        "phi2": {"basal": safe_mean(phi2_b), "approach": safe_mean(phi2_a)},
        "phi3_active_frac": {"basal": safe_mean(phi3_b > 0), "approach": safe_mean(phi3_a > 0)},
    }

    pex = stats_dict['mean_excess3']['p_welch']
    prate = stats_dict['high_level3_rate']['p_welch']
    pex_str = f"{pex:.2e}" if pex is not None else "null"
    prate_str = f"{prate:.2e}" if prate is not None else "null"
    print(f"mean_excess3 basal={stats_dict['mean_excess3']['basal']:.5f}  approach={stats_dict['mean_excess3']['approach']:.5f}  Δ={stats_dict['mean_excess3']['delta']:+.5f}  p={pex_str}")
    print(f"high3_rate  basal={stats_dict['high_level3_rate']['basal']:.4f}  approach={stats_dict['high_level3_rate']['approach']:.4f}  Δ={stats_dict['high_level3_rate']['delta']:+.4f}  p={prate_str}")
    print(f"phi1 basal={phi_means['phi1']['basal']:.4f} app={phi_means['phi1']['approach']:.4f}")
    print(f"phi2 basal={phi_means['phi2']['basal']:.4f} app={phi_means['phi2']['approach']:.4f}")

    # Plots
    fig_dir = os.path.join(FIG_ROOT, record)
    os.makedirs(fig_dir, exist_ok=True)

    ex_b = get_valid(excess3, basal_mask)
    ex_a = get_valid(excess3, app_mask)

    fig16 = plot_excess3(t_sym, excess3, event_hr, basal_start, basal_end, approach_start, approach_end,
                         record, fig_dir,
                         stats_dict["high_level3_rate"]["basal"],
                         stats_dict["high_level3_rate"]["approach"],
                         high_thresh=high_t)

    fig17 = plot_excess3_boxplots(ex_b, ex_a, record,
                                  f"Basal\n({basal_start:.1f}-{basal_end:.1f}h)",
                                  "Approach\n(last ~3h)",
                                  fig_dir)

    # JSON
    out = {
        "protocol": "CCTP v0.1 RECD ordinal levels on RR",
        "record": record,
        "source": "SDDB Holter (same rr_*_clean.npz as τ_s pilot)",
        "params": {**recd["params"], "theta3": theta, "high_thresh": high_t},
        "windows": {
            "basal": [basal_start, basal_end],
            "approach": [approach_start, approach_end]
        },
        "n_sym": int(len(t_sym)),
        "mean_excess3": stats_dict["mean_excess3"],
        "high_level3_rate": stats_dict["high_level3_rate"],
        "phi_means": phi_means,
        "figures": [fig16, fig17],
        "timestamp": datetime.now().isoformat()
    }
    out_path = os.path.join(RES_DIR, f"recd_rr_{record}.json")
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2)
    print(f"Saved {out_path}")

    return {
        "record": record,
        "json": out_path,
        "mean_excess3_delta": stats_dict["mean_excess3"]["delta"],
        "p_ex": stats_dict["mean_excess3"]["p_welch"],
        "high_rate_delta": stats_dict["high_level3_rate"]["delta"],
        "p_rate": stats_dict["high_level3_rate"]["p_welch"],
    }

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--record", default="all", help="Record id: 30 | 35 | all")
    parser.add_argument("--theta3", type=float, default=None, help="Override theta3 (default 0.10; re-cal 0.08)")
    parser.add_argument("--high-thresh", type=float, default=None, help="Override high_thresh for high_level3_rate (default 1.75; re-cal 0.65)")
    args = parser.parse_args()

    theta = args.theta3 if args.theta3 is not None else THETA3
    high_t = args.high_thresh if args.high_thresh is not None else HIGH_THRESH

    if args.record == "all":
        records = ["30", "35"]
    else:
        records = [args.record]

    print("=" * 72)
    print("RECD Ordinal Levels (Φ1/Φ2/Φ3/excess3) on real RR series")
    print("Using exact same cleaned data + same basal/approach windows as τ_s pilot")
    print(f"Effective: theta3={theta} high_thresh={high_t}")
    print("=" * 72)

    results = []
    for rec in records:
        res = run_recd_for_record(rec, theta3=theta, high_thresh=high_t)
        results.append(res)

    # Comparative summary table
    if len(results) >= 2:
        print("\n" + "=" * 72)
        print("Comparative summary (Level 3 excess3 on RR proxy)")
        print(f"{'Record':<8} {'Δ mean_excess3':>16} {'p (Welch)':>12} {'Δ high3_rate':>14} {'p_rate':>10}")
        print("-" * 72)
        for r in results:
            print(f"{r['record']:<8} {r['mean_excess3_delta']:>+16.5f} {r['p_ex']:>12.2e} "
                  f"{r['high_rate_delta']:>+14.4f} {r['p_rate']:>10.2e}")
        print("=" * 72)
        print("Interpretation key: positive Δ = higher Level 3 synergy/irreducibility in approach (pre-event).")
        print("Compare sign/direction vs Δτ_s from prior (terminal 35 vs intermediate 30).")

    print("\nDone. JSONs in results/ ; figures in figures/{30,35}/ (16_*.png, 17_*.png)")
    print("Next: inspect plots, update HANDOFF + report, decide follow-ups (e.g. record 31, full RECD with α(λ)).")

if __name__ == "__main__":
    main()
