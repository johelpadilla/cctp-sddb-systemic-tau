#!/usr/bin/env python3
"""
run_cctp_surrogates.py
Light per-record surrogates (phase-shuffle independent) on Δτ_s for CCTP v0.1.

Design (per user spec):
- Phase-shuffle independiente por componente del proxy bivariado (RR + |ΔRR|).
- 8 surrogates por registro (granularidad 0.125).
- Pipeline idéntico a analyze_cctp_pilot.py:
    mismos W_EWS/W_TAU/STRIDE, mismo proxy z(RR)+z(|dRR|), mismo compute_taus,
    mismas ventanas basal/approach por RECORD, mismas máscaras y extracción de medias.
- Métrica: Δτ_s = mean(approach) - mean(basal)
- Test: p = proporción de surrogates donde |Δ_surr| >= |Δ_obs| (two-sided)
         + p direccional específico según signo del observado.
- Registros: 30 (prioridad, Δ negativo) + 35 (Δ positivo, ya teníamos un run previo).
- Salidas: JSON por registro + histograma simple por registro + tabla comparativa clara.

Uso quirúrgico:
  python3 code/run_cctp_surrogates.py --record 30
  python3 code/run_cctp_surrogates.py --record 35
  python3 code/run_cctp_surrogates.py --record all
"""

import os
import sys
import argparse
import json
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
import warnings
warnings.filterwarnings("ignore")

# === Import systemictau (pip package or SYSTEMICTAU_SRC / monorepo) ===
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _bootstrap import import_systemictau_core
compute_taus, _systemic_tau, HAS_SYSTEMICTAU = import_systemictau_core()
if not HAS_SYSTEMICTAU:
    print("ERROR: could not import systemictau.core. Cannot run surrogates.")
    print("Install with: pip install systemictau  (or set SYSTEMICTAU_SRC)")
    sys.exit(1)

# === Paths ===
BASE = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
DATA_DIR = os.path.join(BASE, "data")
FIG_ROOT = os.path.join(BASE, "figures")
RES_DIR = os.path.join(BASE, "results")
os.makedirs(RES_DIR, exist_ok=True)

# === Surrogate helpers (phase-shuffle per component, IAAFT style) ===
# Copied/adapted from established implementation used in prior work for H0 "no relational structure".
def iaaft_surrogate_1d(x: np.ndarray, n_iters: int = 30, seed: int = None) -> np.ndarray:
    """IAAFT 1D surrogate: preserves amplitude distribution + approx power spectrum."""
    rng = np.random.default_rng(seed)
    x = np.asarray(x).ravel()
    n = len(x)
    xf = np.fft.rfft(x)
    amp = np.abs(xf)
    phase = rng.uniform(0, 2 * np.pi, len(xf))
    phase[0] = 0.0
    if len(xf) > 1:
        phase[-1] = 0.0
    s = np.fft.irfft(amp * np.exp(1j * phase), n=n)
    x_sorted = np.sort(x)
    for _ in range(n_iters):
        sf = np.fft.rfft(s)
        s = np.fft.irfft(amp * np.exp(1j * np.angle(sf)), n=n)
        ranks = np.argsort(np.argsort(s))
        s = x_sorted[ranks]
    return s

def phase_shuffle_independent(X: np.ndarray, seed: int = None) -> np.ndarray:
    """
    Independent phase-shuffle per column.
    Destroys cross-relations and temporal ordinal structure between the two
    components while approximately preserving marginal spectra + distributions.
    This is the null for "Systemic Tau rise is just linear auto-stats".
    """
    rng = np.random.default_rng(seed)
    X = np.asarray(X).copy()
    T, N = X.shape
    for i in range(N):
        X[:, i] = iaaft_surrogate_1d(X[:, i], n_iters=30, seed=rng.integers(0, 2**32))
    return X

# === Core metric extraction (mirrors analyze_cctp_pilot.py exactly) ===
W_TAU = 101
STRIDE = 5

def get_event_and_windows(record, t_hr, vfon_hr):
    total_h = float(t_hr.max())
    if record == "35" or (vfon_hr > total_h - 0.8):
        event_hr = total_h
    else:
        event_hr = vfon_hr
    approach_start = max(0.0, event_hr - 3.0)
    approach_end = event_hr

    # Robust basal: preserve exact prior numbers for 30/35; sensible pre-approach window otherwise
    if record == "35":
        basal_start, basal_end = 6.0, 16.0
    elif record == "30":
        basal_start, basal_end = 0.5, 3.5
    else:
        b_end = max(3.5, approach_start - 4.0)
        basal_start = max(0.5, b_end - 3.0)
        basal_end = b_end
    return event_hr, (basal_start, basal_end), (approach_start, approach_end)

def build_proxy_and_tau_s(rr, W_TAU=W_TAU, STRIDE=STRIDE):
    drr = np.abs(np.diff(rr, prepend=rr[0]))
    rr_z = (rr - np.mean(rr)) / (np.std(rr) + 1e-12)
    drr_z = (drr - np.mean(drr)) / (np.std(drr) + 1e-12)
    X = np.column_stack([rr_z, drr_z])
    taus_global, _ = compute_taus(X, window_size=W_TAU, stride=STRIDE)
    return X, taus_global

def build_tau_s_from_X(X, W_TAU=W_TAU, STRIDE=STRIDE):
    """Compute tau_s series on a (possibly surrogate) proxy matrix. Time axis unchanged."""
    taus_global, _ = compute_taus(X, window_size=W_TAU, stride=STRIDE)
    return taus_global

def extract_delta_tau_s(taus, t_hr, basal_start, basal_end, approach_start, approach_end):
    """Extract mean tau in windows using identical logic to the pilot (NaN filtering on strided output)."""
    basal_mask = (t_hr >= basal_start) & (t_hr <= basal_end)
    app_mask = (t_hr >= approach_start) & (t_hr < approach_end)

    def get_valid(y, m):
        v = y[m]
        return v[~np.isnan(v)]

    basal_tau = get_valid(taus, basal_mask)
    app_tau = get_valid(taus, app_mask)

    if len(basal_tau) < 5 or len(app_tau) < 5:
        return np.nan, len(basal_tau), len(app_tau)

    delta = float(np.nanmean(app_tau) - np.nanmean(basal_tau))
    return delta, len(basal_tau), len(app_tau)

def run_surrogates_for_record(record, n_surr=8, base_seed=20260708):
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

    # Observed
    X, taus_obs = build_proxy_and_tau_s(rr)
    delta_obs, n_b, n_a = extract_delta_tau_s(taus_obs, t_hr, basal_start, basal_end, approach_start, approach_end)
    print(f"Observed Δτ_s = {delta_obs:+.6f}  (n_basal={n_b}, n_app={n_a})")

    # Surrogates
    surr_deltas = []
    for k in range(n_surr):
        seed = base_seed + int(record) * 100 + k
        Xs = phase_shuffle_independent(X, seed=seed)
        taus_s = build_tau_s_from_X(Xs)
        # Note: time axis + masks identical; only the bivariate values are phase-shuffled
        d_s, _, _ = extract_delta_tau_s(taus_s, t_hr, basal_start, basal_end, approach_start, approach_end)
        surr_deltas.append(float(d_s))
        print(f"  surr {k+1}/{n_surr}  Δτ_s={d_s:+.6f} (seed={seed})")

    surr_arr = np.array(surr_deltas)
    abs_obs = abs(delta_obs)

    # Two-sided on magnitude
    n_extreme = int(np.sum(np.abs(surr_arr) >= abs_obs))
    p_twosided = n_extreme / n_surr

    # Direction-specific
    if delta_obs > 0:
        n_dir = int(np.sum(surr_arr >= delta_obs))
    else:
        n_dir = int(np.sum(surr_arr <= delta_obs))
    p_dir = n_dir / n_surr

    print(f"Surrogates Δτ_s range: [{surr_arr.min():+.4f}, {surr_arr.max():+.4f}]")
    print(f"p (two-sided |Δ| >= obs) = {p_twosided:.3f}   ({n_extreme}/{n_surr})")
    print(f"p (direction-specific)   = {p_dir:.3f}   ({n_dir}/{n_surr})")

    # Save per-record JSON
    out = {
        "protocol": "CCTP v0.1 surrogates",
        "record": record,
        "n_surr": n_surr,
        "delta_obs": delta_obs,
        "surr_deltas": surr_deltas,
        "p_two_sided_abs": p_twosided,
        "p_direction_specific": p_dir,
        "windows": {
            "basal": [basal_start, basal_end],
            "approach": [approach_start, approach_end]
        },
        "W_TAU": W_TAU,
        "stride": STRIDE,
        "timestamp": datetime.now().isoformat()
    }
    out_path = os.path.join(RES_DIR, f"surrogate_cctp_{record}.json")
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2)
    print(f"Saved {out_path}")

    # Histogram figure
    fig_dir = os.path.join(FIG_ROOT, record)
    os.makedirs(fig_dir, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7, 4.5))
    bins = np.linspace(min(surr_arr.min(), delta_obs) - 0.005,
                       max(surr_arr.max(), delta_obs) + 0.005, 9)
    ax.hist(surr_arr, bins=bins, color="#5c6bc0", edgecolor="white", alpha=0.85, label="surrogates (phase-shuffle)")
    ax.axvline(delta_obs, color="#d62728", lw=2.5, label=f"observed Δτ_s = {delta_obs:+.4f}")
    ax.axvline(0.0, color="gray", lw=1, ls=":", alpha=0.7)
    ax.set_xlabel("Δτ_s (approach − basal)")
    ax.set_ylabel("Count (n=8)")
    title_dir = "↑ rise" if delta_obs > 0 else "↓ drop"
    ax.set_title(f"Record {record} — Surrogate null for Δτ_s\n{title_dir} | p_dir={p_dir:.3f} (n={n_surr})")
    ax.legend(loc="upper right", fontsize=9)
    ax.grid(True, alpha=0.25, axis="y")
    plt.tight_layout()
    fig_path = os.path.join(fig_dir, "15_surrogate_delta_tau.png")
    plt.savefig(fig_path, dpi=150)
    plt.close()
    print(f"Saved {fig_path}")

    return {
        "record": record,
        "delta_obs": delta_obs,
        "surr_deltas": surr_deltas,
        "p_twosided": p_twosided,
        "p_dir": p_dir,
        "json": out_path,
        "fig": fig_path
    }

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--record", default="all", help="Record id: 30 | 35 | all")
    parser.add_argument("--n_surr", type=int, default=8)
    args = parser.parse_args()

    if args.record == "all":
        records = ["30", "35"]
    else:
        records = [args.record]

    print("=" * 70)
    print("CCTP Surrogates — light phase-shuffle on Δτ_s (8 per record)")
    print(f"Records: {records}   n_surr={args.n_surr}")
    print("=" * 70)

    results = []
    for rec in records:
        res = run_surrogates_for_record(rec, n_surr=args.n_surr)
        results.append(res)

    # Final comparative table
    print("\n" + "=" * 70)
    print("COMPARATIVE SURROGATE RESULTS — Δτ_s")
    print("=" * 70)
    print(f"{'Record':<8} {'Δ_obs':>12} {'min_surr':>10} {'max_surr':>10} {'p_|Δ|':>8} {'p_dir':>8}")
    print("-" * 70)
    for r in results:
        da = r["delta_obs"]
        sa = np.array(r["surr_deltas"])
        print(f"{r['record']:<8} {da:>+12.6f} {sa.min():>10.6f} {sa.max():>10.6f} "
              f"{r['p_twosided']:>8.3f} {r['p_dir']:>8.3f}")
    print("-" * 70)
    print("p_|Δ|  = fraction of surrogates with |Δ_surr| >= |Δ_obs| (two-sided magnitude test)")
    print("p_dir  = fraction where surrogate is at least as extreme in the observed direction")
    print("Interpretation: p <= 0.125 (1/8) means the observed Δτ_s is among the most extreme under the null.")
    print("=" * 70)

    # Save a small combined summary
    combined = {
        "protocol": "CCTP v0.1 surrogates (light)",
        "n_surr_per_record": args.n_surr,
        "results": results,
        "timestamp": datetime.now().isoformat()
    }
    combined_path = os.path.join(RES_DIR, "surrogate_cctp_combined.json")
    with open(combined_path, "w") as f:
        json.dump(combined, f, indent=2)
    print(f"\nCombined summary saved: {combined_path}")

if __name__ == "__main__":
    main()
