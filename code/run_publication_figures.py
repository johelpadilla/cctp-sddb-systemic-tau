#!/usr/bin/env python3
"""
Publication-ready batch stratification + lead-time / detector figures.

Reads machine-written results tables so figures cannot silently diverge from CSV.
Outputs under figures/batch/ and figures/publication/.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

BASE = Path(__file__).resolve().parent.parent
RES = BASE / "results"
FIG = BASE / "figures"
FIG_BATCH = FIG / "batch"
FIG_PUB = FIG / "publication"


def load_csv(path: Path) -> list:
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def fnum(v, default=np.nan):
    try:
        if v in (None, ""):
            return default
        return float(v)
    except (TypeError, ValueError):
        return default


def fig_stratification(cohort_csv: Path, out: Path):
    rows = load_csv(cohort_csv)
    if not rows:
        raise SystemExit(f"empty cohort: {cohort_csv}")
    recs = [r["record"] for r in rows]
    dt = np.array([fnum(r.get("delta_tau")) for r in rows])
    de = np.array([fnum(r.get("delta_excess3")) for r in rows])
    substrates = [r.get("substrate", "unknown") for r in rows]
    colors = {"sinus": "#1f77b4", "AF": "#d62728", "paced": "#ff7f0e", "unknown": "#7f7f7f"}
    c = [colors.get(s, "#7f7f7f") for s in substrates]

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.2))
    x = np.arange(len(recs))
    axes[0].bar(x, dt, color=c, edgecolor="k", lw=0.4)
    axes[0].axhline(0, color="k", lw=0.7)
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(recs)
    axes[0].set_ylabel(r"$\Delta\tau_s$ (approach − basal)")
    axes[0].set_title("Systemic Tau by record (colored by substrate)")
    axes[0].grid(True, axis="y", alpha=0.25)

    axes[1].bar(x, de, color=c, edgecolor="k", lw=0.4)
    axes[1].axhline(0, color="k", lw=0.7)
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(recs)
    axes[1].set_ylabel(r"$\Delta$excess3 (approach − basal)")
    axes[1].set_title("RECD excess3 by record (colored by substrate)")
    axes[1].grid(True, axis="y", alpha=0.25)

    # legend proxies
    from matplotlib.patches import Patch

    handles = [Patch(facecolor=colors[k], edgecolor="k", label=k) for k in ("sinus", "AF", "paced") if k in substrates]
    if handles:
        fig.legend(handles=handles, loc="upper center", ncol=3, frameon=False, bbox_to_anchor=(0.5, 1.02))
    fig.suptitle("CCTP / SDDB analytic cohort — stratified relational deltas", y=1.08, fontsize=12)
    fig.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out}")


def fig_leadtime(lead_csv: Path, cum_csv: Path, out: Path):
    leads = load_csv(lead_csv)
    cum = load_csv(cum_csv) if cum_csv.exists() else []
    metrics = ["tau_s", "excess3", "var", "ar1"]
    labels = {
        "tau_s": r"$\tau_s$",
        "excess3": "excess3",
        "var": "variance",
        "ar1": "AR(1)",
    }
    colors = {
        "tau_s": "#2ca02c",
        "excess3": "#9467bd",
        "var": "#1f77b4",
        "ar1": "#ff7f0e",
    }

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.4))

    # left: lead-time strip / box by metric
    data = []
    labs = []
    cols = []
    for m in metrics:
        vals = [
            fnum(r.get("lead_time_h"))
            for r in leads
            if r.get("metric") == m and fnum(r.get("alarmed")) > 0.5
        ]
        vals = [v for v in vals if np.isfinite(v)]
        data.append(vals if vals else [np.nan])
        labs.append(labels[m])
        cols.append(colors[m])
    bp = axes[0].boxplot(data, tick_labels=labs, patch_artist=True, showfliers=True)
    for patch, col in zip(bp["boxes"], cols):
        patch.set_facecolor(col)
        patch.set_alpha(0.55)
    axes[0].set_ylabel("Lead time to VF (hours)")
    axes[0].set_title("Lead time among detections (abs-z ≥ 2, ≥3 consecutive)")
    axes[0].grid(True, axis="y", alpha=0.25)

    # right: cumulative detection rate vs horizon
    if cum:
        for m in metrics:
            xs, ys = [], []
            for r in cum:
                if r.get("metric") == m:
                    xs.append(fnum(r.get("horizon_h")))
                    ys.append(fnum(r.get("detection_rate")))
            if xs:
                order = np.argsort(xs)
                xs = np.array(xs)[order]
                ys = np.array(ys)[order]
                axes[1].plot(xs, ys, "o-", label=labels[m], color=colors[m], lw=1.8)
        axes[1].set_xlabel("Minimum lead-time horizon (hours)")
        axes[1].set_ylabel("Detection rate (fraction of records)")
        axes[1].set_ylim(-0.05, 1.05)
        axes[1].set_title("Cumulative detection vs time-to-VF")
        axes[1].legend(frameon=False, fontsize=9)
        axes[1].grid(True, alpha=0.25)
    else:
        axes[1].text(0.5, 0.5, "No cumulative table", ha="center", va="center")
        axes[1].set_axis_off()

    fig.suptitle("Detector performance: relational vs classic EWS (SDDB analytic cohort)", fontsize=12)
    fig.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out}")


def fig_concordance_scatter(h2h_csv: Path, out: Path):
    rows = load_csv(h2h_csv)
    dt = np.array([fnum(r["delta_tau"]) for r in rows])
    de = np.array([fnum(r["delta_excess3"]) for r in rows])
    dv = np.array([fnum(r["delta_var"]) for r in rows])
    recs = [r["record"] for r in rows]

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2))
    axes[0].axhline(0, color="k", lw=0.6)
    axes[0].axvline(0, color="k", lw=0.6)
    axes[0].scatter(dt, de, s=55, c="#2ca02c", edgecolors="k", lw=0.4, zorder=3)
    for x, y, lab in zip(dt, de, recs):
        if np.isfinite(x) and np.isfinite(y):
            axes[0].annotate(lab, (x, y), textcoords="offset points", xytext=(4, 3), fontsize=8)
    axes[0].set_xlabel(r"$\Delta\tau_s$")
    axes[0].set_ylabel(r"$\Delta$excess3")
    axes[0].set_title("Relational concordance (τ_s vs excess3)")
    axes[0].grid(True, alpha=0.25)

    # normalize var for visual scale
    dv_n = (dv - np.nanmean(dv)) / (np.nanstd(dv) + 1e-12)
    axes[1].axhline(0, color="k", lw=0.6)
    axes[1].axvline(0, color="k", lw=0.6)
    axes[1].scatter(dt, dv_n, s=55, c="#1f77b4", edgecolors="k", lw=0.4, zorder=3)
    for x, y, lab in zip(dt, dv_n, recs):
        if np.isfinite(x) and np.isfinite(y):
            axes[1].annotate(lab, (x, y), textcoords="offset points", xytext=(4, 3), fontsize=8)
    axes[1].set_xlabel(r"$\Delta\tau_s$")
    axes[1].set_ylabel(r"z-scored $\Delta$var")
    axes[1].set_title("Relational vs classic variance EWS")
    axes[1].grid(True, alpha=0.25)

    fig.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cohort", default=str(RES / "cctp_cohort_stratified.csv"))
    ap.add_argument("--leadtime", default=str(RES / "leadtime_per_record.csv"))
    ap.add_argument("--cum", default=str(RES / "leadtime_cumulative.csv"))
    ap.add_argument("--h2h", default=str(RES / "ews_head2head.csv"))
    args = ap.parse_args()

    FIG_BATCH.mkdir(parents=True, exist_ok=True)
    FIG_PUB.mkdir(parents=True, exist_ok=True)

    fig_stratification(Path(args.cohort), FIG_PUB / "fig_stratified_deltas.png")
    # also mirror into batch for continuity
    fig_stratification(Path(args.cohort), FIG_BATCH / "batch_stratified_deltas.png")
    fig_leadtime(Path(args.leadtime), Path(args.cum), FIG_PUB / "fig_leadtime_detector.png")
    fig_leadtime(Path(args.leadtime), Path(args.cum), FIG_BATCH / "batch_leadtime_detector.png")
    if Path(args.h2h).exists():
        fig_concordance_scatter(Path(args.h2h), FIG_PUB / "fig_ews_concordance.png")
        fig_concordance_scatter(Path(args.h2h), FIG_BATCH / "batch_ews_concordance.png")
    print("Publication figures complete.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
