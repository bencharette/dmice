#!/usr/bin/env python3
"""
plot_pivot_mpefit_comparison.py

Plots angular difference between standard MPEFit and DM-Ice Pivot MPEFit
for real DM-Ice coincidence events. Also shows LineFit vs Pivot LineFit
for comparison.

Usage:
    python plot_pivot_mpefit_comparison.py [--csv PATH] [--out DIR]
"""

import os
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

DEFAULT_CSV = os.path.expanduser("~/dmice_work/output/comparison/pivot_mpefit_results.csv")
DEFAULT_OUT = os.path.expanduser("~/dmice_work/output/comparison")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--csv", default=DEFAULT_CSV)
    parser.add_argument("--out", default=DEFAULT_OUT)
    args = parser.parse_args()

    os.makedirs(args.out, exist_ok=True)

    df = pd.read_csv(args.csv)
    print(f"Total events: {len(df)}")

    # Events with valid MPE vs PivotMPE comparison
    mpe_valid = df.dropna(subset=["mpe_vs_pivotmpe_ang_diff_deg",
                                   "lf_vs_pivotlf_ang_diff_deg",
                                   "n_hits"])
    mpe_valid = mpe_valid[mpe_valid["n_hits"] > 0]
    print(f"Events with both fits: {len(mpe_valid)}")

    mpe_diff = mpe_valid["mpe_vs_pivotmpe_ang_diff_deg"].values
    lf_diff  = mpe_valid["lf_vs_pivotlf_ang_diff_deg"].values
    n_hits   = mpe_valid["n_hits"].values
    years    = mpe_valid["year"].astype(str).values

    year_colors = {"2012": "tomato", "2013": "darkorange", "2018": "mediumpurple"}

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    # ── Panel 1: Angular diff distribution — MPE and LF overlaid ─────────────
    ax = axes[0]
    bins = np.linspace(0, 90, 46)

    ax.hist(lf_diff, bins=bins, histtype="stepfilled", color="steelblue",
            alpha=0.45, density=True, label=f"LineFit vs Pivot LF  (median={np.median(lf_diff):.1f}°)")
    ax.hist(mpe_diff, bins=bins, histtype="stepfilled", color="darkorange",
            alpha=0.45, density=True, label=f"MPEFit vs Pivot MPE  (median={np.median(mpe_diff):.1f}°)")
    ax.hist(lf_diff,  bins=bins, histtype="step", color="steelblue",  lw=1.5, density=True)
    ax.hist(mpe_diff, bins=bins, histtype="step", color="darkorange", lw=1.5, density=True)

    ax.axvline(np.median(lf_diff),  color="steelblue",  lw=1.5, ls="--", alpha=0.8)
    ax.axvline(np.median(mpe_diff), color="darkorange", lw=1.5, ls="--", alpha=0.8)

    ax.set_xlabel("Angular difference (°)")
    ax.set_ylabel("Normalised events / bin")
    ax.set_title("Direction shift from DM-Ice pivot\n(real coincidence events)")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # ── Panel 2: MPE angular diff vs n_hits ───────────────────────────────────
    ax = axes[1]
    for yr, col in year_colors.items():
        mask = years == yr
        if mask.sum() == 0:
            continue
        ax.scatter(n_hits[mask], mpe_diff[mask],
                   c=col, s=10, alpha=0.5, label=f"{yr} (n={mask.sum()})", zorder=2)

    # Running median
    sort_idx = np.argsort(n_hits)
    nh_sorted = n_hits[sort_idx]
    md_sorted = mpe_diff[sort_idx]
    window = max(1, len(nh_sorted) // 15)
    med_x, med_y = [], []
    for i in range(0, len(nh_sorted) - window, window // 2):
        med_x.append(np.median(nh_sorted[i:i+window]))
        med_y.append(np.median(md_sorted[i:i+window]))
    ax.plot(med_x, med_y, "k-", lw=2, zorder=5, label="Running median")

    ax.set_xscale("log")
    ax.set_xlabel("Total DOM hits (energy proxy)")
    ax.set_ylabel("MPEFit vs Pivot MPEFit (°)")
    ax.set_title("Pivot shift vs hit multiplicity\n(MPEFit)")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3, which="both")

    # ── Panel 3: LF angular diff vs n_hits ────────────────────────────────────
    ax = axes[2]
    for yr, col in year_colors.items():
        mask = years == yr
        if mask.sum() == 0:
            continue
        ax.scatter(n_hits[mask], lf_diff[mask],
                   c=col, s=10, alpha=0.5, label=f"{yr} (n={mask.sum()})", zorder=2)

    lf_sorted = lf_diff[sort_idx]
    med_y_lf = []
    for i in range(0, len(nh_sorted) - window, window // 2):
        med_y_lf.append(np.median(lf_sorted[i:i+window]))
    ax.plot(med_x, med_y_lf, "k-", lw=2, zorder=5, label="Running median")

    ax.set_xscale("log")
    ax.set_xlabel("Total DOM hits (energy proxy)")
    ax.set_ylabel("LineFit vs Pivot LineFit (°)")
    ax.set_title("Pivot shift vs hit multiplicity\n(LineFit)")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3, which="both")

    fig.suptitle(
        "Effect of DM-Ice pivot on direction reconstruction — "
        "real DM-Ice coincidences 2012/2013/2018",
        fontsize=11
    )
    plt.tight_layout()

    out_path = os.path.join(args.out, "pivot_mpefit_comparison.png")
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out_path}")


if __name__ == "__main__":
    main()
