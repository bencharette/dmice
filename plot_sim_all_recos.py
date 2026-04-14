#!/usr/bin/env python3
"""
plot_sim_all_recos.py

Compares angular accuracy of LineFit, Pivot LineFit, MPEFit, and Pivot MPEFit
on the BLO 200-event binned downgoing muon simulation (MC truth available).

Usage:
    python plot_sim_all_recos.py [--csv PATH] [--out DIR]
"""

import os
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

DEFAULT_CSV = os.path.expanduser("~/dmice_work/output/comparison/sim_all_recos.csv")
DEFAULT_OUT = os.path.expanduser("~/dmice_work/output/comparison")

RECOS = [
    ("ic_lf_ang_err_deg",       "IC LineFit",       "steelblue",  "o",  "-"),
    ("pivot_lf_ang_err_deg",    "Pivot LineFit",    "darkorange", "s",  "-"),
    ("mpe_ang_err_deg",         "MPEFit",           "seagreen",   "^",  "-"),
    ("pivot_mpe_ang_err_deg",   "Pivot MPEFit",     "crimson",    "D",  "-"),
]

def median_iqr(vals):
    return np.median(vals), np.percentile(vals, 25), np.percentile(vals, 75)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--csv", default=DEFAULT_CSV)
    parser.add_argument("--out", default=DEFAULT_OUT)
    args = parser.parse_args()

    os.makedirs(args.out, exist_ok=True)
    df = pd.read_csv(args.csv)
    print(f"Events: {len(df)}")

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    # ── Panel 1: Angular error distribution ──────────────────────────────────
    ax = axes[0]
    bins = np.linspace(0, 90, 46)
    for col, label, color, _, ls in RECOS:
        vals = df[col].dropna().values
        med  = np.median(vals)
        ax.hist(vals, bins=bins, histtype="step", color=color,
                lw=1.8, ls=ls, density=True,
                label=f"{label}  (median={med:.1f}°, n={len(vals)})")
        ax.axvline(med, color=color, lw=1.0, ls="--", alpha=0.6)

    ax.set_xlabel("Angular error vs MC truth (°)")
    ax.set_ylabel("Normalised events / bin")
    ax.set_title("Angular error distribution\n(BLO sim, 200 events)")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # ── Panel 2: Median angular error vs energy ───────────────────────────────
    ax = axes[1]
    log_edges = np.logspace(np.log10(df["mc_energy_GeV"].min() * 0.9),
                            np.log10(df["mc_energy_GeV"].max() * 1.1), 7)

    for col, label, color, marker, ls in RECOS:
        sub = df[["mc_energy_GeV", col]].dropna()
        sub["ebin"] = pd.cut(sub["mc_energy_GeV"], bins=log_edges)
        centers, meds, los, his = [], [], [], []
        for interval in sub["ebin"].cat.categories:
            g = sub[sub["ebin"] == interval][col]
            if len(g) < 3:
                continue
            centers.append(np.sqrt(interval.left * interval.right))
            m, q1, q3 = median_iqr(g.values)
            meds.append(m); los.append(q1); his.append(q3)
        centers = np.array(centers)
        meds    = np.array(meds)
        los     = np.array(los)
        his     = np.array(his)
        ax.semilogx(centers, meds, color=color, marker=marker,
                    ls=ls, lw=1.6, ms=6, label=label)
        ax.fill_between(centers, los, his, color=color, alpha=0.12)

    ax.set_xlabel("Muon energy (GeV)")
    ax.set_ylabel("Median angular error (°)")
    ax.set_title("Angular error vs energy\n(median ± IQR band)")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3, which="both")
    ax.set_ylim(bottom=0)

    # ── Summary table ─────────────────────────────────────────────────────────
    print("\nMedian angular errors:")
    for col, label, _, _, _ in RECOS:
        vals = df[col].dropna().values
        if len(vals):
            print(f"  {label:20s}: {np.median(vals):.2f}°  (n={len(vals)})")

    fig.suptitle(
        "BLO 200-event downgoing muon sim — "
        "LineFit vs MPEFit with/without DM-Ice pivot",
        fontsize=11
    )
    plt.tight_layout()

    out_path = os.path.join(args.out, "sim_all_recos_comparison.png")
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"\nSaved: {out_path}")


if __name__ == "__main__":
    main()
