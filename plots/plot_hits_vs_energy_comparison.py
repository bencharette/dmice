#!/usr/bin/env python3
"""
plot_hits_vs_energy_comparison.py

Compares DOM hit multiplicity between:
  - BLO 200-event binned downgoing muon simulation (n_hits vs MC energy)
  - Real DM-Ice coincidence events 2012/2013/2018 (n_dom_hits, no energy)

Usage:
    python plot_hits_vs_energy_comparison.py [--npz PATH] [--pkl PATH] [--out DIR]
"""

import os
import argparse
import pickle
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

DEFAULT_NPZ = os.path.expanduser("~/dmice_work/output/comparison/muons_binned_200ev.npz")
DEFAULT_PKL = os.path.expanduser("~/dmice_work/output/comparison/linefit_all_years.pkl")
DEFAULT_CSV = os.path.expanduser("~/dmice_work/output/comparison/real_hits_energy_v2.csv")
DEFAULT_OUT = os.path.expanduser("~/dmice_work/output/comparison")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--npz", default=DEFAULT_NPZ)
    parser.add_argument("--pkl", default=DEFAULT_PKL)
    parser.add_argument("--csv", default=DEFAULT_CSV)
    parser.add_argument("--out", default=DEFAULT_OUT)
    args = parser.parse_args()

    os.makedirs(args.out, exist_ok=True)

    # ── Load simulation ───────────────────────────────────────────────────────
    d = np.load(args.npz, allow_pickle=True)
    sim_energy  = d["energy_GeV"]
    sim_n_hits  = d["n_hits"]      # total photon hits across all DOMs
    sim_n_doms  = d["n_doms"]      # unique DOMs hit
    sim_bin_id  = d["bin_id"]
    sim_bin_edges = d["bin_edges"]
    n_bins = len(np.unique(sim_bin_id))
    bin_colors = plt.cm.plasma(np.linspace(0.1, 0.85, n_bins))

    # ── Load real data ────────────────────────────────────────────────────────
    with open(args.pkl, "rb") as f:
        real = pickle.load(f)
    real_hits = real["n_dom_hits"].values
    real_years = real["year"].values

    # Load energy CSV (reprocessed TruncatedEnergy)
    import pandas as pd
    edf = pd.read_csv(args.csv)
    edf = edf[edf["energy_GeV"].notna() & (edf["energy_GeV"] > 0)]
    real_e_GeV  = edf["energy_GeV"].values
    real_e_hits = edf["n_hits"].values
    real_e_year = edf["year"].astype(str).values

    SIM_COLOR  = "steelblue"
    REAL_COLOR = "tomato"

    # ── Figure: 2 panels ──────────────────────────────────────────────────────
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # ── Panel 1: hits vs energy scatter ───────────────────────────────────────
    ax = axes[0]
    ax.scatter(sim_energy, sim_n_hits,
               c=SIM_COLOR, s=18, alpha=0.7, label=f"BLO sim (n={len(sim_energy)})", zorder=3)
    ax.scatter(real_e_GeV, real_e_hits,
               c=REAL_COLOR, s=12, alpha=0.6, label=f"Real data (n={len(real_e_GeV)})", zorder=2)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Muon energy (GeV)")
    ax.set_ylabel("Total DOM hits")
    ax.set_title("Hits vs energy")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3, which="both")

    # ── Panel 2: normalised hit count distribution ────────────────────────────
    ax = axes[1]
    bins = np.logspace(0, 6, 45)
    ax.hist(real_hits, bins=bins, histtype="stepfilled", color=REAL_COLOR,
            alpha=0.5, label=f"Real all years (n={len(real_hits)})", density=True)
    ax.hist(sim_n_hits, bins=bins, histtype="stepfilled", color=SIM_COLOR,
            alpha=0.5, label=f"BLO sim (n={len(sim_n_hits)})", density=True)
    ax.hist(real_hits, bins=bins, histtype="step", color=REAL_COLOR, lw=1.5, density=True)
    ax.hist(sim_n_hits, bins=bins, histtype="step", color=SIM_COLOR, lw=1.5, density=True)
    ax.set_xscale("log")
    ax.set_xlabel("Total DOM hits")
    ax.set_ylabel("Normalised events / bin")
    ax.set_title("Hit count distribution (normalised)")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3, which="both")

    fig.suptitle(
        "DOM hits vs energy — BLO sim (downgoing, 100 GeV–100 TeV) vs "
        "real DM-Ice coincidences 2012/2013/2018",
        fontsize=11
    )
    plt.tight_layout()

    out_path = os.path.join(args.out, "hits_vs_energy_comparison.png")
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out_path}")


if __name__ == "__main__":
    main()
