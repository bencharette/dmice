#!/usr/bin/env python3
"""
plot_sim_distributions.py — Zenith and energy distributions for BLO binned muon sim.

Usage:
    python plot_sim_distributions.py [--npz PATH] [--out DIR]
"""

import os
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

DEFAULT_NPZ = os.path.expanduser("~/dmice_work/output/muons_binned_200ev.npz")
DEFAULT_OUT = os.path.expanduser("~/dmice_work/output/200bin_simplots")

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--npz", default=DEFAULT_NPZ)
    parser.add_argument("--out", default=DEFAULT_OUT)
    args = parser.parse_args()

    os.makedirs(args.out, exist_ok=True)

    d = np.load(args.npz, allow_pickle=True)
    energy_GeV = d["energy_GeV"]
    zenith_blo = d["zenith_rad"]           # BLO: 0=up, π=down
    bin_id     = d["bin_id"]
    bin_edges  = d["bin_edges"]            # GeV

    # Convert BLO zenith to standard (0° = straight down from vertical)
    zen_std_deg = 180.0 - np.degrees(zenith_blo)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # ── Left: Zenith distribution ─────────────────────────────────────────────
    zen_bins = np.linspace(0, 70, 22)
    ax1.hist(zen_std_deg, bins=zen_bins, color="steelblue", edgecolor="white", lw=0.5)
    ax1.set_xlabel("Zenith angle from vertical (°)")
    ax1.set_ylabel("Events / bin")
    ax1.set_title("Zenith distribution")
    ax1.grid(True, alpha=0.3)

    # ── Right: Energy distribution ────────────────────────────────────────────
    log_edges = np.linspace(np.log10(energy_GeV.min() * 0.9),
                            np.log10(energy_GeV.max() * 1.1), 30)
    e_bins = 10 ** log_edges
    ax2.hist(energy_GeV, bins=e_bins, color="steelblue", edgecolor="white", lw=0.5)
    ax2.set_xscale("log")
    ax2.set_xlabel("Muon energy (GeV)")
    ax2.set_ylabel("Events / bin")
    ax2.set_title("Energy distribution")
    ax2.grid(True, alpha=0.3, which="both")

    fig.suptitle(f"BLO 200-event binned downgoing muons  ({len(energy_GeV)} events)", fontsize=12)
    plt.tight_layout()

    out_path = os.path.join(args.out, "sim_distributions.png")
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out_path}")

if __name__ == "__main__":
    main()
