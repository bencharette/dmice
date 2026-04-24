#!/usr/bin/env python3
"""
plot_smt8_efficiency.py — SMT8 HLC trigger efficiency vs muon energy.

Reads a repacked simulation npz and plots trigger efficiency per energy bin
with binomial error bars.

Usage:
    python3 ~/dmice/plot_smt8_efficiency.py [npz_path] [--output out.png]

Defaults:
    npz_path : ~/dmice_work/output/muons_binned_1000ev_repacked.npz
    output   : ~/dmice_work/output/smt8_efficiency_vs_energy.png
"""

import os
import sys
import argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument("npz", nargs="?",
                    default=os.path.expanduser(
                        "~/dmice_work/output/muons_binned_1000ev_repacked.npz"),
                    help="Path to simulation npz (repacked or original)")
parser.add_argument("--output", default=None,
                    help="Output PNG path (default: same dir as npz)")
args = parser.parse_args()

npz_path = os.path.expanduser(args.npz)
out_path = args.output or os.path.join(
    os.path.dirname(npz_path), "smt8_efficiency_vs_energy.png")

d = np.load(npz_path, allow_pickle=True)

bin_edges = d['bin_edges']          # shape (n_bins+1,)
bin_id    = d['bin_id']
smt8      = d['smt8_triggered'].astype(bool)

n_bins = len(bin_edges) - 1
effs, errs, centers, xerrs = [], [], [], []

for b in range(n_bins):
    mask   = bin_id == b
    n_tot  = mask.sum()
    n_trig = smt8[mask].sum()
    eff    = n_trig / n_tot if n_tot else 0.0
    err    = np.sqrt(eff * (1 - eff) / n_tot) if n_tot else 0.0
    lo, hi = bin_edges[b], bin_edges[b + 1]
    cen    = np.sqrt(lo * hi)   # geometric centre
    effs.append(eff * 100)
    errs.append(err * 100)
    centers.append(cen)
    xerrs.append((cen - lo, hi - cen))

centers = np.array(centers)
effs    = np.array(effs)
errs    = np.array(errs)
xerr_lo = np.array([x[0] for x in xerrs])
xerr_hi = np.array([x[1] for x in xerrs])

n_per_bin = int((bin_id == 0).sum())

fig, ax = plt.subplots(figsize=(7, 4.5))
ax.errorbar(centers / 1e3, effs, yerr=errs,
            xerr=[xerr_lo / 1e3, xerr_hi / 1e3],
            fmt='o', color='steelblue', capsize=4, capthick=1.5,
            elinewidth=1.5, markersize=7, label=f'{n_per_bin} events/bin')
ax.plot(centers / 1e3, effs, '-', color='steelblue', alpha=0.5)

ax.axhline(100, color='gray', lw=0.8, ls='--', alpha=0.6)
ax.set_xscale('log')
ax.set_xlabel('Muon energy [TeV]', fontsize=12)
ax.set_ylabel('SMT8 trigger efficiency [%]', fontsize=12)
ax.set_title('SMT8 HLC Trigger Efficiency vs Muon Energy\n'
             '(downgoing, 0–60°, aimed through DM-Ice)', fontsize=11)
ax.set_ylim(-5, 110)
ax.set_xlim(bin_edges[0] / 1e3 * 0.6, bin_edges[-1] / 1e3 * 1.5)
ax.set_yticks([0, 20, 40, 60, 80, 100])
ax.grid(True, alpha=0.3, which='both')
ax.legend(fontsize=10)

for x, y in zip(centers / 1e3, effs):
    ax.annotate(f'{y:.0f}%', (x, y), textcoords='offset points',
                xytext=(0, 10), ha='center', fontsize=9, color='steelblue')

fig.tight_layout()
fig.savefig(out_path, dpi=150)
print(f'Saved: {out_path}')
