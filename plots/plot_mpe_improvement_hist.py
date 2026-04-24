#!/usr/bin/env python3
"""
plot_mpe_improvement_hist.py

Histogram of angular error improvement from pivot-seeded IterMPE vs standard
LineFit-seeded IterMPE, across all energy bins of the 5000-event simulation.
Improvement = std_ang_err - piv_ang_err  (positive = pivot helped)
"""

import os, csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

CSV = os.path.expanduser("~/dmice_work/output/splinempe_pivot_comparison.csv")
OUT = os.path.expanduser("~/dmice_work/output/mpe_improvement_hist.png")

# ── Load data ─────────────────────────────────────────────────────────────────

bin_data = {}   # bin_id -> list of (std_err, piv_err, energy_GeV)

with open(CSV) as f:
    for row in csv.DictReader(f):
        if not row["mpe_std_ang_err"] or not row["mpe_piv_ang_err"]:
            continue
        b   = int(row["bin_id"])
        std = float(row["mpe_std_ang_err"])
        piv = float(row["mpe_piv_ang_err"])
        ene = float(row["mc_energy_GeV"])
        bin_data.setdefault(b, []).append((std, piv, ene))

all_bins = sorted(bin_data.keys())

bin_colors = {0: "steelblue", 1: "darkorange", 2: "mediumseagreen",
              3: "orchid", 4: "tomato"}

# ── Per-bin improvement arrays ────────────────────────────────────────────────

improvements = {}
bin_labels   = {}
for b in all_bins:
    rows = bin_data[b]
    imp  = np.array([s - p for s, p, _ in rows])
    improvements[b] = imp
    med_tev = np.median([e for _, _, e in rows]) / 1e3
    bin_labels[b] = f"bin {b}  median {med_tev:.2f} TeV  n={len(imp)}"

# ── Determine x range: use 1st–99th percentile of all data, no hard clip ──────

all_vals = np.concatenate(list(improvements.values()))
x_lo = max(np.percentile(all_vals, 0.5), -90)
x_hi = min(np.percentile(all_vals, 99.5), 90)

# ── Plot: one panel per bin ───────────────────────────────────────────────────

fig, axes = plt.subplots(len(all_bins), 1,
                         figsize=(8, 3.2 * len(all_bins)),
                         sharex=True)

for ax, b in zip(axes, all_bins):
    imp   = improvements[b]
    color = bin_colors.get(b, "gray")

    # Events outside display range
    n_outside = int(np.sum(imp < x_lo) + np.sum(imp > x_hi))
    imp_disp  = imp[(imp >= x_lo) & (imp <= x_hi)]

    bin_edges = np.linspace(x_lo, x_hi, 61)
    ax.hist(imp_disp, bins=bin_edges, color=color, alpha=0.82,
            edgecolor="white", linewidth=0.4)

    med = np.median(imp)
    ax.axvline(0,   color="k",     lw=1.2, ls="--", zorder=5, label="no change")
    ax.axvline(med, color=color,   lw=2.0, ls="-",  zorder=6,
               label=f"median = {med:+.2f}°")

    frac_better  = 100 * np.mean(imp > 0.1)   # pivot meaningfully better
    frac_worse   = 100 * np.mean(imp < -0.1)  # pivot meaningfully worse
    title = (f"{bin_labels[b]}   |   pivot better: {frac_better:.0f}%   "
             f"worse: {frac_worse:.0f}%")
    if n_outside:
        title += f"   [{n_outside} outside range]"
    ax.set_title(title, fontsize=9)
    ax.set_ylabel("Events")
    ax.legend(fontsize=8, framealpha=0.7, loc="upper right")

axes[-1].set_xlabel(
    "Angular error improvement  [deg]\n"
    r"(std IterMPE $\Delta\Psi$  −  pivot IterMPE $\Delta\Psi$,"
    "  positive = pivot seeding helped)"
)

fig.suptitle("IterMPE angular error improvement from pivot seeding\n"
             "5000-event noisy simulation, all energy bins", fontsize=11)
fig.tight_layout()
fig.savefig(OUT, dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"Saved: {OUT}")

# ── Print summary table ───────────────────────────────────────────────────────
print(f"\n{'bin':>3}  {'med E (TeV)':>11}  {'n':>5}  {'median imp':>10}  "
      f"{'% better':>8}  {'% worse':>7}  {'n outside':>9}")
for b in all_bins:
    imp     = improvements[b]
    rows    = bin_data[b]
    med_tev = np.median([e for _, _, e in rows]) / 1e3
    n_out   = int(np.sum(imp < x_lo) + np.sum(imp > x_hi))
    print(f"{b:>3}  {med_tev:>11.2f}  {len(imp):>5}  "
          f"{np.median(imp):>+10.3f}°  "
          f"{100*np.mean(imp>0.1):>7.1f}%  "
          f"{100*np.mean(imp<-0.1):>6.1f}%  "
          f"{n_out:>9}")
