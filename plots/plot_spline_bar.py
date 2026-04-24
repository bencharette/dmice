#!/usr/bin/env python3
"""
plot_spline_bar.py

Bar chart of median angular error for all reconstructions — det1 SplineMPE run.
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

CSV     = os.path.expanduser("~/dmice_work/output/comparison/sim_all_recos_det1_spline_mpe.csv")
OUT_DIR = os.path.expanduser("~/dmice_work/output")

RECOS = [
    ("ic_lf_ang_err_deg",         "LineFit",               "#aaaaaa"),
    ("pivot_lf_ang_err_deg",      "Pivot LineFit",          "steelblue"),
    ("mpe_ang_err_deg",           "MPEFit",                 "seagreen"),
    ("iter_mpe_ang_err_deg",      "IterMPE\n(3-iter MPE+cap)","darkorange"),
    ("spline_std_ang_err_deg",    "SplineMPE\n(LF seed)",   "#e05050"),
    ("spline_iter_ang_err_deg",   "SplineMPE\n(Iter seed)", "#9b59b6"),
    ("spline_piv_ang_err_deg",    "SplineMPE\n(Pivot seed)","crimson"),
]

df = pd.read_csv(CSV)

labels, medians, q25s, q75s = [], [], [], []
for col, label, _ in RECOS:
    vals = df[col].dropna().values
    labels.append(label)
    medians.append(np.median(vals))
    q25s.append(np.percentile(vals, 25))
    q75s.append(np.percentile(vals, 75))

medians = np.array(medians)
errs_lo = medians - np.array(q25s)
errs_hi = np.array(q75s) - medians

colors = [c for _, _, c in RECOS]

fig, ax = plt.subplots(figsize=(10, 5))

bars = ax.barh(range(len(labels)), medians, color=colors,
               xerr=[errs_lo, errs_hi], error_kw=dict(ecolor="black", capsize=4, lw=1.5),
               height=0.6, zorder=3)

# value labels
for i, (med, lo, hi) in enumerate(zip(medians, errs_lo, errs_hi)):
    ax.text(med + hi + 0.05, i, f"{med:.2f}°", va="center", fontsize=10, fontweight="bold")

# highlight best
best_idx = int(np.argmin(medians))
bars[best_idx].set_edgecolor("gold")
bars[best_idx].set_linewidth(2.5)
ax.text(medians[best_idx] + errs_hi[best_idx] + 0.05,
        best_idx + 0.35, "★ best", fontsize=9, color="goldenrod", fontweight="bold")

ax.set_yticks(range(len(labels)))
ax.set_yticklabels(labels, fontsize=10)
ax.set_xlabel("Median angular error (°)", fontsize=12)
ax.set_title("Reconstruction comparison — det1 (str. 87, 1000 events)\nMedian ± IQR",
             fontsize=13, fontweight="bold")
ax.set_xlim(0, max(medians + errs_hi) * 1.25)
ax.invert_yaxis()
ax.grid(True, axis="x", alpha=0.3)
ax.set_axisbelow(True)

# divider between baseline and SplineMPE group
ax.axhline(3.5, color="gray", lw=0.8, ls="--", alpha=0.5)
ax.text(ax.get_xlim()[1] * 0.97, 3.65, "SplineMPE", ha="right",
        fontsize=9, color="gray", style="italic")

plt.tight_layout()
out = os.path.join(OUT_DIR, "spline_bar_det1_mpe.png")
fig.savefig(out, dpi=150, bbox_inches="tight")
plt.close()
print(f"Saved: {out}")
