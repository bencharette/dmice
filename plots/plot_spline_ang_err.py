#!/usr/bin/env python3
"""
plot_spline_ang_err.py

Angular error vs muon energy for det1, comparing SplineMPE variants
against LineFit, PivotLineFit, MPEFit, and IterMPE.

Usage:
    python3 ~/dmice/plot_spline_ang_err.py
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

CSV     = os.path.expanduser("~/dmice_work/output/comparison/sim_all_recos_det1_spline.csv")
OUT_DIR = os.path.expanduser("~/dmice_work/output")

RECOS = [
    # (column,                    label,                  color,         marker, ls,    lw,  zorder)
    ("ic_lf_ang_err_deg",        "LineFit",              "grey",        "o",    "--",  1.5, 1),
    ("pivot_lf_ang_err_deg",     "Pivot LineFit",        "steelblue",   "s",    "--",  1.5, 1),
    ("mpe_ang_err_deg",          "MPEFit",               "seagreen",    "^",    "-",   2.0, 2),
    ("iter_mpe_ang_err_deg",     "IterMPE",              "darkorange",  "D",    "-",   2.0, 2),
    ("spline_std_ang_err_deg",   "SplineMPE (LF seed)",  "crimson",     "v",    "-",   2.5, 3),
    ("spline_piv_ang_err_deg",   "SplineMPE (Piv seed)", "purple",      "P",    "-",   2.5, 3),
    ("spline_iter_ang_err_deg",  "SplineMPE (Iter seed)","darkorchid",  "*",    "-",   2.5, 3),
]

def median_iqr(vals):
    return np.median(vals), np.percentile(vals, 25), np.percentile(vals, 75)

df = pd.read_csv(CSV)
print(f"Loaded {len(df)} events")
print(f"  {'Fit':<28}  {'Valid':>5}  {'Median ang err':>14}")

log_edges = np.logspace(np.log10(df["mc_energy_GeV"].min() * 0.9),
                        np.log10(df["mc_energy_GeV"].max() * 1.1), 7)

fig, ax = plt.subplots(figsize=(10, 6))
fig.suptitle("Angular error vs muon energy — det1 (near-centre, str. 87)\nSplineMPE vs Pandel fits",
             fontsize=13, fontweight="bold")

for col, label, color, marker, ls, lw, zorder in RECOS:
    sub = df[["mc_energy_GeV", col]].dropna()
    if sub.empty:
        print(f"  {label:<28}  {'0':>5}  {'n/a':>14}")
        continue

    sub = sub.copy()
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
    med_all = np.median(df[col].dropna())

    ax.semilogx(centers, meds, color=color, marker=marker,
                ls=ls, lw=lw, ms=7, zorder=zorder,
                label=f"{label} ({med_all:.2f}°)")
    ax.fill_between(centers, los, his, color=color, alpha=0.10, zorder=zorder)

    n_valid = sub[col].count() if col in sub else df[col].notna().sum()
    print(f"  {label:<28}  {n_valid:>5}  {med_all:>13.2f}°")

ax.set_xlabel("Muon energy (GeV)", fontsize=12)
ax.set_ylabel("Median angular error (°)", fontsize=12)
ax.legend(fontsize=9, loc="upper right")
ax.grid(True, alpha=0.3, which="both")
ax.set_ylim(bottom=0)

plt.tight_layout()
out = os.path.join(OUT_DIR, "spline_ang_err_vs_energy_det1.png")
fig.savefig(out, dpi=150, bbox_inches="tight")
plt.close()
print(f"\nSaved: {out}")
