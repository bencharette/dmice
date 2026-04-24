#!/usr/bin/env python3
"""
plot_ang_err_vs_energy_detectors.py

Angular error vs muon energy for det1, det2, and det_center side by side.

Usage:
    python3 plot_ang_err_vs_energy_detectors.py
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT_DIR  = os.path.expanduser("~/dmice_work/output")
CSVS = {
    "det2 (edge, str. 88)":       os.path.expanduser("~/dmice_work/output/comparison/sim_all_recos_det2.csv"),
    "det1 (near-centre, str. 87)": os.path.expanduser("~/dmice_work/output/comparison/sim_all_recos_det1.csv"),
    "det_center (true centre)":    os.path.expanduser("~/dmice_work/output/comparison/sim_all_recos_det_center.csv"),
}

RECOS = [
    ("ic_lf_ang_err_deg",     "IC LineFit",    "steelblue",  "o",  "-"),
    ("pivot_lf_ang_err_deg",  "Pivot LineFit", "darkorange", "s",  "-"),
    ("mpe_ang_err_deg",       "MPEFit",        "seagreen",   "^",  "-"),
    ("pivot_mpe_ang_err_deg", "Pivot MPEFit",  "crimson",    "D",  "-"),
]

def median_iqr(vals):
    return np.median(vals), np.percentile(vals, 25), np.percentile(vals, 75)

fig, axes = plt.subplots(1, 3, figsize=(16, 5), sharey=True)
fig.suptitle("Angular error vs muon energy — DM-Ice detector position comparison",
             fontsize=13, fontweight="bold")

for ax, (det_label, csv_path) in zip(axes, CSVS.items()):
    df = pd.read_csv(csv_path)
    log_edges = np.logspace(np.log10(df["mc_energy_GeV"].min() * 0.9),
                            np.log10(df["mc_energy_GeV"].max() * 1.1), 7)

    for col, label, color, marker, ls in RECOS:
        sub = df[["mc_energy_GeV", col]].dropna()
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
                    ls=ls, lw=2, ms=7,
                    label=f"{label} ({med_all:.2f}°)")
        ax.fill_between(centers, los, his, color=color, alpha=0.12)

    ax.set_xlabel("Muon energy (GeV)", fontsize=11)
    ax.set_title(det_label, fontsize=11, fontweight="bold")
    ax.legend(fontsize=8, loc="upper right")
    ax.grid(True, alpha=0.3, which="both")
    ax.set_ylim(bottom=0)

axes[0].set_ylabel("Median angular error (°)", fontsize=11)

plt.tight_layout()
out = os.path.join(OUT_DIR, "ang_err_vs_energy_all_detectors.png")
fig.savefig(out, dpi=150, bbox_inches="tight")
plt.close()
print(f"Saved: {out}")
