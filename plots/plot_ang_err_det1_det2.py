#!/usr/bin/env python3
"""
plot_ang_err_det1_det2.py

Angular error vs muon energy for det1 and det2, LineFit and Pivot LineFit only.
Reads the 5000-event splinempe_pivot_comparison.csv and merges target_det from the NPZ.

Usage (on Cobalt, no IceTray needed):
    python3 plot_ang_err_det1_det2.py
"""

import os
import math
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

CSV_PATH = os.path.expanduser(
    "~/dmice_work/output/splinempe_pivot_comparison.csv"
)
NPZ_PATH = os.path.expanduser(
    "~/dmice_work/output/muons_binned_5000ev_repacked_injected.npz"
)
OUT_DIR = os.path.expanduser("~/dmice_work/output")

RECOS = [
    ("lf_ang_err",     "IC LineFit",    "steelblue",  "o"),
    ("piv_lf_ang_err", "Pivot LineFit", "darkorange", "s"),
]

N_BINS = 8

def median_iqr(vals):
    return np.median(vals), np.percentile(vals, 25), np.percentile(vals, 75)

df = pd.read_csv(CSV_PATH)

# Attach target_det from NPZ using event_id as index
d = np.load(NPZ_PATH, allow_pickle=True)
target_det = d["target_det"]
df["target_det"] = df["event_id"].apply(lambda i: int(target_det[i]))

det_map = {
    "det1 (str. 87)": 0,
    "det2 (str. 88)": 1,
}

fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
fig.suptitle(
    "Angular error vs muon energy — DM-Ice det1 vs det2 (5000-event sim)",
    fontsize=13, fontweight="bold",
)

for ax, (det_label, det_id) in zip(axes, det_map.items()):
    sub = df[df["target_det"] == det_id].copy()

    # Drop events with no IC hits (n_doms_ic == 0)
    sub = sub[sub["n_doms_ic"] > 0]

    log_edges = np.logspace(
        np.log10(sub["mc_energy_GeV"].min() * 0.9),
        np.log10(sub["mc_energy_GeV"].max() * 1.1),
        N_BINS + 1,
    )

    for col, label, color, marker in RECOS:
        valid = sub[["mc_energy_GeV", col]].dropna()
        valid = valid[valid[col] >= 0]
        valid["ebin"] = pd.cut(valid["mc_energy_GeV"], bins=log_edges)
        centers, meds, los, his = [], [], [], []
        for interval in valid["ebin"].cat.categories:
            g = valid[valid["ebin"] == interval][col]
            if len(g) < 3:
                continue
            centers.append(math.sqrt(interval.left * interval.right))
            m, q1, q3 = median_iqr(g.values)
            meds.append(m)
            los.append(q1)
            his.append(q3)
        centers = np.array(centers)
        meds    = np.array(meds)
        med_all = np.median(valid[col].dropna())
        ax.semilogx(
            centers, meds,
            color=color, marker=marker, ls="-", lw=2, ms=7,
            label=f"{label} ({med_all:.2f}°)",
        )
        ax.fill_between(centers, los, his, color=color, alpha=0.12)

    ax.set_xlabel("Muon energy (GeV)", fontsize=11)
    ax.set_title(det_label, fontsize=11, fontweight="bold")
    ax.legend(fontsize=10, loc="upper right")
    ax.grid(True, alpha=0.3, which="both")
    ax.set_ylim(bottom=0)

axes[0].set_ylabel("Median angular error (°)", fontsize=11)

plt.tight_layout()
out = os.path.join(OUT_DIR, "ang_err_det1_det2_lf_pivot.png")
fig.savefig(out, dpi=150, bbox_inches="tight")
plt.close()
print(f"Saved: {out}")

# Print summary stats
for det_label, det_id in det_map.items():
    sub = df[(df["target_det"] == det_id) & (df["n_doms_ic"] > 0)]
    print(f"\n{det_label} ({len(sub)} events with IC hits):")
    for col, label, _, _ in RECOS:
        vals = sub[col].dropna()
        vals = vals[vals >= 0]
        if len(vals):
            print(f"  {label:<20} median={np.median(vals):.2f}°  n={len(vals)}")
