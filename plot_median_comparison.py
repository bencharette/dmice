"""
plot_median_comparison.py

Single plot showing median angular error vs energy for 4 methods:
  MPEFit standard, MPEFit pivot, SPEFit standard, SPEFit pivot
Median lines only — no scatter, no bands.
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

CSV = os.path.expanduser("~/dmice_work/output/splinempe_pivot_comparison.csv")
OUT = os.path.expanduser("~/dmice_work/output/median_comparison.png")

data    = np.genfromtxt(CSV, delimiter=",", names=True, dtype=None, encoding=None)
energy  = data["mc_energy_GeV"].astype(float)
log_e   = np.log10(energy)

def safe_col(col):
    out = np.full(len(data), np.nan)
    for i, v in enumerate(data[col].astype(str)):
        try: out[i] = float(v)
        except: pass
    return out

log_edges   = np.linspace(np.log10(100), np.log10(1e5), 21)
log_centers = 0.5 * (log_edges[:-1] + log_edges[1:])
e_centers   = 10 ** log_centers

def medians(col):
    meds = []
    for lo, hi in zip(log_edges[:-1], log_edges[1:]):
        mask = (log_e >= lo) & (log_e < hi) & np.isfinite(col)
        vals = col[mask]
        meds.append(np.median(vals) if len(vals) >= 5 else np.nan)
    return np.array(meds)

methods = [
    ("MPEFit  (std seed)",   "mpe_std_ang_err", "steelblue",   "--"),
    ("MPEFit  (pivot seed)", "mpe_piv_ang_err", "steelblue",   "-"),
    ("SPEFit  (std seed)",   "spe_std_ang_err", "firebrick",   "--"),
    ("SPEFit  (pivot seed)", "spe_piv_ang_err", "firebrick",   "-"),
]

fig, ax = plt.subplots(figsize=(8, 5))

for label, col, color, ls in methods:
    vals = safe_col(col)
    med  = medians(vals)
    n    = int(np.sum(np.isfinite(vals)))
    ok   = np.isfinite(med)
    overall = np.nanmedian(vals)
    ax.plot(e_centers[ok], med[ok],
            color=color, ls=ls, lw=2.5,
            label=f"{label}  ({overall:.2f}°, n={n})")

ax.set_xscale("log")
ax.set_xlabel("True muon energy  [GeV]", fontsize=13)
ax.set_ylabel("Median angular error  [°]", fontsize=13)
ax.set_title("MPEFit vs SPEFit — standard vs DM-Ice pivot seed", fontsize=13)
ax.set_ylim(0, 8)
ax.set_xlim(100, 1e5)
ax.set_xticks([100, 1000, 10000, 100000])
ax.set_xticklabels(["100 GeV", "1 TeV", "10 TeV", "100 TeV"])
ax.legend(fontsize=10, loc="upper right")
ax.grid(True, which="both", alpha=0.25, ls="--")

plt.tight_layout()
fig.savefig(OUT, dpi=150, bbox_inches="tight")
plt.close()
print(f"Saved: {OUT}")
