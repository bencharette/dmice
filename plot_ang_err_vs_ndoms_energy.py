#!/usr/bin/env python3
"""
Plot angular error vs n_doms and vs energy from phase1 simulation CSV.
Compares IC-only analytic LineFit vs DM-Ice Pivot (iterative) LineFit.
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import os
CSV          = os.environ.get("CSV",          os.path.expanduser("~/dmice/output/blo_linefit_results.csv"))
OUT_COMBINED = os.environ.get("OUT_COMBINED", os.path.expanduser("~/dmice/output/ang_err_combined.png"))

df = pd.read_csv(CSV)
print(f"Total events: {len(df)}")
print(f"Columns: {list(df.columns)}")

# Drop rows with NaN in key columns
df = df.rename(columns={"ic_ang_err_deg": "ic_analytic_ang_err_deg"})
df = df.dropna(subset=["ic_analytic_ang_err_deg", "cfit_iter_ang_err_deg",
                        "n_doms", "mc_energy_GeV"])


def median_and_iqr(group_vals):
    return np.median(group_vals), np.percentile(group_vals, 25), np.percentile(group_vals, 75)


# ── Plot 1: Angular error vs n_doms ──────────────────────────────────────────
bins_ndoms = [0, 20, 50, 100, 200, 350, 500, 750, 1000]
labels_ndoms = [f"{bins_ndoms[i]}–{bins_ndoms[i+1]}" for i in range(len(bins_ndoms)-1)]
df["ndoms_bin"] = pd.cut(df["n_doms"], bins=bins_ndoms, labels=labels_ndoms)

ic_med, ic_lo, ic_hi = [], [], []
piv_med, piv_lo, piv_hi = [], [], []
bin_centers = []
bin_counts = []

for label in labels_ndoms:
    sub = df[df["ndoms_bin"] == label]
    if len(sub) < 5:
        continue
    lo = bins_ndoms[labels_ndoms.index(label)]
    hi = bins_ndoms[labels_ndoms.index(label) + 1]
    bin_centers.append((lo + hi) / 2)
    bin_counts.append(len(sub))

    m, q1, q3 = median_and_iqr(sub["ic_analytic_ang_err_deg"])
    ic_med.append(m); ic_lo.append(q1); ic_hi.append(q3)

    m, q1, q3 = median_and_iqr(sub["cfit_iter_ang_err_deg"])
    piv_med.append(m); piv_lo.append(q1); piv_hi.append(q3)

bin_centers = np.array(bin_centers)
ic_med = np.array(ic_med); ic_lo = np.array(ic_lo); ic_hi = np.array(ic_hi)
piv_med = np.array(piv_med); piv_lo = np.array(piv_lo); piv_hi = np.array(piv_hi)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

# ── Left: Angular error vs n_doms ─────────────────────────────────────────────
ax1.plot(bin_centers, ic_med, 'o-', color='steelblue', label='IC-only LineFit')
ax1.fill_between(bin_centers, ic_lo, ic_hi, alpha=0.2, color='steelblue')
ax1.plot(bin_centers, piv_med, 's-', color='darkorange', label='DM-Ice Pivot (iterative)')
ax1.fill_between(bin_centers, piv_lo, piv_hi, alpha=0.2, color='darkorange')
ax1.set_xlabel("Number of hit DOMs")
ax1.set_ylabel("Angular error vs MC truth (deg)")
ax1.set_title("Angular error vs N$_{DOMs}$")
ax1.legend()
ax1.set_ylim(0, 20)
ax1.grid(True, alpha=0.3)
for xc, n in zip(bin_centers, bin_counts):
    ax1.text(xc, 5, str(n), ha='center', va='bottom', fontsize=7, color='gray')

# ── Right: Angular error vs energy ────────────────────────────────────────────
log_edges = np.logspace(np.log10(df["mc_energy_GeV"].min()),
                        np.log10(df["mc_energy_GeV"].max()), 9)
df["energy_bin"] = pd.cut(df["mc_energy_GeV"], bins=log_edges)

ic_med, ic_lo, ic_hi = [], [], []
piv_med, piv_lo, piv_hi = [], [], []
bin_centers_e = []
bin_counts_e = []

for interval in df["energy_bin"].cat.categories:
    sub = df[df["energy_bin"] == interval]
    if len(sub) < 5:
        continue
    center = np.sqrt(interval.left * interval.right)
    bin_centers_e.append(center)
    bin_counts_e.append(len(sub))
    m, q1, q3 = median_and_iqr(sub["ic_analytic_ang_err_deg"])
    ic_med.append(m); ic_lo.append(q1); ic_hi.append(q3)
    m, q1, q3 = median_and_iqr(sub["cfit_iter_ang_err_deg"])
    piv_med.append(m); piv_lo.append(q1); piv_hi.append(q3)

bin_centers_e = np.array(bin_centers_e)
ic_med = np.array(ic_med); ic_lo = np.array(ic_lo); ic_hi = np.array(ic_hi)
piv_med = np.array(piv_med); piv_lo = np.array(piv_lo); piv_hi = np.array(piv_hi)

ax2.semilogx(bin_centers_e, ic_med, 'o-', color='steelblue', label='IC-only LineFit')
ax2.fill_between(bin_centers_e, ic_lo, ic_hi, alpha=0.2, color='steelblue')
ax2.semilogx(bin_centers_e, piv_med, 's-', color='darkorange', label='DM-Ice Pivot (iterative)')
ax2.fill_between(bin_centers_e, piv_lo, piv_hi, alpha=0.2, color='darkorange')
ax2.set_xlabel("Muon energy (GeV)")
ax2.set_ylabel("Angular error vs MC truth (deg)")
ax2.set_title("Angular error vs Energy")
ax2.legend()
ax2.set_ylim(0, 20)
ax2.grid(True, alpha=0.3, which='both')
for xc, n in zip(bin_centers_e, bin_counts_e):
    ax2.text(xc, 5, str(n), ha='center', va='bottom', fontsize=7, color='gray')

fig.suptitle("BLO 200-event binned downgoing muons — IC-only vs DM-Ice Pivot LineFit", fontsize=12)
plt.tight_layout()
plt.savefig(OUT_COMBINED, dpi=150)
plt.close()
print(f"Saved: {OUT_COMBINED}")
