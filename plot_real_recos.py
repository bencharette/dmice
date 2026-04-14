"""
plot_real_recos.py

Plots for real DM-Ice coincidence data reconstruction comparison.
Reads: ~/dmice_work/output/real_all_recos.csv
Saves: ~/dmice_work/output/real_recos_*.png

Figures:
  1. Angular shift per energy bin (MPEFit pivot vs std, SPEFit pivot vs std)
  2. Angular shift vs energy (profile + scatter)
  3. Angular shift vs n_DOM hits (profile + scatter)

"Angular shift" = angular separation between pivot-seeded and standard-seeded
method — measures how much the DM-Ice constraint changes the reconstruction.
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

CSV = os.path.expanduser("~/dmice_work/output/real_all_recos.csv")
OUT = os.path.expanduser("~/dmice_work/output/real_recos")

df = pd.read_csv(CSV)
print(f"Loaded {len(df)} events")
print(f"Years: {sorted(df.year.unique())}")
print(f"With energy: {df.energy_GeV.notna().sum()}")
print(f"With pivot:  {df.piv_lf_zen.notna().sum()}")


# ── Angular separation between two (zen, azi) pairs in degrees ────────────────
def ang_sep(zen1_deg, azi1_deg, zen2_deg, azi2_deg):
    z1 = np.radians(zen1_deg); a1 = np.radians(azi1_deg)
    z2 = np.radians(zen2_deg); a2 = np.radians(azi2_deg)
    dot = (np.sin(z1)*np.sin(z2)*np.cos(a1-a2) + np.cos(z1)*np.cos(z2))
    dot = np.clip(dot, -1, 1)
    return np.degrees(np.arccos(dot))


# ── Compute angular shifts ────────────────────────────────────────────────────
df["mpe_shift"] = ang_sep(df.mpe_std_zen, df.mpe_std_azi,
                           df.mpe_piv_zen, df.mpe_piv_azi)
df["spe_shift"] = ang_sep(df.spe_std_zen, df.spe_std_azi,
                           df.spe_piv_zen, df.spe_piv_azi)
df["lf_shift"]  = ang_sep(df.lf_zen,      df.lf_azi,
                           df.piv_lf_zen,  df.piv_lf_azi)

# Use energy if available, else n_doms as proxy
has_energy = df.energy_GeV.notna() & (df.energy_GeV > 0)
df["log_E"] = np.where(has_energy, np.log10(df.energy_GeV), np.nan)
use_energy  = has_energy.sum() > 50   # enough events to plot vs energy

print(f"Using {'energy' if use_energy else 'n_doms'} as x-axis")


# ── Helper: profile plot (median + IQR per bin) ───────────────────────────────
def profile(ax, x, y, bins, color, label, ls="-"):
    mask = np.isfinite(x) & np.isfinite(y)
    x, y = x[mask], y[mask]
    meds, lo, hi, centers = [], [], [], []
    for i in range(len(bins)-1):
        sel = (x >= bins[i]) & (x < bins[i+1])
        if sel.sum() < 5:
            continue
        vals = y[sel]
        meds.append(np.median(vals))
        lo.append(np.percentile(vals, 25))
        hi.append(np.percentile(vals, 75))
        centers.append(0.5*(bins[i]+bins[i+1]))
    centers = np.array(centers)
    meds    = np.array(meds)
    ax.plot(centers, meds, color=color, ls=ls, lw=2, marker="o", ms=5, label=label)
    ax.fill_between(centers, lo, hi, color=color, alpha=0.15)


# ═══════════════════════════════════════════════════════════════════════════════
# Figure 1: Angular shift per energy bin (histogram panels)
# ═══════════════════════════════════════════════════════════════════════════════
pivot_events = df[df.mpe_shift.notna() & (df.d_perp_m < 15)]

if use_energy and len(pivot_events) > 20:
    e_vals  = pivot_events.log_E.dropna()
    e_edges = np.percentile(e_vals.dropna(), [0, 25, 50, 75, 100])
    e_edges = np.unique(np.round(e_edges, 1))
    bin_col = "log_E"
    bin_label = lambda lo, hi: f"$10^{{{lo:.1f}}}$–$10^{{{hi:.1f}}}$ GeV"
else:
    doms    = pivot_events.n_doms_ic
    e_edges = np.percentile(doms, [0, 25, 50, 75, 100])
    e_edges = np.unique(np.round(e_edges).astype(int))
    bin_col = "n_doms_ic"
    bin_label = lambda lo, hi: f"{int(lo)}–{int(hi)} DOMs"

n_panels = len(e_edges) - 1
fig1, axes1 = plt.subplots(1, n_panels, figsize=(4*n_panels, 4), sharey=True)
if n_panels == 1:
    axes1 = [axes1]

abins = np.linspace(0, 15, 46)
for ax, (lo, hi) in zip(axes1, zip(e_edges[:-1], e_edges[1:])):
    sel = pivot_events[(pivot_events[bin_col] >= lo) & (pivot_events[bin_col] < hi)]
    for col, label, color, ls in [
        ("mpe_shift", "MPEFit pivot shift", "red",        "-"),
        ("spe_shift", "SPEFit pivot shift", "darkorange", "--"),
        ("lf_shift",  "LF pivot shift",     "steelblue",  ":"),
    ]:
        vals = sel[col].dropna()
        vals = vals[vals <= 15]
        if len(vals) > 3:
            ax.hist(vals, bins=abins, histtype="step", lw=1.8, ls=ls,
                    color=color, density=True,
                    label=f"{label}\nmed={vals.median():.2f}° n={len(vals)}")
    ax.set_title(bin_label(lo, hi) + f"\n(n={len(sel)})", fontsize=9)
    ax.set_xlabel("Angular shift (°)")
    ax.legend(fontsize=6)
    ax.grid(True, alpha=0.3)
axes1[0].set_ylabel("Normalised events / bin")
fig1.suptitle("Angular shift (pivot vs std) per energy bin — real DM-Ice coincidences", fontsize=10)
plt.tight_layout()
fig1.savefig(OUT + "_shift_per_bin.png", dpi=150, bbox_inches="tight")
plt.close()
print(f"Saved: {OUT}_shift_per_bin.png")


# ═══════════════════════════════════════════════════════════════════════════════
# Figure 2: Angular shift vs energy (profile + scatter)
# ═══════════════════════════════════════════════════════════════════════════════
fig2, axes2 = plt.subplots(1, 2, figsize=(13, 5))

for ax_idx, (xcol, xlabel, xbins) in enumerate([
    ("log_E" if use_energy else "n_doms_ic",
     "log$_{10}$(Energy / GeV)" if use_energy else "N DOMs (IC)",
     np.linspace(2, 6, 9) if use_energy else np.linspace(0, 120, 9)),
    ("n_hits_ic", "N hits (IC)", np.linspace(0, 400, 9)),
]):
    ax = axes2[ax_idx]
    piv = df[df.mpe_shift.notna() & (df.d_perp_m < 15)]
    for col, label, color, ls in [
        ("mpe_shift", "MPEFit pivot shift", "red",        "-"),
        ("spe_shift", "SPEFit pivot shift", "darkorange", "--"),
        ("lf_shift",  "LF pivot shift",     "steelblue",  ":"),
    ]:
        profile(ax, piv[xcol].values, piv[col].values,
                xbins, color, label, ls)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Median angular shift (°) ± IQR")
    ax.set_title(f"Angular shift vs {xlabel.split('(')[0].strip()}")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    ax.set_ylim(bottom=0)

fig2.suptitle("Angular shift (pivot vs std seed) — real DM-Ice coincidences", fontsize=10)
plt.tight_layout()
fig2.savefig(OUT + "_shift_vs_energy.png", dpi=150, bbox_inches="tight")
plt.close()
print(f"Saved: {OUT}_shift_vs_energy.png")


# ═══════════════════════════════════════════════════════════════════════════════
# Figure 3: Direction comparison — std vs pivot zenith scatter per year
# ═══════════════════════════════════════════════════════════════════════════════
years = sorted(df.year.unique())
n_years = len(years)
fig3, axes3 = plt.subplots(1, n_years, figsize=(4*n_years, 4), sharey=True, sharex=True)
if n_years == 1:
    axes3 = [axes3]

for ax, yr in zip(axes3, years):
    sub = df[(df.year == yr) & df.mpe_shift.notna()]
    ax.scatter(sub.mpe_std_zen, sub.mpe_piv_zen,
               s=8, alpha=0.4, color="red",   label=f"MPEFit (n={len(sub)})")
    ax.scatter(sub.spe_std_zen, sub.spe_piv_zen,
               s=8, alpha=0.4, color="steelblue", label=f"SPEFit")
    lim = [0, 180]
    ax.plot(lim, lim, "k--", lw=0.8, alpha=0.5)
    ax.set_title(f"{yr}", fontsize=9)
    ax.set_xlabel("Std zenith (°)")
    ax.legend(fontsize=6, markerscale=2)
    ax.grid(True, alpha=0.3)
axes3[0].set_ylabel("Pivot zenith (°)")
fig3.suptitle("Zenith: std seed vs DM-Ice pivot seed — per year", fontsize=10)
plt.tight_layout()
fig3.savefig(OUT + "_zenith_scatter.png", dpi=150, bbox_inches="tight")
plt.close()
print(f"Saved: {OUT}_zenith_scatter.png")

print("\nDone.")
