"""
plot_energy_vs_ang_err.py

3 figures (LineFit / MPEFit / SPEFit) showing median angular error vs energy,
with one curve per DM-Ice information level:

  1. No DM-Ice         — baseline
  2. DM spatial only   — position constraint, no timing
  3. DM spatial+time   — pivot seed (position + NaI hit time)
  4. DM combined LL    — IC Pandel + DM-Ice Gaussian in same minimiser (no re-seed)

Reads: ~/dmice_work/output/splinempe_pivot_comparison.csv
Saves: ~/dmice_work/output/energy_vs_ang_err_{lf,mpe,spe}.png

Run locally (no IceTray needed):
    python3 plot_energy_vs_ang_err.py
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

CSV = os.path.expanduser("~/dmice_work/output/splinempe_pivot_comparison.csv")
OUT = os.path.expanduser("~/dmice_work/output/energy_vs_ang_err")

df = pd.read_csv(CSV)
print(f"Loaded {len(df)} events")

# Use only events with DM-Ice hit for a fair comparison
has_dm = df["has_dm_hit"] == 1
df_dm  = df[has_dm & df["mc_energy_GeV"].notna() & (df["mc_energy_GeV"] > 0)].copy()
df_dm["log_E"] = np.log10(df_dm["mc_energy_GeV"])
print(f"  With DM-Ice hit: {len(df_dm)} events")

ebins   = np.linspace(df_dm["log_E"].min() - 0.05, df_dm["log_E"].max() + 0.05, 9)
centers = 0.5 * (ebins[:-1] + ebins[1:])


def profile(df_in, col):
    """Median angular error per energy bin. Returns (centers, medians, lo, hi)."""
    meds, lo, hi, ctrs = [], [], [], []
    if col not in df_in.columns:
        return [], [], [], []
    for i in range(len(ebins) - 1):
        mask = (df_in["log_E"] >= ebins[i]) & (df_in["log_E"] < ebins[i+1])
        v = df_in.loc[mask, col].dropna()
        v = v[v <= 15]
        if len(v) >= 5:
            meds.append(np.median(v))
            lo.append(np.percentile(v, 25))
            hi.append(np.percentile(v, 75))
            ctrs.append(centers[i])
    return ctrs, meds, lo, hi


def make_panel(ax, curves, title):
    """Plot one panel with multiple curves."""
    for col, label, color, ls in curves:
        ctrs, meds, lo, hi = profile(df_dm, col)
        if not ctrs:
            continue
        ax.plot(ctrs, meds, color=color, ls=ls, lw=2, marker="o", ms=5, label=label)
        ax.fill_between(ctrs, lo, hi, color=color, alpha=0.10)
    ax.set_xlabel("log$_{10}$(Energy / GeV)", fontsize=11)
    ax.set_ylabel("Median angular error (°)", fontsize=11)
    ax.set_title(title, fontsize=12)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    ax.set_ylim(bottom=0)


# ── Figure 1: LineFit progression ─────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(7, 5))
make_panel(ax, [
    ("lf_ang_err",       "LineFit — no DM-Ice",          "gray",   ":"),
    ("piv_spat_ang_err", "LineFit — DM spatial only",     "purple", "--"),
    ("piv_lf_ang_err",   "LineFit — DM spatial + time",   "red",    "-"),
], title="LineFit: incremental DM-Ice information")
plt.tight_layout()
fig.savefig(OUT + "_lf.png", dpi=150, bbox_inches="tight")
plt.close()
print(f"Saved: {OUT}_lf.png")

# ── Figure 2: MPEFit progression ──────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(7, 5))
make_panel(ax, [
    ("mpe_std_ang_err",  "MPEFit — no DM-Ice",                    "gray",       ":"),
    ("mpe_spat_ang_err", "MPEFit — DM spatial seed",               "purple",     "--"),
    ("mpe_piv_ang_err",  "MPEFit — DM spatial+time seed",          "red",        "-"),
    ("mpe_dm_ang_err",   "MPEFit — DM combined LL (no re-seed)",   "darkred",    "-."),
], title="MPEFit: incremental DM-Ice information")
plt.tight_layout()
fig.savefig(OUT + "_mpe.png", dpi=150, bbox_inches="tight")
plt.close()
print(f"Saved: {OUT}_mpe.png")

# ── Figure 3: SPEFit progression ──────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(7, 5))
make_panel(ax, [
    ("spe_std_ang_err",  "SPEFit — no DM-Ice",                    "gray",        ":"),
    ("spe_spat_ang_err", "SPEFit — DM spatial seed",               "steelblue",  "--"),
    ("spe_piv_ang_err",  "SPEFit — DM spatial+time seed",          "darkorange", "-"),
    ("spe_dm_ang_err",   "SPEFit — DM combined LL (no re-seed)",   "saddlebrown","-."),
], title="SPEFit: incremental DM-Ice information")
plt.tight_layout()
fig.savefig(OUT + "_spe.png", dpi=150, bbox_inches="tight")
plt.close()
print(f"Saved: {OUT}_spe.png")

# ── Summary table ─────────────────────────────────────────────────────────────
cols = [
    ("lf_ang_err",       "LineFit (no DM)"),
    ("piv_spat_ang_err", "LineFit (spatial)"),
    ("piv_lf_ang_err",   "LineFit (spatial+time)"),
    ("mpe_std_ang_err",  "MPEFit (no DM)"),
    ("mpe_spat_ang_err", "MPEFit (spatial)"),
    ("mpe_piv_ang_err",  "MPEFit (spatial+time)"),
    ("mpe_dm_ang_err",   "MPEFit (combined LL)"),
    ("spe_std_ang_err",  "SPEFit (no DM)"),
    ("spe_spat_ang_err", "SPEFit (spatial)"),
    ("spe_piv_ang_err",  "SPEFit (spatial+time)"),
    ("spe_dm_ang_err",   "SPEFit (combined LL)"),
]
print(f"\n{'Method':<35} {'Median error':>14} {'N':>6}")
print("-" * 58)
for col, label in cols:
    if col not in df_dm.columns:
        print(f"  {label:<33} {'(missing)':>14}")
        continue
    v = df_dm[col].dropna()
    v = v[v <= 15]
    print(f"  {label:<33} {np.median(v):>12.3f}°  {len(v):>5}")
