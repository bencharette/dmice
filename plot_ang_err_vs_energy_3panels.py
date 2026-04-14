"""
plot_ang_err_vs_energy_3panels.py

Three separate plots of angular error vs muon energy from the 5000-event BLO sim:
  1. SplineMPE (standard seed)
  2. MPEFit (DM-Ice pivot seed)
  3. Pivot LineFit

Each plot shows:
  - Grey scatter of individual events
  - Binned median line
  - 25th–75th percentile shaded band
  - 10th–90th percentile shaded band (lighter)

Run (no IceTray needed):
  python3 ~/dmice/plot_ang_err_vs_energy_3panels.py
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

CSV = os.path.expanduser("~/dmice_work/output/splinempe_pivot_comparison.csv")
OUT = os.path.expanduser("~/dmice_work/output")

# ── Load CSV ──────────────────────────────────────────────────────────────────
data = np.genfromtxt(CSV, delimiter=",", names=True, dtype=None, encoding=None)
energy  = data["mc_energy_GeV"].astype(float)
log_e   = np.log10(energy)

def safe_col(col):
    vals = data[col].astype(str)
    out  = np.full(len(vals), np.nan)
    for i, v in enumerate(vals):
        try:
            out[i] = float(v)
        except (ValueError, TypeError):
            pass
    return out

mpe_std  = safe_col("mpe_std_ang_err")   # SplineMPE standard
mpe_piv  = safe_col("mpe_piv_ang_err")   # MPEFit pivot
piv_lf   = safe_col("piv_lf_ang_err")    # Pivot LineFit

# ── Energy bins (log-spaced, 20 bins across 100 GeV–100 TeV) ─────────────────
log_edges = np.linspace(np.log10(100), np.log10(1e5), 21)
log_centers = 0.5 * (log_edges[:-1] + log_edges[1:])
e_centers = 10 ** log_centers

def bin_stats(err_col):
    medians, p25, p75, p10, p90, ns = [], [], [], [], [], []
    for lo, hi in zip(log_edges[:-1], log_edges[1:]):
        mask = (log_e >= lo) & (log_e < hi) & np.isfinite(err_col)
        vals = err_col[mask]
        if len(vals) >= 5:
            medians.append(np.median(vals))
            p25.append(np.percentile(vals, 25))
            p75.append(np.percentile(vals, 75))
            p10.append(np.percentile(vals, 10))
            p90.append(np.percentile(vals, 90))
            ns.append(len(vals))
        else:
            for lst in (medians, p25, p75, p10, p90):
                lst.append(np.nan)
            ns.append(0)
    return (np.array(medians), np.array(p25), np.array(p75),
            np.array(p10),  np.array(p90),  np.array(ns))

# ── Plot config ───────────────────────────────────────────────────────────────
spe_piv  = safe_col("spe_piv_ang_err")    # SPEFit pivot

methods = [
    ("SplineMPE  (standard seed)",      mpe_std, "steelblue",   "splinempe_ang_err_vs_energy.png"),
    ("MPEFit  (DM-Ice pivot seed)",      mpe_piv, "firebrick",   "mpefit_pivot_ang_err_vs_energy.png"),
    ("Pivot LineFit  (DM-Ice anchor)",   piv_lf,  "darkorange",  "pivot_lf_ang_err_vs_energy.png"),
    ("SPEFit  (DM-Ice pivot seed)",      spe_piv, "mediumorchid","spefit_pivot_ang_err_vs_energy.png"),
]

for title, col, color, fname in methods:
    med, q25, q75, q10, q90, ns = bin_stats(col)

    # count valid events
    n_valid = int(np.sum(np.isfinite(col)))
    med_overall = np.nanmedian(col)

    fig, ax = plt.subplots(figsize=(7, 5))

    # scatter (individual events, capped at 30°)
    valid = np.isfinite(col)
    ax.scatter(energy[valid], col[valid],
               s=3, alpha=0.08, color=color, rasterized=True, zorder=1)

    # 25–75 band
    valid_bins = np.isfinite(med)
    ax.fill_between(e_centers[valid_bins],
                    q25[valid_bins], q75[valid_bins],
                    color=color, alpha=0.30, label="25th–75th pct", zorder=3)

    # median line
    ax.plot(e_centers[valid_bins], med[valid_bins],
            color=color, lw=2.5, label=f"Median  (overall {med_overall:.2f}°)", zorder=4)

    ax.set_xscale("log")
    ax.set_xlabel("True muon energy  [GeV]", fontsize=13)
    ax.set_ylabel("Angular error  [°]", fontsize=13)
    ax.set_title(title, fontsize=13)
    ax.set_ylim(0, 20)
    ax.set_xlim(100, 1e5)
    ax.set_xticks([100, 1000, 10000, 100000])
    ax.set_xticklabels(["100 GeV", "1 TeV", "10 TeV", "100 TeV"])
    ax.legend(fontsize=10, loc="upper right")
    ax.grid(True, which="both", alpha=0.25, ls="--")
    ax.annotate(f"n = {n_valid} events", xy=(0.03, 0.04),
                xycoords="axes fraction", fontsize=9, color="grey")

    plt.tight_layout()
    out_path = os.path.join(OUT, fname)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved: {out_path}  (n={n_valid}, median={med_overall:.2f}°)")
