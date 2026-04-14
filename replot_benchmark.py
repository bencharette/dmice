"""Regenerate benchmark plots from existing CSV (no IceTray needed)."""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

CSV  = os.path.expanduser("~/dmice_work/output/splinempe_pivot_comparison.csv")
OUT  = os.path.expanduser("~/dmice_work/output/splinempe_pivot_comparison")

df = pd.read_csv(CSV)
print(f"Loaded {len(df)} events from {CSV}")

abins = np.linspace(0, 15, 61)

# ── Figure: MPEFit vs SPEFit — before and after DM-Ice combined likelihood ────

# ── Figure: MPEFit vs SPEFit (std + pivot) per energy bin ────────────────────
style_4 = [
    ("mpe_std_ang_err", "MPEFit (std)",   "forestgreen", "--"),
    ("mpe_piv_ang_err", "MPEFit (pivot)", "red",         "-"),
    ("spe_std_ang_err", "SPEFit (std)",   "steelblue",   "--"),
    ("spe_piv_ang_err", "SPEFit (pivot)", "darkorange",  "-"),
]

bins_list = sorted(df.bin_id.dropna().unique())

def plot_comparison(pairs, title, suffix):
    """One figure with one panel per energy bin, showing the given (col, label, color, ls) pairs."""
    fig, axes = plt.subplots(1, len(bins_list), figsize=(4*len(bins_list), 4), sharey=True)
    if len(bins_list) == 1:
        axes = [axes]
    for ax, bid in zip(axes, bins_list):
        sub = df[df.bin_id == bid]
        e_med = sub.mc_energy_GeV.median()
        for col, label, color, ls in pairs:
            if col not in df: continue
            vals = sub[col].dropna()
            if len(vals) > 3:
                ax.hist(vals[vals <= 15], bins=abins, histtype="step", lw=2, ls=ls,
                        color=color, density=True,
                        label=f"{label} {vals.median():.2f}°")
        ax.set_title(f"Bin {int(bid)}\n{e_med:.0f} GeV (n={len(sub)})", fontsize=9)
        ax.set_xlabel("Angular error (°)"); ax.legend(fontsize=8); ax.grid(True, alpha=0.3)
    axes[0].set_ylabel("Normalised events / bin")
    fig.suptitle(title, fontsize=11)
    plt.tight_layout()
    path = OUT + suffix
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved: {path}")

# MPEFit std vs LF-pivot vs MPE-pivot
plot_comparison(
    [("mpe_std_ang_err",  "MPEFit (std)",       "forestgreen", "--"),
     ("mpe_piv_ang_err",  "MPEFit (LF-pivot)",  "red",         "-"),
     ("mpe_piv2_ang_err", "MPEFit (MPE-pivot)", "darkred",     "-.")],
    title="MPEFit — std vs LF-pivot vs MPE-pivot seed — per energy bin",
    suffix="_mpe_vs_pivot.png",
)

# SPEFit std vs pivot
plot_comparison(
    [("spe_std_ang_err", "SPEFit (std)",   "steelblue",  "--"),
     ("spe_piv_ang_err", "SPEFit (pivot)", "darkorange", "-")],
    title="SPEFit — standard seed vs DM-Ice pivot seed — per energy bin",
    suffix="_spe_vs_pivot.png",
)

# All three MPEFit variants head-to-head (overall)
plot_comparison(
    [("mpe_std_ang_err",  "MPEFit (std)",       "forestgreen", "--"),
     ("mpe_piv_ang_err",  "MPEFit (LF-pivot)",  "red",         "-"),
     ("mpe_piv2_ang_err", "MPEFit (MPE-pivot)", "darkred",     "-.")],
    title="MPEFit pivot comparison — LF-pivot vs MPE-pivot seed",
    suffix="_mpe_pivot_comparison.png",
)

# MPEFit before vs after DM-Ice combined likelihood
plot_comparison(
    [("mpe_std_ang_err", "MPEFit (std)",         "forestgreen", "--"),
     ("mpe_piv_ang_err", "MPEFit (LF-pivot)",    "red",         "-"),
     ("mpe_dm_ang_err",  "MPEFit (DM combined)", "darkred",     "-.")],
    title="MPEFit — std vs pivot seed vs DM-Ice combined likelihood",
    suffix="_mpe_before_after_dm.png",
)

# SPEFit before vs after DM-Ice combined likelihood
plot_comparison(
    [("spe_std_ang_err", "SPEFit (std)",         "steelblue",   "--"),
     ("spe_piv_ang_err", "SPEFit (LF-pivot)",    "darkorange",  "-"),
     ("spe_dm_ang_err",  "SPEFit (DM combined)", "saddlebrown", "-.")],
    title="SPEFit — std vs pivot seed vs DM-Ice combined likelihood",
    suffix="_spe_before_after_dm.png",
)

# ── Figure: Energy vs angular error (profile) — all methods ──────────────────
def plot_energy_vs_ang_err(suffix="_energy_vs_ang_err.png"):
    """Energy vs angular error profile for all methods — the postdoc's plot."""
    methods = [
        ("lf_ang_err",       "LineFit (std)",          "gray",        ":"),
        ("mpe_std_ang_err",  "MPEFit (std)",            "forestgreen", "--"),
        ("spe_std_ang_err",  "SPEFit (std)",            "steelblue",   "--"),
        ("piv_lf_ang_err",   "LineFit (DM spatial)",    "purple",      "-"),
        ("mpe_piv_ang_err",  "MPEFit (DM spatial+time)","red",         "-"),
        ("spe_piv_ang_err",  "SPEFit (DM spatial+time)","darkorange",  "-"),
        ("mpe_dm_ang_err",   "MPEFit (DM combined LL)", "darkred",     "-."),
        ("spe_dm_ang_err",   "SPEFit (DM combined LL)", "saddlebrown", "-."),
    ]

    # Use events with DM-Ice hits only for fair comparison
    has_dm = df["has_dm_hit"] == 1 if "has_dm_hit" in df else pd.Series([True]*len(df))
    sub = df[has_dm & df["mc_energy_GeV"].notna()]
    sub = sub[sub["mc_energy_GeV"] > 0].copy()
    sub["log_E"] = np.log10(sub["mc_energy_GeV"])
    ebins = np.linspace(sub["log_E"].min(), sub["log_E"].max(), 10)

    fig, axes = plt.subplots(1, 3, figsize=(15, 5), sharey=False)

    for ax, (groups, title) in zip(axes, [
        ([methods[0], methods[1], methods[3], methods[5]],   "LineFit"),
        ([methods[0], methods[2], methods[3], methods[6]],   "MPEFit"),
        ([methods[0], methods[1], methods[3], methods[4], methods[6]], "SPEFit"),
    ]):
        pass  # placeholder — reuse profile logic below

    # Simpler: one figure, profile per method
    fig, ax = plt.subplots(figsize=(10, 6))
    for col, label, color, ls in methods:
        if col not in sub.columns:
            continue
        vals = sub[[col, "log_E"]].dropna()
        meds, centers = [], []
        for i in range(len(ebins)-1):
            mask = (vals["log_E"] >= ebins[i]) & (vals["log_E"] < ebins[i+1])
            v = vals.loc[mask, col]
            v = v[v <= 15]
            if len(v) >= 5:
                meds.append(np.median(v))
                centers.append(0.5*(ebins[i]+ebins[i+1]))
        if centers:
            ax.plot(centers, meds, color=color, ls=ls, lw=2, marker="o", ms=4, label=label)

    ax.set_xlabel("log$_{10}$(Energy / GeV)")
    ax.set_ylabel("Median angular error (°)")
    ax.set_title("Angular error vs energy — all reconstruction methods (DM-Ice hit events)")
    ax.legend(fontsize=8, ncol=2)
    ax.grid(True, alpha=0.3)
    ax.set_ylim(bottom=0)
    plt.tight_layout()
    path = OUT + suffix
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved: {path}")

plot_energy_vs_ang_err()
