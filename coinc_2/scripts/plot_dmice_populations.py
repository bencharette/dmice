#!/usr/bin/env python3
"""
plot_dmice_populations.py

Run locally. Reads dmice_waveforms.csv produced by extract_dmice_waveforms.py
and reproduces:
  Fig 7.11 — tau vs pulse height scatter plot (muon/alpha/gamma separation)
  Fig 7.15 — energy spectrum (all events vs muons)

Input:  ~/dmice/output/dmice_waveforms.csv
Output: ~/dmice/output/coinc_2/fig7_11_tau_vs_height.png
        ~/dmice/output/coinc_2/fig7_15_energy_spectrum.png
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

CSV = os.path.expanduser("~/dmice/output/dmice_waveforms.csv")
OUTDIR = os.path.expanduser("~/dmice/output/coinc_2/")
os.makedirs(OUTDIR, exist_ok=True)

print("Loading CSV...")
df = pd.read_csv(CSV)
print(f"  {len(df)} total events, {df['is_muon'].sum()} muons")

# Use PMT-1a (pmt1a) for det1 plots — single PMT gives better separation (Hubbard §7.2.1)
det1 = df[(df["detector"] == "det1") & (df["pmt"] == "pmt1a") & (df["hv_era"] == "lowHV")]
det2 = df[(df["detector"] == "det2") & (df["pmt"] == "pmt2a") & (df["hv_era"] == "lowHV")]

# ── Figure 7.11 — tau vs pulse height ────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(13, 5))
fig.suptitle("τ vs Pulse Height — DM-Ice Event Population (lowHV)", fontsize=13)

for ax, sub, title, h_cut, h_psd, tau_cut in [
    (axes[0], det1, "Det-1 (PMT-1a)", 650, 325, 177),
    (axes[1], det2, "Det-2 (PMT-2a)", 400, 324, 184),
]:
    muons  = sub[sub["is_muon"]]
    others = sub[~sub["is_muon"]]

    # Classify non-muons into alphas vs gammas by pulse height
    # Alphas: pulse height > h_psd (but failed tau cut), Gammas: lower height
    alphas = others[others["pulse_height"] > h_psd]
    gammas = others[others["pulse_height"] <= h_psd]

    ax.scatter(gammas["pulse_height"], gammas["tau"], s=0.3, c="blue",  alpha=0.3, label="Gammas",  rasterized=True)
    ax.scatter(alphas["pulse_height"], alphas["tau"], s=0.3, c="red",   alpha=0.3, label="Alphas",  rasterized=True)
    ax.scatter(muons["pulse_height"],  muons["tau"],  s=0.5, c="black", alpha=0.5, label="Muons",   rasterized=True)

    ax.axvline(h_cut,  color="gray", lw=0.8, ls="--", alpha=0.7)
    ax.axvline(h_psd,  color="gray", lw=0.8, ls=":",  alpha=0.7)
    ax.axhline(tau_cut, color="gray", lw=0.8, ls=":",  alpha=0.7)

    ax.set_xlabel("PMT Pulse Height [ADC counts]", fontsize=11)
    ax.set_ylabel("τ [ns]", fontsize=11)
    ax.set_title(title, fontsize=11)
    ax.legend(fontsize=9, markerscale=5)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

fig.tight_layout()
out1 = os.path.join(OUTDIR, "fig7_11_tau_vs_height.png")
fig.savefig(out1, dpi=150)
print(f"Saved: {out1}")
plt.close()

# ── Figure 7.15 — energy spectrum ────────────────────────────────────────────
# Use Det-1 only (Det-2 PMTs saturate in muon regime — Hubbard p.126)
det1_all = df[(df["detector"] == "det1") & (df["pmt"] == "pmt1a") & (df["hv_era"] == "lowHV")]
det1_mu  = det1_all[det1_all["is_muon"]]

# Energy in ADC·bin — convert to keV requires calibration; plot in ADC units for now
fig, ax = plt.subplots(figsize=(9, 5))

bins = np.linspace(0, det1_all["energy"].quantile(0.999), 200)

ax.hist(det1_all["energy"], bins=bins, color="black", histtype="step", lw=1.2, label="All data")
ax.hist(det1_mu["energy"],  bins=bins, color="blue",  histtype="step", lw=1.2, label="Muons")

ax.set_yscale("log")
ax.set_xlabel("Waveform Integral — sum_128 [ADC·bin]", fontsize=11)
ax.set_ylabel("Counts", fontsize=11)
ax.set_title("DM-Ice Det-1 Energy Spectrum (PMT-1a, lowHV)", fontsize=12)
ax.legend(fontsize=10)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

# Mark approximate region boundaries
ymax = ax.get_ylim()[1]
ax.axvline(det1_all["energy"].quantile(0.30), color="blue",  lw=0.8, ls="--", alpha=0.5)
ax.axvline(det1_all["energy"].quantile(0.80), color="blue",  lw=0.8, ls="--", alpha=0.5)
ax.text(det1_all["energy"].quantile(0.05),  ymax*0.5, "Gamma\ndominated", fontsize=8, ha="center", color="gray")
ax.text(det1_all["energy"].quantile(0.55),  ymax*0.5, "Alpha\ndominated",  fontsize=8, ha="center", color="gray")
ax.text(det1_all["energy"].quantile(0.92),  ymax*0.5, "Muon\ndominated",   fontsize=8, ha="center", color="gray")

fig.tight_layout()
out2 = os.path.join(OUTDIR, "fig7_15_energy_spectrum.png")
fig.savefig(out2, dpi=150)
print(f"Saved: {out2}")
plt.close()
