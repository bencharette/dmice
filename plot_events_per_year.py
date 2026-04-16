#!/usr/bin/env python3
"""
plot_events_per_year.py

Plot DM-Ice coincidence events per year: raw vs IceCube-rate-normalized.

Run on Cobalt (no IceTray needed, just matplotlib):
    python3 ~/dmice/plot_events_per_year.py
"""

import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

OUT_DIR = os.path.expanduser("~/dmice_work/output")
os.makedirs(OUT_DIR, exist_ok=True)

# From the time-diff analysis (0-10 µs signal bin minus accidental baseline)
years = [2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021]

signal     = [855, 972, 728, 571, 372, 442, 427, 430, 350, 316]
ic_rate_hz = [1446, 1402, 1102, 868, 522, 733, 732, 536, 610, 574]
step3      = [2006, 1997, 1973, 1844, 1579, 1877, 1809, 1864, 1558, 1482]

IC_REF = 1000.0   # reference IceCube L2 rate [Hz]
days   = [365, 365, 365, 365, 365, 365, 365, 365, 366, 365]

raw_per_day  = [s/d      for s, d in zip(signal, days)]
norm_per_day = [s * IC_REF / (r * d)
                for s, r, d in zip(signal, ic_rate_hz, days)]

raw_per_year  = signal
norm_per_year = [s * IC_REF / r for s, r in zip(signal, ic_rate_hz)]

# ── Figure ────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(13, 5))
fig.suptitle("DM-Ice Coincidence Events 2012–2021", fontsize=13, fontweight='bold')

x = np.array(years)
w = 0.35

for ax, raw, norm, ylabel, title in [
    (axes[0], raw_per_year,  norm_per_year,
     "Events per year",      "Events per year"),
    (axes[1], raw_per_day,   norm_per_day,
     "Events per day",       "Events per day"),
]:
    bars_raw  = ax.bar(x - w/2, raw,  width=w, label="Raw",
                       color="steelblue", alpha=0.85)
    bars_norm = ax.bar(x + w/2, norm, width=w, label=f"Normalised to {IC_REF:.0f} Hz IC L2",
                       color="darkorange", alpha=0.85)

    # Mean line for normalised
    mean_norm = np.mean(norm)
    ax.axhline(mean_norm, color="darkorange", ls="--", lw=1.5, alpha=0.8,
               label=f"Norm. mean = {mean_norm:.1f}")

    ax.set_xlabel("Year")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.set_xticks(years)
    ax.set_xticklabels(years, rotation=45, ha='right')
    ax.legend(fontsize=8)
    ax.grid(axis='y', alpha=0.3)

    # Annotate IC rate on raw bars
    for xi, r in zip(x, ic_rate_hz):
        ax.text(xi - w/2, 0.02, f"{r:.0f}Hz",
                ha='center', va='bottom', fontsize=6, rotation=90,
                color='white', fontweight='bold')

axes[0].set_ylim(0, max(raw_per_year) * 1.15)
axes[1].set_ylim(0, max(raw_per_day + norm_per_day) * 1.15)

plt.tight_layout()
out = os.path.join(OUT_DIR, "events_per_year_normalized.png")
fig.savefig(out, dpi=150, bbox_inches='tight')
plt.close()
print(f"Saved: {out}")

# Print table
print(f"\n{'Year':<6} {'IC Hz':>7} {'Raw/yr':>8} {'Norm/yr':>9} {'Raw/day':>8} {'Norm/day':>9}")
print("-" * 55)
for y, s, r, rpd, npd, rpy, npy in zip(
        years, signal, ic_rate_hz, raw_per_day, norm_per_day, raw_per_year, norm_per_year):
    print(f"{y:<6} {r:>7.0f} {rpy:>8.0f} {npy:>9.1f} {rpd:>8.2f} {npd:>9.2f}")
print(f"\nNormalised mean: {np.mean(norm_per_day):.2f} events/day  "
      f"std: {np.std(norm_per_day):.2f}")
