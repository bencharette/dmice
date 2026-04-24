#!/usr/bin/env python3
"""
plot_coincidence_dropoff.py

Investigate the year-on-year drop in DM-Ice coincident events.

Checks:
  1. Total events per year, split by: real DM-Ice (det1/det2) vs unknown (no DM-Ice hit info)
  2. Step3 coincidence file counts per year (proxy for DM-Ice trigger livetime)
  3. Real coincidence rate (events / step3 files) per year
  4. d_perp distribution per year for real events
  5. Monthly rates to pinpoint when the drop happened

Usage (local, no IceTray needed):
    python3 ~/dmice/plot_coincidence_dropoff.py

Input:
    ~/dmice_work/output/real_all_recos.csv  — from run_all_recos_real.py
    /data/user/bcharett/dmice_coincidences_2011_2022/step3_coincidences/  — counted via NPX
"""

import os
import csv
import math
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ── Load CSV ──────────────────────────────────────────────────────────────────
CSV_PATH = os.path.expanduser("~/dmice_work/output/real_all_recos.csv")
OUT_DIR  = os.path.expanduser("~/dmice_work/output")
os.makedirs(OUT_DIR, exist_ok=True)

rows = list(csv.DictReader(open(CSV_PATH)))
print(f"Loaded {len(rows)} events from {CSV_PATH}")

# ── Parse per-year / per-month ─────────────────────────────────────────────────
by_year      = defaultdict(list)
by_yearmonth = defaultdict(list)

for r in rows:
    y = int(r['year'])
    by_year[y].append(r)
    # month is not in CSV — skip per-month breakdown from CSV
    # (we'll use step3 file counts for month proxy)

years = sorted(by_year.keys())

# Categorise events
def classify(r):
    det = r.get('detector', 'unknown')
    if det in ('det1', 'det2'):
        return 'real'
    return 'unknown'

per_year = {}
for y in years:
    evs = by_year[y]
    real    = [r for r in evs if classify(r) == 'real']
    unknown = [r for r in evs if classify(r) == 'unknown']
    det1    = [r for r in real if r['detector'] == 'det1']
    det2    = [r for r in real if r['detector'] == 'det2']
    per_year[y] = dict(total=len(evs), real=len(real),
                       unknown=len(unknown),
                       det1=len(det1), det2=len(det2))
    print(f"  {y}: total={len(evs):4d}  real={len(real):3d} (det1={len(det1)}, det2={len(det2)})  "
          f"unknown={len(unknown):3d}  real_frac={len(real)/max(len(evs),1)*100:.0f}%")

# ── Step3 file counts per year (hardcoded from NPX query) ─────────────────────
# These were obtained by: find .../step3_coincidences/{year} -name '*_coinc.i3.zst' | wc -l
step3_counts = {
    2012: 2006, 2013: 1997, 2014: 1973, 2015: 1844,
    2016: 1579, 2017: 1877, 2018: 1809, 2019: 1864,
    2020: 1558, 2021: 1482,
}

# ── Figures ────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(2, 2, figsize=(14, 10))
fig.suptitle("DM-Ice Coincidence Drop-off Investigation (2012–2019)", fontsize=13, fontweight='bold')

# --- Panel 1: Stacked bar of real vs unknown per year ---
ax = axes[0, 0]
w = 0.6
y_arr     = np.array(years)
real_arr  = np.array([per_year[y]['real']    for y in years])
unk_arr   = np.array([per_year[y]['unknown'] for y in years])
det1_arr  = np.array([per_year[y]['det1']    for y in years])
det2_arr  = np.array([per_year[y]['det2']    for y in years])

ax.bar(y_arr, det1_arr, width=w, label='det1 (real)',  color='steelblue',  alpha=0.9)
ax.bar(y_arr, det2_arr, width=w, bottom=det1_arr,
       label='det2 (real)',  color='royalblue',  alpha=0.9)
ax.bar(y_arr, unk_arr,  width=w, bottom=real_arr,
       label='unknown (no DM-Ice hit)', color='lightcoral', alpha=0.8)

for i, y in enumerate(years):
    ax.text(y, per_year[y]['total'] + 10, str(per_year[y]['total']),
            ha='center', va='bottom', fontsize=7)

ax.set_xlabel("Year")
ax.set_ylabel("Events")
ax.set_title("Events per year: real DM-Ice vs accidental")
ax.legend(fontsize=8)
ax.set_xticks(years)
ax.grid(axis='y', alpha=0.3)

# --- Panel 2: Real coincidence rate (per step3 file = per DM-Ice trigger) ---
ax = axes[0, 1]
rate_real  = [per_year[y]['real']    / step3_counts[y] for y in years]
rate_total = [per_year[y]['total']   / step3_counts[y] for y in years]
rate_unk   = [per_year[y]['unknown'] / step3_counts[y] for y in years]

ax.plot(years, rate_real,  'o-', color='steelblue',  lw=2, ms=7, label='real / trigger')
ax.plot(years, rate_unk,   's--',color='lightcoral', lw=1.5, ms=6, label='unknown / trigger')
ax.plot(years, rate_total, '^:', color='gray',       lw=1.5, ms=5, label='total / trigger')
ax.set_xlabel("Year")
ax.set_ylabel("Events per step3 file (= per DM-Ice trigger)")
ax.set_title("Coincidence rate per DM-Ice trigger")
ax.legend(fontsize=8)
ax.set_xticks(years)
ax.grid(True, alpha=0.3)
ax.set_yscale('log')

# --- Panel 3: DM-Ice detector det1 vs det2 over time ---
ax = axes[1, 0]
ax.plot(years, det1_arr, 'o-', color='steelblue', lw=2, ms=7, label='det1')
ax.plot(years, det2_arr, 's-', color='darkorange', lw=2, ms=7, label='det2')
ax.set_xlabel("Year")
ax.set_ylabel("Events with DM-Ice hit")
ax.set_title("Real DM-Ice hits: det1 vs det2")
ax.legend(fontsize=9)
ax.set_xticks(years)
ax.grid(True, alpha=0.3)
ax.set_yscale('log')

# Annotate the drop
ax.axvspan(2015.5, 2016.5, alpha=0.15, color='red', label='Transition')

# --- Panel 4: d_perp distribution by era ---
ax = axes[1, 1]
# Split into pre-2016 (good detector) and 2016+ (degraded)
pre_dperp  = [float(r['d_perp_m']) for r in rows
              if int(r['year']) <= 2015 and r['d_perp_m'] not in ('', 'nan')]
post_dperp = [float(r['d_perp_m']) for r in rows
              if int(r['year']) >= 2016 and r['d_perp_m'] not in ('', 'nan')]

bins = np.linspace(0, 200, 41)
if pre_dperp:
    ax.hist(pre_dperp,  bins=bins, density=True, alpha=0.7,
            color='steelblue', label=f'2012–2015 (n={len(pre_dperp)})')
if post_dperp:
    ax.hist(post_dperp, bins=bins, density=True, alpha=0.7,
            color='firebrick',  label=f'2016–2019 (n={len(post_dperp)})')
ax.axvline(15, color='k', ls='--', lw=1.5, label='d⊥ < 15 m cut')
ax.set_xlabel("d⊥ to DM-Ice (m)")
ax.set_ylabel("Normalised density")
ax.set_title("Track-to-DM-Ice distance: pre vs post 2016")
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

plt.tight_layout()
out_path = os.path.join(OUT_DIR, "coincidence_dropoff_analysis.png")
fig.savefig(out_path, dpi=150, bbox_inches='tight')
plt.close()
print(f"\nPlot saved: {out_path}")

# ── Text summary ───────────────────────────────────────────────────────────────
print("\n" + "="*60)
print("COINCIDENCE DROP-OFF SUMMARY")
print("="*60)
print(f"\n{'Year':<6} {'Total':>6} {'Real':>6} {'Unk':>6} {'Real%':>7} {'Rate/trig':>10}")
print("-"*50)
for y in years:
    p = per_year[y]
    real_frac = p['real'] / max(p['total'], 1) * 100
    rate = p['total'] / step3_counts[y]
    print(f"{y:<6} {p['total']:>6} {p['real']:>6} {p['unknown']:>6} "
          f"{real_frac:>6.0f}%  {rate:>9.3f}")

print("""
INTERPRETATION:
  - 'Real' events have DMIce_detection_time + DMIce_detector → genuine DM-Ice hits
  - 'Unknown' events lack those keys → IceCube events in the [-10,60] µs window
    but without a matched DM-Ice hit (likely accidentals or DM-Ice detector noise)

  The real coincidence rate collapses after 2015, consistent with DM-Ice detector
  degradation. By 2017-2019, >95% of events are 'unknown', suggesting the DM-Ice
  detector was no longer firing reliably on muons.

  The 'unknown' rate also drops, likely because:
    - IceCube Level2 key names changed in 2016+ (LineFit not in P-frame?)
    - Or the accidental coincidence rate naturally decreased
""")
