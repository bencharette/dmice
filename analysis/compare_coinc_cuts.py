#!/usr/bin/env python3
"""
compare_coinc_cuts.py

Compare three event selection strategies on the master coincidence file:
  A. Current cut:   d_perp(LineFit) < D_PERP_MAX  (geometric)
  B. Gaussian cut:  |Δt − μ| < N_SIGMA * σ          (timing-based)
  C. Combined:      A AND B

Also compares what fraction of the broad coincidence window is signal vs accidental.

Run on Cobalt:
  /cvmfs/icecube.opensciencegrid.org/py3-v4.3.0/RHEL_9_x86_64/metaprojects/icetray/v1.12.1/env-shell.sh \\
    python3 -u ~/dmice/compare_coinc_cuts.py [--year 2012]

Output:
  ~/dmice_work/output/coinc_cut_comparison.png
  ~/dmice_work/output/coinc_cut_comparison.csv
"""

import os, sys, math, argparse
import numpy as np

# ── Parameters ───────────────────────────────────────────────────────────────

MU_NS      = 280.0    # NaI timing model offset [ns]
SIGMA_NS   = 81.0     # NaI timing model spread  [ns]
N_SIGMA    = 3.0      # Gaussian cut threshold (3σ = ±243 ns)
D_PERP_MAX = 15.0     # geometric cut [m]
C_M_NS     = 0.2998   # speed of light [m/ns]

DMICE_POS_IC = {
    "det1": np.array([ 31.25,  -72.93, -511.05]),
    "det2": np.array([-334.80, -424.50, -511.26]),
}

I3_PATH  = "/data/user/bcharett/dmice_coincidences_2011_2022/all_dmice_coincidences_2011_2022_fixed.i3"
OUT_DIR  = os.path.expanduser("~/dmice_work/output")

# ── Args ─────────────────────────────────────────────────────────────────────

parser = argparse.ArgumentParser()
parser.add_argument("--year", type=int, default=None,
                    help="Restrict to a single year (e.g. 2012). Default: all years.")
parser.add_argument("--n-sigma", type=float, default=N_SIGMA)
args = parser.parse_args()
N_SIGMA = args.n_sigma

# ── IceTray imports ───────────────────────────────────────────────────────────

from icecube import icetray, dataio, dataclasses

# ── Helpers ───────────────────────────────────────────────────────────────────

IC_STRINGS   = set(range(1, 87))
MUON_STREAMS = {'', 'in_ice', 'InIceSplit'}

# Reconstruction key priority lists — try best first
MPE_KEYS = ["MPEFit", "PoleMuonLlhFit"]
LF_KEYS  = ["LineFit", "PoleMuonLinefit"]

def get_track(frame, keys):
    """Return first OK particle from keys list, or None."""
    for k in keys:
        if k in frame:
            p = frame[k]
            if hasattr(p, "fit_status") and p.fit_status == dataclasses.I3Particle.FitStatus.OK:
                return p, k
    return None, None

def d_perp_and_tgeo(lf, dm_pos):
    """Return (d_perp [m], t_geo [ns]) for LineFit track → DM-Ice position."""
    r     = dm_pos - np.array([lf.pos.x, lf.pos.y, lf.pos.z])
    d_hat = np.array([lf.dir.x, lf.dir.y, lf.dir.z])
    s     = float(np.dot(r, d_hat))
    perp2 = max(0.0, float(np.dot(r, r)) - s**2)
    d_perp = math.sqrt(perp2)
    t_pca  = lf.time + s / C_M_NS
    # For NaI: d_perp is always ~0 in reality, but use full formula for LineFit estimate
    if d_perp < 0.1:
        t_geo = t_pca
    else:
        # geometric Cherenkov propagation (only used for LineFit quality estimate)
        t_geo = t_pca + d_perp / (C_M_NS * math.sin(math.acos(1.0/1.3195)))
    return d_perp, t_geo

# ── Process i3 file ───────────────────────────────────────────────────────────

print(f"Reading: {I3_PATH}")
if args.year:
    print(f"Filtering to year: {args.year}")
print(f"Gaussian cut: |Δt − {MU_NS:.0f}| < {N_SIGMA:.0f}σ = ±{N_SIGMA*SIGMA_NS:.0f} ns")
print(f"Geometric cut: d_perp < {D_PERP_MAX:.0f} m")

records  = []
n_total  = 0
n_no_reco = 0
n_no_dm  = 0
seen     = set()

f = dataio.I3File(I3_PATH)
while f.more():
    frame = f.pop_frame()
    if frame.Stop != icetray.I3Frame.Physics:
        continue

    hdr    = frame["I3EventHeader"]
    stream = getattr(hdr, "sub_event_stream", "")
    if stream not in MUON_STREAMS:
        continue

    year = hdr.start_time.utc_year
    if args.year and year != args.year:
        continue

    # Dedup
    key = (hdr.run_id, hdr.event_id, stream)
    if key in seen:
        continue
    seen.add(key)
    n_total += 1

    # Best reconstruction: MPEFit preferred, LineFit fallback
    mpe, mpe_key = get_track(frame, MPE_KEYS)
    lf,  lf_key  = get_track(frame, LF_KEYS)
    best = mpe if mpe is not None else lf
    best_key = mpe_key if mpe is not None else lf_key

    if best is None:
        n_no_reco += 1
        records.append(dict(year=year, reco_used="none",
                            d_perp_lf=np.nan, delta_t_lf=np.nan,
                            d_perp_mpe=np.nan, delta_t_mpe=np.nan,
                            zen_deg=np.nan, coinc_time_us=np.nan))
        continue

    # DM-Ice hit time
    if "DMIce_detection_time" not in frame:
        n_no_dm += 1
        records.append(dict(year=year, reco_used=best_key,
                            d_perp_lf=np.nan, delta_t_lf=np.nan,
                            d_perp_mpe=np.nan, delta_t_mpe=np.nan,
                            zen_deg=np.degrees(best.dir.zenith), coinc_time_us=np.nan))
        continue

    det_str = str(frame["DMIce_detector"]) if "DMIce_detector" in frame else "det1"
    det_key = "det1" if "det1" in det_str else "det2"
    dm_pos  = DMICE_POS_IC[det_key]

    event_start_daq = hdr.start_time.utc_daq_time
    dm_t_ns  = (frame["DMIce_detection_time"].value - event_start_daq) * 0.1
    coinc_us = dm_t_ns / 1000.0

    # Compute d_perp and Δt for both LF and MPE independently
    dp_lf, dt_lf   = (d_perp_and_tgeo(lf,  dm_pos) if lf  is not None
                      else (np.nan, np.nan))
    dp_mpe, dt_mpe = (d_perp_and_tgeo(mpe, dm_pos) if mpe is not None
                      else (np.nan, np.nan))

    delta_t_lf  = dm_t_ns - dt_lf  if not np.isnan(dt_lf)  else np.nan
    delta_t_mpe = dm_t_ns - dt_mpe if not np.isnan(dt_mpe) else np.nan

    records.append(dict(
        year          = year,
        detector      = det_key,
        reco_used     = best_key,
        has_mpe       = mpe is not None,
        d_perp_lf     = dp_lf,
        delta_t_lf    = delta_t_lf,
        d_perp_mpe    = dp_mpe,
        delta_t_mpe   = delta_t_mpe,
        zen_deg       = np.degrees(best.dir.zenith),
        coinc_time_us = coinc_us,
    ))

f.close()
print(f"\nProcessed {n_total} events  (no reco: {n_no_reco}, no DM time: {n_no_dm})")

# ── Analysis ──────────────────────────────────────────────────────────────────

import pandas as pd
df = pd.DataFrame(records)

n_mpe = df["has_mpe"].sum() if "has_mpe" in df else 0
print(f"Events with MPEFit: {n_mpe}/{n_total} ({100*n_mpe/max(n_total,1):.1f}%)")

# LF-based cuts
df["pass_geom_lf"]  = df["d_perp_lf"] < D_PERP_MAX
df["pass_gauss_lf"] = (df["delta_t_lf"] - MU_NS).abs() < N_SIGMA * SIGMA_NS
df["pass_both_lf"]  = df["pass_geom_lf"] & df["pass_gauss_lf"]

# MPE-based cuts
df["pass_geom_mpe"]  = df["d_perp_mpe"] < D_PERP_MAX
df["pass_gauss_mpe"] = (df["delta_t_mpe"] - MU_NS).abs() < N_SIGMA * SIGMA_NS
df["pass_both_mpe"]  = df["pass_geom_mpe"] & df["pass_gauss_mpe"]

n = len(df)

def pct(x): return f"{100*x/n:>7.1f}%"

print(f"\n{'Cut':<40} {'N (LF)':>8} {'N (MPE)':>8}")
print("-" * 60)
print(f"{'All events':<40} {n:>8} {n:>8}")
for label, col_lf, col_mpe in [
    ("Geometric (d⊥<15m)",          "pass_geom_lf",  "pass_geom_mpe"),
    (f"Gaussian (|Δt−280|<{N_SIGMA*SIGMA_NS:.0f}ns)", "pass_gauss_lf", "pass_gauss_mpe"),
    ("Both",                         "pass_both_lf",  "pass_both_mpe"),
]:
    nl = df[col_lf].sum();  nm = df[col_mpe].sum()
    print(f"  {label:<38} {nl:>8} {nm:>8}")

# Overlap for each method
for method, col_geom, col_gauss in [("LF", "pass_geom_lf", "pass_gauss_lf"),
                                      ("MPE","pass_geom_mpe","pass_gauss_mpe")]:
    geom_only  = df[col_geom]  & ~df[col_gauss]
    gauss_only = df[col_gauss] & ~df[col_geom]
    both       = df[col_geom]  & df[col_gauss]
    print(f"\n[{method}] Geometric only: {geom_only.sum()}  |  "
          f"Gaussian only: {gauss_only.sum()}  |  Both: {both.sum()}")

# Per-year breakdown
print("\nPer-year (MPE-based cuts):")
print(f"{'Year':>6} {'N':>6} {'MPE avail':>10} {'Geom':>6} {'Gauss':>6} {'Both':>6}")
for yr, g in df.groupby("year"):
    nm = g["has_mpe"].sum() if "has_mpe" in g else 0
    pg = g["pass_geom_mpe"].sum()
    pb = g["pass_gauss_mpe"].sum()
    pc = g["pass_both_mpe"].sum()
    print(f"{yr:>6} {len(g):>6} {nm:>10} {pg:>6} {pb:>6} {pc:>6}")

# Save CSV
csv_path = os.path.join(OUT_DIR, "coinc_cut_comparison.csv")
df.to_csv(csv_path, index=False)
print(f"\nCSV: {csv_path}")

# ── Plots ─────────────────────────────────────────────────────────────────────

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import norm as sp_norm

fig, axes = plt.subplots(2, 3, figsize=(15, 9))

# 1. Δt distribution: LF vs MPE, with model overlay
ax = axes[0, 0]
dt_lf  = df["delta_t_lf"].dropna()
dt_mpe = df["delta_t_mpe"].dropna()
lo = min(np.percentile(dt_lf, 2) if len(dt_lf) else MU_NS - 8*SIGMA_NS,
         np.percentile(dt_mpe, 2) if len(dt_mpe) else MU_NS - 8*SIGMA_NS)
hi = max(np.percentile(dt_lf, 98) if len(dt_lf) else MU_NS + 8*SIGMA_NS,
         np.percentile(dt_mpe, 98) if len(dt_mpe) else MU_NS + 8*SIGMA_NS)
lo = max(lo, MU_NS - 12*SIGMA_NS); hi = min(hi, MU_NS + 12*SIGMA_NS)
bins = np.linspace(lo, hi, 60)
ax.hist(dt_lf,  bins=bins, density=True, color="steelblue",  alpha=0.55, label=f"LineFit (n={len(dt_lf)})")
ax.hist(dt_mpe, bins=bins, density=True, color="darkorange", alpha=0.65, label=f"MPEFit  (n={len(dt_mpe)})")
x = np.linspace(lo, hi, 400)
ax.plot(x, sp_norm.pdf(x, MU_NS, SIGMA_NS), "r-", lw=2.5, label=f"Model N({MU_NS:.0f},{SIGMA_NS:.0f})")
ax.axvspan(MU_NS - N_SIGMA*SIGMA_NS, MU_NS + N_SIGMA*SIGMA_NS, alpha=0.08, color="red",
           label=f"±{N_SIGMA:.0f}σ = ±{N_SIGMA*SIGMA_NS:.0f} ns")
ax.set_xlabel("Δt = t_DM − t_geo(reco) [ns]")
ax.set_ylabel("Normalised density")
ax.set_title("Timing residual: LineFit vs MPEFit\n(all events with DM-Ice hit)")
ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

# 2. d_perp: LF vs MPE
ax = axes[0, 1]
dp_lf  = df["d_perp_lf"].dropna()
dp_mpe = df["d_perp_mpe"].dropna()
bins_d = np.linspace(0, 200, 60)
ax.hist(dp_lf,  bins=bins_d, density=True, color="steelblue",  alpha=0.55, label=f"LineFit (n={len(dp_lf)})")
ax.hist(dp_mpe, bins=bins_d, density=True, color="darkorange", alpha=0.65, label=f"MPEFit  (n={len(dp_mpe)})")
ax.axvline(D_PERP_MAX, color="red", ls="--", lw=1.5, label=f"d⊥ cut = {D_PERP_MAX:.0f} m")
ax.set_xlabel("d⊥ (m)"); ax.set_ylabel("Normalised density")
ax.set_title("d⊥ distribution: LineFit vs MPEFit")
ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

# 3. Scatter: d_perp_mpe vs Δt_mpe, coloured by cut
ax = axes[0, 2]
has = df["d_perp_mpe"].notna() & df["delta_t_mpe"].notna()
sub = df[has]
c = np.where(sub["pass_both_mpe"],  "green",
    np.where(sub["pass_geom_mpe"],  "steelblue",
    np.where(sub["pass_gauss_mpe"], "darkorange", "lightgray")))
ax.scatter(sub["d_perp_mpe"].clip(0, 150),
           sub["delta_t_mpe"].clip(MU_NS - 6*SIGMA_NS, MU_NS + 6*SIGMA_NS),
           c=c, s=5, alpha=0.5)
ax.axvline(D_PERP_MAX, color="steelblue", lw=1.5, ls="--", label=f"d⊥<{D_PERP_MAX:.0f}m")
ax.axhline(MU_NS + N_SIGMA*SIGMA_NS, color="darkorange", lw=1.5, ls="--")
ax.axhline(MU_NS - N_SIGMA*SIGMA_NS, color="darkorange", lw=1.5, ls="--",
           label=f"|Δt−μ|<{N_SIGMA:.0f}σ")
ax.set_xlabel("d⊥ MPEFit (m)"); ax.set_ylabel("Δt MPEFit [ns]")
ax.set_title("MPE: d⊥ vs Δt\ngreen=both, blue=geom, orange=gauss, gray=none")
ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

# 4. Bar: event counts for LF vs MPE cuts
ax = axes[1, 0]
cuts = ["Geom\n(LF)", "Gauss\n(LF)", "Both\n(LF)", "Geom\n(MPE)", "Gauss\n(MPE)", "Both\n(MPE)"]
cols = ["pass_geom_lf","pass_gauss_lf","pass_both_lf","pass_geom_mpe","pass_gauss_mpe","pass_both_mpe"]
counts = [df[c].sum() for c in cols]
colors = ["steelblue","steelblue","steelblue","darkorange","darkorange","darkorange"]
bars = ax.bar(cuts, counts, color=colors, alpha=0.8)
for b, cnt in zip(bars, counts):
    ax.text(b.get_x() + b.get_width()/2, cnt + 1, str(cnt),
            ha="center", fontsize=10, fontweight="bold")
ax.axvline(2.5, color="gray", lw=1, ls="--")
ax.set_ylabel("Events passing cut")
ax.set_title(f"Event yield: LF-based vs MPE-based cuts\n({n} total)")
ax.grid(axis="y", alpha=0.3)

# 5. Per-year counts (MPE-based)
ax = axes[1, 1]
years = sorted(df["year"].unique())
x = np.arange(len(years)); w = 0.25
gm = [df[df["year"]==yr]["pass_geom_mpe"].sum()  for yr in years]
gb = [df[df["year"]==yr]["pass_gauss_mpe"].sum() for yr in years]
gc = [df[df["year"]==yr]["pass_both_mpe"].sum()  for yr in years]
ax.bar(x-w, gm, w, color="steelblue",  alpha=0.85, label=f"Geom MPE (d⊥<{D_PERP_MAX:.0f}m)")
ax.bar(x,   gb, w, color="darkorange", alpha=0.85, label=f"Gauss MPE (±{N_SIGMA:.0f}σ)")
ax.bar(x+w, gc, w, color="green",      alpha=0.85, label="Both (MPE)")
ax.set_xticks(x); ax.set_xticklabels(years, rotation=30, ha="right")
ax.set_ylabel("Events"); ax.set_title("MPE-based cut yield per year")
ax.legend(fontsize=8); ax.grid(axis="y", alpha=0.3)

# 6. Zenith distribution by MPE cut group
ax = axes[1, 2]
has_zen = df["zen_deg"].notna()
bins_z  = np.linspace(0, 90, 30)
ax.hist(df[has_zen]["zen_deg"],                           bins=bins_z, density=True,
        color="lightgray",   alpha=0.9, label=f"All ({has_zen.sum()})", zorder=1)
ax.hist(df[has_zen & df["pass_geom_mpe"]]["zen_deg"],     bins=bins_z, density=True,
        color="steelblue",   alpha=0.6, label=f"Geom MPE ({df['pass_geom_mpe'].sum()})", zorder=2)
ax.hist(df[has_zen & df["pass_gauss_mpe"]]["zen_deg"],    bins=bins_z, density=True,
        color="darkorange",  alpha=0.6, label=f"Gauss MPE ({df['pass_gauss_mpe'].sum()})", zorder=3)
ax.hist(df[has_zen & df["pass_both_mpe"]]["zen_deg"],     bins=bins_z, density=True,
        color="green",       alpha=0.8, label=f"Both ({df['pass_both_mpe'].sum()})", zorder=4)
ax.set_xlabel("Reconstructed zenith (°)"); ax.set_ylabel("Normalised density")
ax.set_title("Zenith by MPE cut group\n(downgoing = small zenith)")
ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

ng_mpe = df["pass_geom_mpe"].sum()
nb_mpe = df["pass_gauss_mpe"].sum()
nc_mpe = df["pass_both_mpe"].sum()
year_str = str(args.year) if args.year else "all years"
fig.suptitle(
    f"Coincidence cut comparison (MPEFit) — {year_str} ({n} events)\n"
    f"Geom: d⊥<{D_PERP_MAX:.0f}m → {ng_mpe} | "
    f"Gauss: |Δt−{MU_NS:.0f}|<{N_SIGMA:.0f}σ={N_SIGMA*SIGMA_NS:.0f}ns → {nb_mpe} | "
    f"Both: {nc_mpe}",
    fontsize=10
)
plt.tight_layout()

out_png = os.path.join(OUT_DIR, "coinc_cut_comparison.png")
fig.savefig(out_png, dpi=150, bbox_inches="tight")
plt.close()
print(f"Plot: {out_png}")
