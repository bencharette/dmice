#!/usr/bin/env python3
"""
score_real_coincidences.py

Apply the DM-Ice timing model to real coincidence events from the
all_dmice_coincidences_2011_2022.i3.zst file.

Sanity check: events where the LineFit track passes close to DM-Ice (d_perp < D_PERP_MAX)
should have Δt consistent with the model Gaussian N(μ, σ).

Run on Cobalt:
  /cvmfs/icecube.opensciencegrid.org/py3-v4.3.0/RHEL_9_x86_64/metaprojects/icetray/v1.12.1/env-shell.sh \
    python3 ~/dmice/score_real_coincidences.py

Output:
  ~/dmice_work/output/real_coincidence_scoring.png
  ~/dmice_work/output/real_coincidence_scoring.csv
"""

import os, math
import numpy as np

# ── Constants ────────────────────────────────────────────────────────────────

C_M_NS   = 0.2998          # speed of light in vacuum [m/ns]
N_ICE    = 1.3195          # ice refractive index at Cherenkov frequency
THETA_C  = math.acos(1.0 / N_ICE)   # ≈ 40.8°

# DM-Ice positions in IceCube coordinates [m]
DMICE_POS_IC = {
    "det1": np.array([ 31.25,  -72.93, -511.05]),
    "det2": np.array([-334.80, -424.50, -511.26]),
}

MODEL_PATH = os.path.expanduser("~/dmice_work/output/dmice_timing_model.npz")
I3_PATH    = "/data/user/bcharett/dmice_coincidences_2011_2022/all_dmice_coincidences_2011_2022.i3.zst"
OUT_DIR    = os.path.expanduser("~/dmice_work/output")
D_PERP_MAX = 15.0   # m — "through-detector" cut

# ── Geometry helpers ──────────────────────────────────────────────────────────

def t_geometric(track_pos, track_dir, t0_ns, dm_pos):
    """Expected first Cherenkov photon arrival at dm_pos [ns], IceCube coords."""
    r      = np.asarray(dm_pos) - np.asarray(track_pos)
    d_hat  = np.asarray(track_dir, dtype=float)
    d_hat  = d_hat / np.linalg.norm(d_hat)
    s      = np.dot(r, d_hat)
    d_perp = math.sqrt(max(0.0, np.dot(r, r) - s**2))
    t_pca  = t0_ns + s / C_M_NS
    if d_perp < 0.01:
        t_geo = t_pca
    else:
        t_geo = t_pca + d_perp / (C_M_NS * math.sin(THETA_C))
    return t_geo, d_perp


# ── Load timing model ─────────────────────────────────────────────────────────

m = np.load(MODEL_PATH, allow_pickle=True)
mu_ns    = float(m["mu_ns"])
sigma_ns = float(m["sigma_ns"])
eps0     = float(m["efficiency"])
print(f"Model: μ={mu_ns:+.1f} ns  σ={sigma_ns:.1f} ns  ε={eps0:.3f}  d⊥_max={D_PERP_MAX} m")

# ── Process I3 file ───────────────────────────────────────────────────────────

from icecube import icetray, dataio, dataclasses, recclasses, simclasses

TRACK_KEY = "LineFit"   # could also try MPEFit

records = []
n_total = 0
n_no_lf = 0
n_no_dm = 0

f = dataio.I3File(I3_PATH)
while f.more():
    frame = f.pop_frame()
    if frame.Stop != icetray.I3Frame.Physics:
        continue
    n_total += 1

    # ── Get LineFit track ─────────────────────────────────────────────────
    if TRACK_KEY not in frame:
        n_no_lf += 1
        continue
    lf  = frame[TRACK_KEY]
    # Track reference point + time (event-local ns — same as pulse times)
    lf_pos = np.array([lf.pos.x, lf.pos.y, lf.pos.z])
    lf_dir = np.array([lf.dir.x, lf.dir.y, lf.dir.z])   # travel direction
    lf_t0  = lf.time   # ns, event-local

    # ── DM-Ice hit ────────────────────────────────────────────────────────
    if "DMIce_detection_time" not in frame or "DMIce_detector" not in frame:
        n_no_dm += 1
        continue
    # DMIce_detection_time is in 0.1-ns DAQ ticks (same as utc_daq_time).
    # Subtract event start (ticks) then convert to ns.
    event_start_daq = frame["I3EventHeader"].start_time.utc_daq_time
    dm_t_ns = (frame["DMIce_detection_time"].value - event_start_daq) * 0.1  # event-local ns

    det_str  = str(frame["DMIce_detector"])           # e.g. 'I3String("det1")'
    det_key  = "det1" if "det1" in det_str else "det2"
    dm_pos   = DMICE_POS_IC[det_key]

    # ── Compute t_geo and d_perp ──────────────────────────────────────────
    t_geo, d_perp = t_geometric(lf_pos, lf_dir, lf_t0, dm_pos)
    delta_t = dm_t_ns - t_geo   # timing residual [ns]

    # ── Score ─────────────────────────────────────────────────────────────
    if d_perp < D_PERP_MAX:
        log_l_timing = (-0.5 * ((delta_t - mu_ns) / sigma_ns)**2
                        - math.log(sigma_ns * math.sqrt(2 * math.pi)))
        log_l_eff    = math.log(max(eps0, 1e-9))
        log_l        = log_l_timing + log_l_eff
    else:
        log_l = float("nan")

    year = frame["I3EventHeader"].start_time.utc_year

    records.append(dict(
        year     = year,
        det      = det_key,
        d_perp_m = d_perp,
        delta_t_ns = delta_t,
        log_l    = log_l,
        lf_speed = lf.speed,
    ))

f.close()
print(f"\nProcessed {n_total} P-frames  (no LF: {n_no_lf}, no DM hit: {n_no_dm})")
print(f"Records: {len(records)}")

import pandas as pd
df = pd.DataFrame(records)
print(df.describe())

# Save CSV
csv_path = os.path.join(OUT_DIR, "real_coincidence_scoring.csv")
df.to_csv(csv_path, index=False)
print(f"\nCSV saved: {csv_path}")

# ── Plots ─────────────────────────────────────────────────────────────────────

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import norm as sp_norm

fig, axes = plt.subplots(1, 3, figsize=(15, 5))

# 1. d_perp distribution
ax = axes[0]
ax.hist(df.d_perp_m.clip(0, 300), bins=60, color="steelblue", alpha=0.7)
ax.axvline(D_PERP_MAX, color="red", ls="--", label=f"d⊥ cut = {D_PERP_MAX:.0f} m")
ax.set_xlabel("d⊥ (m)")
ax.set_ylabel("Events")
ax.set_title(f"LineFit d⊥ to DM-Ice\n({len(df)} real coincidences)")
ax.legend()
ax.grid(True, alpha=0.3)

# 2. Δt distribution for through-detector events vs model
sub = df[df.d_perp_m < D_PERP_MAX].copy()
ax = axes[1]
if len(sub) > 0:
    dt_lo, dt_hi = np.percentile(sub.delta_t_ns, [1, 99])
    pad = max(200, 3 * sigma_ns)
    bins = np.linspace(min(dt_lo, mu_ns - 4*sigma_ns) - pad,
                       max(dt_hi, mu_ns + 4*sigma_ns) + pad, 40)
    ax.hist(sub.delta_t_ns, bins=bins, density=True, color="steelblue",
            alpha=0.7, label=f"Real data (n={len(sub)})")
    x_fine = np.linspace(bins[0], bins[-1], 400)
    ax.plot(x_fine, sp_norm.pdf(x_fine, mu_ns, sigma_ns),
            "r-", lw=2.5, label=f"Model: N({mu_ns:+.0f}, {sigma_ns:.0f}) ns")
    ax.axvline(0,      color="k",   ls="--", lw=0.8, alpha=0.5, label="t_geo=0")
    ax.axvline(mu_ns,  color="red", ls=":",  lw=1.2)
    ax.set_xlabel("Δt = t_DM − t_geo [ns]")
    ax.set_ylabel("Normalised events / bin")
    ax.set_title(f"Timing residual: real data vs model\n(d⊥ < {D_PERP_MAX:.0f} m, n={len(sub)})")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
else:
    ax.text(0.5, 0.5, "No events pass d⊥ cut", transform=ax.transAxes, ha="center")

# 3. log L distribution (through-detector events)
ax = axes[2]
valid = df[df.log_l.notna()]
ax.hist(valid.log_l, bins=30, color="darkorange", alpha=0.7, label=f"Through-det (n={len(valid)})")
ax.axvline(log_l_at_mu := (-0.5 * 0 - math.log(sigma_ns * math.sqrt(2*math.pi)) + math.log(eps0)),
           color="red", ls="--", label=f"Model peak: {log_l_at_mu:.1f}")
ax.set_xlabel("log L_DM-Ice")
ax.set_ylabel("Events")
ax.set_title("DM-Ice log likelihood\n(through-detector events only)")
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

fig.suptitle(
    f"DM-Ice model sanity check — real coincidence data (2011–2022)\n"
    f"{len(df)} events total | {len(sub)} through-detector (d⊥<{D_PERP_MAX:.0f}m) | "
    f"Model: μ={mu_ns:+.0f}ns σ={sigma_ns:.0f}ns",
    fontsize=10
)
plt.tight_layout()
plot_path = os.path.join(OUT_DIR, "real_coincidence_scoring.png")
fig.savefig(plot_path, dpi=150, bbox_inches="tight")
plt.close()
print(f"Plot saved: {plot_path}")

# ── Summary stats ─────────────────────────────────────────────────────────────
if len(sub) > 0:
    print(f"\nThrough-detector events (d⊥ < {D_PERP_MAX:.0f} m): n={len(sub)}")
    print(f"  Δt median: {sub.delta_t_ns.median():+.0f} ns  (model μ: {mu_ns:+.0f} ns)")
    print(f"  Δt std:    {sub.delta_t_ns.std():.0f} ns  (model σ: {sigma_ns:.0f} ns)")
    print(f"  Δt IQR:    {sub.delta_t_ns.quantile(0.25):.0f} to {sub.delta_t_ns.quantile(0.75):.0f} ns")
