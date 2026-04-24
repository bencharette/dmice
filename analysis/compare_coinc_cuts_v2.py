#!/usr/bin/env python3
"""
compare_coinc_cuts_v2.py

Three discriminants compared side-by-side on the master coincidence file:

  A. MPEFit Gaussian cut (original):
       Δt_mpe = dm_t_ns − t_PCA(MPEFit → crystal)
       pass: |Δt_mpe − 280| < 3σ = 243 ns

  B. First-hit timing anchor (new):
       Δt_1st = dm_t_ns − (t_first_IC_hit + s_first_to_crystal / c)
       pass: |Δt_1st − 280| < 3σ = 243 ns
       (uses earliest IC hit as anchor; avoids MPEFit vertex-time bias)

  C. LF vs PivotLF angular difference (new):
       Δθ = angle(LineFit, PivotLineFit)
       pass: Δθ < threshold (scan 5°, 10°, 20°)

Run on cobalt:
  /cvmfs/icecube.opensciencegrid.org/py3-v4.3.0/RHEL_9_x86_64/metaprojects/icetray/v1.12.1/env-shell.sh \\
    python3 -u ~/dmice/compare_coinc_cuts_v2.py [--year 2012]

Output:
  ~/dmice_work/output/coinc_cuts_v2.csv
  ~/dmice_work/output/coinc_cuts_v2.png
"""

import os, sys, math, argparse
import numpy as np

# ── Parameters ────────────────────────────────────────────────────────────────

MU_NS      = 280.0
SIGMA_NS   = 81.0
N_SIGMA    = 3.0
C_M_NS     = 0.2998
N_ICE      = 1.3195

DMICE_POS_IC = {
    "det1": np.array([ 31.25,  -72.93, -511.05]),
    "det2": np.array([-334.80, -424.50, -511.26]),
}

GCD_FILE = ("/cvmfs/icecube.opensciencegrid.org/data/GCD/"
            "GeoCalibDetectorStatus_2013.56429_V1.i3.gz")
I3_FILE  = ("/data/user/bcharett/dmice_coincidences_2011_2022/"
            "all_dmice_coincidences_2011_2022_fixed.i3")
OUT_DIR  = os.path.expanduser("~/dmice_work/output")

PULSE_PRIORITY = [
    "SplitInIcePulses", "OnlineL2_CleanedMuonPulses",
    "OfflinePulses", "SRTInIcePulses",
    "ReextractedInIcePulses", "InIcePulses",
]
MPE_KEYS     = ["MPEFit", "PoleMuonLlhFit"]
LF_KEYS      = ["LineFit", "PoleMuonLinefit"]
IC_STRINGS   = set(range(1, 87))
MUON_STREAMS = {'', 'in_ice', 'InIceSplit'}

# ── Args ──────────────────────────────────────────────────────────────────────

parser = argparse.ArgumentParser()
parser.add_argument("--year", type=int, default=None)
parser.add_argument("--gcd",  default=GCD_FILE)
args = parser.parse_args()

# ── IceTray ───────────────────────────────────────────────────────────────────

from icecube import icetray, dataio, dataclasses

# ── Helpers ───────────────────────────────────────────────────────────────────

def get_track(frame, keys):
    for k in keys:
        if k in frame:
            p = frame[k]
            if hasattr(p,"fit_status") and p.fit_status == dataclasses.I3Particle.FitStatus.OK:
                return p, k
    return None, None

def d_hat(p):
    return np.array([p.dir.x, p.dir.y, p.dir.z])

def t_pca(track, dm_pos):
    """t at closest approach of track to dm_pos."""
    r = dm_pos - np.array([track.pos.x, track.pos.y, track.pos.z])
    s = float(np.dot(r, d_hat(track)))
    return track.time + s / C_M_NS

def ang_diff_deg(p1, p2):
    """Angular difference between two I3Particles in degrees."""
    dh1 = d_hat(p1); dh2 = d_hat(p2)
    dot = max(-1.0, min(1.0, float(np.dot(dh1, dh2))))
    return math.degrees(math.acos(dot))

def get_om_positions(frame):
    geo = frame["I3Geometry"].omgeo
    pos = {}
    for omk, omg in geo.items():
        if omk.string in IC_STRINGS:
            pos[(omk.string, omk.om)] = np.array([omg.position.x,
                                                    omg.position.y,
                                                    omg.position.z])
    return pos

def get_pulses(frame):
    for key in PULSE_PRIORITY:
        if key not in frame:
            continue
        try:
            pmap = dataclasses.I3RecoPulseSeriesMap.from_frame(frame, key)
        except Exception:
            try:
                pmap = frame[key].apply(frame)
            except Exception:
                pmap = frame[key]
        return pmap, key
    return None, None

def compute_firsthit_delta_t(pmap, om_pos, track, dm_pos, dm_t_ns):
    """
    Δt using the earliest IC hit as timing anchor.
    Returns (delta_t_ns, d_perp_first_m) or (nan, nan) if no hits.
    """
    best_t   = math.inf
    best_pos = None

    for omk, pulses in pmap.items():
        if omk.string not in IC_STRINGS:
            continue
        key_t = (omk.string, omk.om)
        if key_t not in om_pos:
            continue
        if not pulses:
            continue
        t0 = min(p.time for p in pulses)
        if t0 < best_t:
            best_t   = t0
            best_pos = om_pos[key_t]

    if best_pos is None or not math.isfinite(best_t):
        return np.nan, np.nan

    # Project first-hit DOM onto track direction
    dh   = d_hat(track)
    r    = best_pos - np.array([track.pos.x, track.pos.y, track.pos.z])
    s    = float(np.dot(r, dh))
    perp = r - s * dh
    d_perp_first = float(np.linalg.norm(perp))

    # Cherenkov geometric delay for first-hit DOM
    # (how late the first photon arrives relative to track passage)
    cos_ch  = 1.0 / N_ICE
    sin_ch  = math.sqrt(1.0 - cos_ch**2)
    t_geo_first = d_perp_first / (C_M_NS * sin_ch) if d_perp_first > 0.1 else 0.0

    # Time when muon was at the first-hit DOM position on the track
    t_muon_at_first_dom = best_t - t_geo_first

    # Along-track distance from first-hit DOM to crystal
    s_dom_to_crystal = float(np.dot(dm_pos - best_pos, dh))

    # Expected muon transit time at crystal
    t_at_crystal = t_muon_at_first_dom + s_dom_to_crystal / C_M_NS

    delta_t = dm_t_ns - t_at_crystal
    return delta_t, d_perp_first

def pivot_linefit(xs, ys, zs, ts, ws, dm_pos, dm_t_corrected, seed_dir):
    """Pivot LineFit anchored at DM-Ice crystal."""
    dts  = [t - dm_t_corrected for t in ts]
    drxs = [x - dm_pos[0] for x in xs]
    drys = [y - dm_pos[1] for y in ys]
    drzs = [z - dm_pos[2] for z in zs]

    den = sum(w*dt*dt for w, dt in zip(ws, dts))
    if den < 1e-10:
        return None

    vx = sum(w*dt*dr for w, dt, dr in zip(ws, dts, drxs)) / den
    vy = sum(w*dt*dr for w, dt, dr in zip(ws, dts, drys)) / den
    vz = sum(w*dt*dr for w, dt, dr in zip(ws, dts, drzs)) / den
    spd = math.sqrt(vx*vx + vy*vy + vz*vz)
    if spd < 1e-10:
        return None

    # Disambiguate direction vs seed
    if vx*seed_dir[0] + vy*seed_dir[1] + vz*seed_dir[2] < 0:
        vx, vy, vz = -vx, -vy, -vz

    return (vx/spd, vy/spd, vz/spd)

# ── Load geometry ─────────────────────────────────────────────────────────────

om_pos = {}
gcd_f  = dataio.I3File(args.gcd)
while gcd_f.more():
    fr = gcd_f.pop_frame()
    if fr.Stop == icetray.I3Frame.Geometry:
        om_pos = get_om_positions(fr)
        print(f"Geometry: {len(om_pos)} IC OMs")
        break
gcd_f.close()

# ── Scan ──────────────────────────────────────────────────────────────────────

print(f"Reading: {I3_FILE}")
records = []
seen    = set()
n_total = 0

f = dataio.I3File(I3_FILE)
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

    uid = (hdr.run_id, hdr.event_id, stream)
    if uid in seen:
        continue
    seen.add(uid)
    n_total += 1

    # Tracks
    mpe, _ = get_track(frame, MPE_KEYS)
    lf,  _ = get_track(frame, LF_KEYS)
    if lf is None:
        continue

    # DM-Ice timing
    if "DMIce_detection_time" not in frame:
        continue

    det_str = str(frame["DMIce_detector"]) if "DMIce_detector" in frame else "det1"
    det_key = "det1" if "det1" in det_str else "det2"
    dm_pos  = DMICE_POS_IC[det_key]

    event_start_daq = hdr.start_time.utc_daq_time
    dm_t_ns  = (frame["DMIce_detection_time"].value - event_start_daq) * 0.1
    dm_t_cor = dm_t_ns - MU_NS

    # ── Discriminant A: MPEFit Gaussian ──────────────────────────────────────
    dt_mpe = np.nan
    if mpe is not None:
        dt_mpe = dm_t_ns - t_pca(mpe, dm_pos)

    # ── Discriminant B: first-hit anchor ─────────────────────────────────────
    pmap, _ = get_pulses(frame)
    dt_1st  = np.nan
    if pmap is not None:
        dt_1st, _ = compute_firsthit_delta_t(pmap, om_pos, lf, dm_pos, dm_t_ns)

    # ── Discriminant C: LF vs PivotLF angle ──────────────────────────────────
    da_lf_piv = np.nan
    if pmap is not None:
        # Build PivotLineFit inline using IC pulses
        xs, ys, zs, ts, ws = [], [], [], [], []
        for omk, pulses in pmap.items():
            if omk.string not in IC_STRINGS:
                continue
            key_t = (omk.string, omk.om)
            if key_t not in om_pos or not pulses:
                continue
            x, y, z = om_pos[key_t]
            t0 = min(p.time for p in pulses)
            q  = sum(p.charge for p in pulses)
            xs.append(x); ys.append(y); zs.append(z)
            ts.append(t0); ws.append(q)

        if len(xs) >= 3:
            seed = [lf.dir.x, lf.dir.y, lf.dir.z]
            piv_dir = pivot_linefit(xs, ys, zs, ts, ws, dm_pos, dm_t_cor, seed)
            if piv_dir is not None:
                lf_d = np.array([lf.dir.x, lf.dir.y, lf.dir.z])
                dot  = max(-1.0, min(1.0, float(np.dot(lf_d, np.array(piv_dir)))))
                da_lf_piv = math.degrees(math.acos(dot))

    # ── Cuts ─────────────────────────────────────────────────────────────────
    pass_mpe = (not np.isnan(dt_mpe) and
                abs(dt_mpe - MU_NS) < N_SIGMA * SIGMA_NS)
    pass_1st = (not np.isnan(dt_1st) and
                abs(dt_1st - MU_NS) < N_SIGMA * SIGMA_NS)
    pass_ang5  = (not np.isnan(da_lf_piv) and da_lf_piv < 5.0)
    pass_ang10 = (not np.isnan(da_lf_piv) and da_lf_piv < 10.0)
    pass_ang20 = (not np.isnan(da_lf_piv) and da_lf_piv < 20.0)

    records.append(dict(
        year        = year,
        run_id      = hdr.run_id,
        event_id    = hdr.event_id,
        detector    = det_key,
        dm_t_ns     = dm_t_ns,
        lf_zen      = math.degrees(lf.dir.zenith),
        dt_mpe      = dt_mpe,
        dt_1st      = dt_1st,
        da_lf_piv   = da_lf_piv,
        has_mpe     = mpe is not None,
        pass_mpe    = pass_mpe,
        pass_1st    = pass_1st,
        pass_ang5   = pass_ang5,
        pass_ang10  = pass_ang10,
        pass_ang20  = pass_ang20,
        pass_both_mpe_ang10  = pass_mpe  and pass_ang10,
        pass_both_1st_ang10  = pass_1st  and pass_ang10,
    ))

f.close()
print(f"\nProcessed: {n_total} events  Records: {len(records)}")

# ── Summary ───────────────────────────────────────────────────────────────────

n = len(records)
def s(col): return sum(1 for r in records if r.get(col))

print(f"\n{'Discriminant':<35} {'N pass':>8} {'%':>6}")
print("-"*52)
for col, lab in [
    ('pass_mpe',            'MPEFit Gaussian (|Δt−280|<243 ns)'),
    ('pass_1st',            'First-hit anchor (|Δt−280|<243 ns)'),
    ('pass_ang5',           'LF vs PivotLF Δθ < 5°'),
    ('pass_ang10',          'LF vs PivotLF Δθ < 10°'),
    ('pass_ang20',          'LF vs PivotLF Δθ < 20°'),
    ('pass_both_mpe_ang10', 'MPEFit Gauss AND Δθ<10°'),
    ('pass_both_1st_ang10', 'First-hit AND Δθ<10°'),
]:
    np_ = s(col)
    print(f"  {lab:<33} {np_:>8}  {100*np_/n:>5.1f}%")

# Per-year
print("\nPer-year breakdown:")
from itertools import groupby
for yr in sorted(set(r['year'] for r in records)):
    yr_r = [r for r in records if r['year']==yr]
    print(f"  {yr}: N={len(yr_r):4d}  "
          f"mpe={s.__class__.__call__(s, 'pass_mpe') if False else sum(r['pass_mpe'] for r in yr_r):3d}  "
          f"1st={sum(r['pass_1st'] for r in yr_r):3d}  "
          f"Δθ<5={sum(r['pass_ang5'] for r in yr_r):3d}  "
          f"Δθ<10={sum(r['pass_ang10'] for r in yr_r):3d}")

# ── Save CSV ──────────────────────────────────────────────────────────────────

import csv as csvmod
out_csv = os.path.join(OUT_DIR, "coinc_cuts_v2.csv")
fields  = list(records[0].keys())
with open(out_csv, 'w', newline='') as fh:
    w = csvmod.DictWriter(fh, fieldnames=fields)
    w.writeheader(); w.writerows(records)
print(f"\nCSV: {out_csv}")

# ── Plots ─────────────────────────────────────────────────────────────────────

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

dt_mpe_all  = [r['dt_mpe']    for r in records if not np.isnan(r['dt_mpe'])]
dt_1st_all  = [r['dt_1st']    for r in records if not np.isnan(r['dt_1st'])]
da_all      = [r['da_lf_piv'] for r in records if not np.isnan(r['da_lf_piv'])]

fig, axes = plt.subplots(1, 3, figsize=(18, 5))

# Panel 1: Δt_mpe vs Δt_1st
ax = axes[0]
has_both = [r for r in records if not np.isnan(r['dt_mpe']) and not np.isnan(r['dt_1st'])]
xs_ = [r['dt_mpe'] for r in has_both]
ys_ = [r['dt_1st'] for r in has_both]
ax.scatter(xs_, ys_, s=3, alpha=0.4, color='steelblue')
ax.axvline(MU_NS - N_SIGMA*SIGMA_NS, color='orange', lw=1, ls='--')
ax.axvline(MU_NS + N_SIGMA*SIGMA_NS, color='orange', lw=1, ls='--', label='MPE 3σ window')
ax.axhline(MU_NS - N_SIGMA*SIGMA_NS, color='lime',   lw=1, ls='--')
ax.axhline(MU_NS + N_SIGMA*SIGMA_NS, color='lime',   lw=1, ls='--', label='1st-hit 3σ window')
ax.axline((0,0), slope=1, color='white', lw=0.8, ls=':', alpha=0.5, label='y=x')
ax.set_xlim(-500, 2000); ax.set_ylim(-500, 2000)
ax.set_xlabel('Δt_mpe [ns]'); ax.set_ylabel('Δt_1st-hit [ns]')
ax.set_title(f'MPE vs first-hit Δt\n(n={len(has_both)})')
ax.legend(fontsize=8); ax.grid(alpha=0.2)

# Panel 2: Δt_1st distribution
ax = axes[1]
win = (MU_NS - N_SIGMA*SIGMA_NS, MU_NS + N_SIGMA*SIGMA_NS)
bins = np.linspace(min(dt_1st_all + [-500]), min(max(dt_1st_all + [3000]), 5000), 80)
ax.hist(dt_1st_all, bins=bins, color='steelblue', alpha=0.8)
ax.axvspan(*win, alpha=0.15, color='lime', label=f'3σ window [{win[0]:.0f},{win[1]:.0f}]ns')
ax.axvline(MU_NS, color='red', lw=1.5, ls='--', label='μ=280 ns')
ax.set_xlabel('Δt_first-hit [ns]')
ax.set_ylabel('Events')
ax.set_title(f'First-hit anchor Δt distribution\n{sum(r["pass_1st"] for r in records)} pass 3σ')
ax.legend(fontsize=8); ax.grid(alpha=0.2)
ax.set_xlim(-500, 3000)

# Panel 3: LF vs PivotLF Δθ
ax = axes[2]
bins_a = np.linspace(0, 90, 45)
ax.hist(da_all, bins=bins_a, color='darkorange', alpha=0.85)
for thr, col, lab in [(5,'lime','5°'), (10,'yellow','10°'), (20,'red','20°')]:
    np_ = sum(1 for d in da_all if d < thr)
    ax.axvline(thr, color=col, lw=1.5, ls='--',
               label=f'<{lab}: {np_} ({100*np_/len(da_all):.0f}%)')
ax.set_xlabel('Δθ(LF, PivotLF) [°]')
ax.set_ylabel('Events')
ax.set_title(f'LF vs PivotLF angular difference\n(computed from IC pulse series)')
ax.legend(fontsize=8); ax.grid(alpha=0.2)

year_str = str(args.year) if args.year else "all years"
fig.suptitle(f'Three coincidence discriminants — {year_str} ({n} events)', fontsize=12)
plt.tight_layout()
out_png = os.path.join(OUT_DIR, "coinc_cuts_v2.png")
fig.savefig(out_png, dpi=150, bbox_inches='tight')
plt.close()
print(f"Plot: {out_png}")
print("Done.")
