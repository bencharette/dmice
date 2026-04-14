#!/usr/bin/env python3
"""
run_pivot_mpefit.py

For each DM-Ice coincidence event:
  1. Compute Pivot LineFit (anchor time to DM-Ice transit, using LineFit dir as seed)
  2. Re-run MPEFit seeded from Pivot LineFit direction
  3. Compare angular difference between standard MPEFit and Pivot-seeded MPEFit

Output CSV columns:
    year, run, event, detector,
    lf_zenith_deg, lf_azimuth_deg,
    mpe_zenith_deg, mpe_azimuth_deg,
    pivot_lf_zenith_deg, pivot_lf_azimuth_deg,
    pivot_mpe_zenith_deg, pivot_mpe_azimuth_deg,
    mpe_vs_pivotmpe_ang_diff_deg,
    lf_vs_pivotlf_ang_diff_deg,
    n_doms, n_hits

Usage (on Cobalt with IceTray env):
    /cvmfs/.../env-shell.sh python3 run_pivot_mpefit.py
"""

import os, csv, math
import numpy as np

GCD_FILE   = "/cvmfs/icecube.opensciencegrid.org/data/GCD/GeoCalibDetectorStatus_2013.56429_V1.i3.gz"
IN_FILE    = "/data/user/bcharett/dmice_coincidences_2011_2022/all_dmice_coincidences_2011_2022.i3.zst"
OUT_CSV    = os.path.expanduser("~/dmice_work/output/comparison/pivot_mpefit_results.csv")
SPLINE_DIR = "/cvmfs/icecube.opensciencegrid.org/data/photon-tables/splines"

# DM-Ice positions [m] in IceCube coords
DMICE_POS = {
    "det1": np.array([ 31.25,  -72.93, -511.05]),
    "det2": np.array([-334.80, -424.50, -511.26]),
}
C_M_NS = 0.2998

os.makedirs(os.path.dirname(OUT_CSV), exist_ok=True)

from icecube import icetray, dataio, dataclasses, recclasses
from icecube import lilliput, gulliver, gulliver_modules
import icecube.lilliput.segments
from icecube.icetray import I3Tray

# ── Pulse key priority (same as run_truncated_energy.py) ─────────────────────
PULSE_PRIORITY  = ["OfflinePulses", "SRTInIcePulses", "ReextractedInIcePulses", "InIcePulses"]
UNIFIED_PULSES  = "UnifiedPulses"
PIVOT_LF_KEY    = "DMIce_PivotLineFit"
PIVOT_MPE_KEY   = "DMIce_PivotMPEFit"

# ── Pivot LineFit (pure Python) ───────────────────────────────────────────────
def _wm(vals, ws):
    W = sum(ws)
    return sum(v * w for v, w in zip(vals, ws)) / W if W else 0.0

def pivot_linefit(xs, ys, zs, ts, ws, dm_pos, seed_dir):
    """LineFit anchored to DM-Ice transit time estimated from seed_dir."""
    cx, cy, cz = _wm(xs, ws), _wm(ys, ws), _wm(zs, ws)
    tb = _wm(ts, ws)
    d_proj = ((dm_pos[0]-cx)*seed_dir[0] + (dm_pos[1]-cy)*seed_dir[1]
              + (dm_pos[2]-cz)*seed_dir[2])
    t_dm = tb + d_proj / C_M_NS
    dts  = [t - t_dm for t in ts]
    drxs = [x - dm_pos[0] for x in xs]
    drys = [y - dm_pos[1] for y in ys]
    drzs = [z - dm_pos[2] for z in zs]
    den  = sum(w * dt * dt for w, dt in zip(ws, dts))
    if not den:
        return None
    vx = sum(w * dt * dr for w, dt, dr in zip(ws, dts, drxs)) / den
    vy = sum(w * dt * dr for w, dt, dr in zip(ws, dts, drys)) / den
    vz = sum(w * dt * dr for w, dt, dr in zip(ws, dts, drzs)) / den
    spd = math.sqrt(vx*vx + vy*vy + vz*vz)
    return (vx/spd, vy/spd, vz/spd) if spd else None

def dir_to_zen_azi(dx, dy, dz):
    zen = math.degrees(math.acos(max(-1.0, min(1.0, dz))))
    azi = math.degrees(math.atan2(dy, dx)) % 360.0
    return zen, azi

def ang_diff_deg(d1, d2):
    dot = max(-1.0, min(1.0, d1[0]*d2[0] + d1[1]*d2[1] + d1[2]*d2[2]))
    return math.degrees(math.acos(abs(dot)))

# ── IceTray modules ───────────────────────────────────────────────────────────
rows = []

def unify_pulses(frame):
    for key in PULSE_PRIORITY:
        if key in frame:
            frame[UNIFIED_PULSES] = frame[key]
            return

def compute_pivot_lf(frame):
    if "LineFit" not in frame or UNIFIED_PULSES not in frame:
        return

    lf = frame["LineFit"]
    lf_dir = (lf.dir.x, lf.dir.y, lf.dir.z)

    det_str = frame["DMIce_detector"].value if "DMIce_detector" in frame else "det1"
    dm_pos  = DMICE_POS.get(det_str, DMICE_POS["det1"])

    try:
        pulses = dataclasses.I3RecoPulseSeriesMap.from_frame(frame, UNIFIED_PULSES)
    except Exception:
        return

    xs, ys, zs, ts, ws = [], [], [], [], []
    for omk, pulse_list in pulses:
        for p in pulse_list:
            xs.append(0.0); ys.append(0.0); zs.append(0.0)  # placeholder
            ts.append(p.time)
            ws.append(p.charge)

    # Get DOM positions from geometry
    try:
        geo = frame["I3Geometry"].omgeo
        xs.clear(); ys.clear(); zs.clear(); ts.clear(); ws.clear()
        for omk, pulse_list in pulses:
            if omk not in geo:
                continue
            pos = geo[omk].position
            for p in pulse_list:
                xs.append(pos.x); ys.append(pos.y); zs.append(pos.z)
                ts.append(p.time); ws.append(p.charge)
    except Exception:
        return

    if len(xs) < 4:
        return

    piv = pivot_linefit(xs, ys, zs, ts, ws, dm_pos, lf_dir)
    if piv is None:
        return

    p_particle = dataclasses.I3Particle()
    p_particle.dir = dataclasses.I3Direction(piv[0], piv[1], piv[2])
    p_particle.fit_status = dataclasses.I3Particle.FitStatus.OK
    frame[PIVOT_LF_KEY] = p_particle

def extract(frame):
    hdr  = frame["I3EventHeader"]
    year = str(hdr.start_time.utc_year)
    run  = hdr.run_id
    evt  = hdr.event_id
    det  = frame["DMIce_detector"].value if "DMIce_detector" in frame else "unknown"

    def get_dir(key):
        try:
            p = frame[key]
            return (p.dir.x, p.dir.y, p.dir.z)
        except Exception:
            return None

    lf_dir       = get_dir("LineFit")
    mpe_dir      = get_dir("MPEFit")
    pivot_lf_dir = get_dir(PIVOT_LF_KEY)
    pivot_mpe_dir= get_dir(PIVOT_MPE_KEY)

    def zen_azi(d):
        if d is None: return float("nan"), float("nan")
        return dir_to_zen_azi(*d)

    lf_zen,  lf_azi  = zen_azi(lf_dir)
    mpe_zen, mpe_azi = zen_azi(mpe_dir)
    plf_zen, plf_azi = zen_azi(pivot_lf_dir)
    pmp_zen, pmp_azi = zen_azi(pivot_mpe_dir)

    mpe_vs_piv = ang_diff_deg(mpe_dir, pivot_mpe_dir) if mpe_dir and pivot_mpe_dir else float("nan")
    lf_vs_plf  = ang_diff_deg(lf_dir, pivot_lf_dir)  if lf_dir  and pivot_lf_dir  else float("nan")

    n_doms = n_hits = 0
    try:
        pulses = dataclasses.I3RecoPulseSeriesMap.from_frame(frame, UNIFIED_PULSES)
        n_doms = len(pulses)
        n_hits = sum(len(v) for v in pulses.values())
    except Exception:
        pass

    rows.append(dict(
        year=year, run=run, event=evt, detector=det,
        lf_zenith_deg=lf_zen, lf_azimuth_deg=lf_azi,
        mpe_zenith_deg=mpe_zen, mpe_azimuth_deg=mpe_azi,
        pivot_lf_zenith_deg=plf_zen, pivot_lf_azimuth_deg=plf_azi,
        pivot_mpe_zenith_deg=pmp_zen, pivot_mpe_azimuth_deg=pmp_azi,
        mpe_vs_pivotmpe_ang_diff_deg=mpe_vs_piv,
        lf_vs_pivotlf_ang_diff_deg=lf_vs_plf,
        n_doms=n_doms, n_hits=n_hits,
    ))

# ── Build tray ────────────────────────────────────────────────────────────────
tray = I3Tray()

tray.AddModule("I3Reader", FilenameList=[GCD_FILE, IN_FILE])

tray.AddModule(unify_pulses,     Streams=[icetray.I3Frame.Physics])
tray.AddModule(compute_pivot_lf, Streams=[icetray.I3Frame.Physics])

# Re-run MPEFit seeded from Pivot LineFit (Pandel MPE likelihood)
tray.Add(icecube.lilliput.segments.I3SinglePandelFitter,
    fitname = PIVOT_MPE_KEY,
    domllh  = "MPE",
    pulses  = UNIFIED_PULSES,
    seeds   = [PIVOT_LF_KEY],
    If      = lambda f: PIVOT_LF_KEY in f,
)

tray.AddModule(extract, Streams=[icetray.I3Frame.Physics])

tray.Execute()
tray.Finish()

fieldnames = list(rows[0].keys()) if rows else []
with open(OUT_CSV, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fieldnames)
    w.writeheader()
    w.writerows(rows)

n_piv = sum(1 for r in rows if not math.isnan(r["pivot_mpe_zenith_deg"]))
print(f"Done: {len(rows)} events, {n_piv} with PivotMPEFit")
print(f"Saved: {OUT_CSV}")
