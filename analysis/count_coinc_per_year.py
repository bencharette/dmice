#!/usr/bin/env python3
"""
count_coinc_per_year.py

Quickly count coincident events per year from the master i3 file,
applying the era-aware fixes (InIceSplit/in_ice only, PoleMuonLinefit fallback).
No reconstructions run — just counting and d_perp computation.

Run on Cobalt:
  /cvmfs/.../env-shell.sh python3 -u ~/dmice/count_coinc_per_year.py
"""

import os, math
import numpy as np
from collections import defaultdict
from icecube import icetray, dataio, dataclasses

C_M_NS = 0.2998
DMICE_POS = {
    "det1": np.array([ 31.25,  -72.93, -511.05]),
    "det2": np.array([-334.80, -424.50, -511.26]),
}
MU_NS   = 280.0
D_CUTS  = [15.0, 50.0, 100.0]   # show counts at multiple thresholds
IC_STRINGS = set(range(1, 87))
PULSE_PRIORITY = ["OfflinePulses", "SRTInIcePulses", "ReextractedInIcePulses", "InIcePulses"]
MUON_STREAMS   = {'', 'in_ice', 'InIceSplit'}

IN_FILE = "/data/user/bcharett/dmice_coincidences_2011_2022/all_dmice_coincidences_2011_2022_fixed.i3"

def d_perp(lf_pos, lf_dir, dm_pos):
    r = np.asarray(dm_pos) - np.asarray(lf_pos)
    d = np.asarray(lf_dir) / np.linalg.norm(lf_dir)
    s = float(np.dot(r, d))
    return math.sqrt(max(0.0, float(np.dot(r, r)) - s**2))

counts      = defaultdict(int)   # year → total InIceSplit events
has_dmice   = defaultdict(int)   # year → has DM-Ice keys
has_lf      = defaultdict(int)   # year → has LineFit seed
close       = {d: defaultdict(int) for d in D_CUTS}  # year → d_perp < threshold
seen        = set()
geo_omgeo   = {}
n_skipped_stream = 0
n_skipped_dedup  = 0

f = dataio.I3File(IN_FILE)
n = 0
while f.more():
    frame = f.pop_frame()

    if frame.Stop == icetray.I3Frame.Geometry:
        geo = frame["I3Geometry"]
        geo_omgeo = {k: v for k, v in geo.omgeo.items()}
        continue

    if frame.Stop != icetray.I3Frame.Physics:
        continue

    hdr    = frame["I3EventHeader"]
    stream = getattr(hdr, 'sub_event_stream', '')

    # Skip non-muon streams
    if stream not in MUON_STREAMS:
        n_skipped_stream += 1
        continue

    uid = (hdr.run_id, hdr.event_id, stream)
    if uid in seen:
        n_skipped_dedup += 1
        continue
    seen.add(uid)

    year = hdr.start_time.utc_year

    # Pulse key check
    pulse_key = next((k for k in PULSE_PRIORITY if k in frame), None)
    if pulse_key is None:
        continue
    try:
        all_pulses = dataclasses.I3RecoPulseSeriesMap.from_frame(frame, pulse_key)
    except Exception:
        continue
    n_ic = sum(1 for omk, _ in all_pulses.items() if omk.string in IC_STRINGS)
    if n_ic < 4:
        continue

    counts[year] += 1
    n += 1
    if n % 500 == 0:
        print(f"  {n} events processed...", flush=True)

    # DM-Ice keys
    if "DMIce_detection_time" not in frame or "DMIce_detector" not in frame:
        continue
    has_dmice[year] += 1

    # LineFit seed
    lf_key = next((k for k in ("LineFit", "PoleMuonLinefit") if k in frame), None)
    if lf_key is None:
        continue
    has_lf[year] += 1

    lf = frame[lf_key]
    lf_pos = np.array([lf.pos.x, lf.pos.y, lf.pos.z])
    lf_dir = np.array([lf.dir.x, lf.dir.y, lf.dir.z])

    det_raw = str(frame["DMIce_detector"])
    det_key = "det1" if "det1" in det_raw else "det2"
    dp = d_perp(lf_pos, lf_dir, DMICE_POS[det_key])
    for d in D_CUTS:
        if dp < d:
            close[d][year] += 1

f.close()

print(f"\nSkipped: {n_skipped_stream} non-muon stream frames, {n_skipped_dedup} duplicates")
header = f"{'Year':<6} {'Total':>6} {'HasDMIce':>9} {'HasLF':>7}" + \
         "".join(f"  d<{int(d)}m" for d in D_CUTS)
print(f"\n{header}")
print("-" * len(header))
for y in sorted(counts):
    row = f"{y:<6} {counts[y]:>6} {has_dmice[y]:>9} {has_lf[y]:>7}"
    for d in D_CUTS:
        row += f"  {close[d][y]:>5}"
    print(row)
row = f"{'Total':<6} {sum(counts.values()):>6} {sum(has_dmice.values()):>9} {sum(has_lf.values()):>7}"
for d in D_CUTS:
    row += f"  {sum(close[d].values()):>5}"
print(row)
