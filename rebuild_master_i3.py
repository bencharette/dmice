#!/usr/bin/env python3
"""
rebuild_master_i3.py

Rebuild the master coincidence i3 file from scratch by merging all step3
coincidence files for 2012-2021. Writes to a new output path to avoid
corrupting the existing file.

Run on Cobalt inside IceTray environment:
  /cvmfs/icecube.opensciencegrid.org/py3-v4.3.0/RHEL_9_x86_64/metaprojects/icetray/v1.12.1/env-shell.sh \
      python3 ~/dmice/rebuild_master_i3.py
"""

import os
import glob
from icecube import icetray, dataio

BASEDIR   = "/data/user/bcharett/dmice_coincidences_2011_2022"
STEP3_DIR = os.path.join(BASEDIR, "step3_coincidences")
OUT_FILE  = os.path.join(BASEDIR, "all_dmice_coincidences_2012_2021.i3.zst")

# Gather all coinc files across all years
coinc_files = sorted(glob.glob(os.path.join(STEP3_DIR, "*", "*", "*", "*_coinc.i3.zst")))
coinc_files = [f for f in coinc_files if os.path.getsize(f) > 5000]
print(f"Found {len(coinc_files)} coincidence files")

years = {}
for f in coinc_files:
    y = f.split(STEP3_DIR+"/")[1].split("/")[0]
    years[y] = years.get(y, 0) + 1
for y in sorted(years):
    print(f"  {y}: {years[y]} files")

print(f"\nOutput: {OUT_FILE}")

outfile = dataio.I3File(OUT_FILE, "w")
written_gcd = False
n_physics = 0

for idx, inpath in enumerate(coinc_files):
    if (idx + 1) % 500 == 0:
        print(f"  [{idx+1}/{len(coinc_files)}]  physics so far: {n_physics}")
    try:
        infile = dataio.I3File(inpath, "r")
    except Exception as e:
        print(f"  WARNING: could not open {os.path.basename(inpath)}: {e}")
        continue
    for frame in infile:
        stop = frame.Stop
        if stop in (icetray.I3Frame.Geometry,
                    icetray.I3Frame.Calibration,
                    icetray.I3Frame.DetectorStatus):
            if not written_gcd:
                outfile.push(frame)
            continue
        if stop == icetray.I3Frame.TrayInfo:
            continue
        if stop in (icetray.I3Frame.DAQ, icetray.I3Frame.Physics):
            outfile.push(frame)
            if stop == icetray.I3Frame.Physics:
                n_physics += 1
    infile.close()
    written_gcd = True

outfile.close()
size_mb = os.path.getsize(OUT_FILE) / 1e6
print(f"\nDone. {n_physics} physics frames, {size_mb:.1f} MB")
print(f"Output: {OUT_FILE}")
