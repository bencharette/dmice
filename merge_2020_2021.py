#!/usr/bin/env python3
"""
Merge 2020-2021 coincidence i3 files and optionally append to the master file.
Run on cobalt inside the IceTray environment:

    /cvmfs/icecube.opensciencegrid.org/py3-v4.3.0/RHEL_9_x86_64/metaprojects/icetray/v1.12.1/env-shell.sh \\
        python3 ~/dmice/merge_2020_2021.py

Options:
    --output PATH       Where to write merged 2020-2021 file
                        (default: $BASEDIR/dmice_coincidences_2020_2021.i3.zst)
    --append            After merging, also append to the master
                        all_dmice_coincidences_2011_2022.i3.zst file
"""

import os
import sys
import glob

from icecube import icetray, dataio

USER      = os.environ.get("USER", "bcharett")
BASEDIR   = "/data/user/{}/dmice_coincidences_2011_2022".format(USER)
STEP3_DIR = os.path.join(BASEDIR, "step3_coincidences")
MASTER    = os.path.join(BASEDIR, "all_dmice_coincidences_2011_2022.i3.zst")
DEFAULT_OUT = os.path.join(BASEDIR, "dmice_coincidences_2020_2021.i3.zst")


def merge_files(input_files, output_path):
    print("Merging {} files -> {}".format(len(input_files), output_path))
    outfile = dataio.I3File(output_path, "w")
    written_gcd = False
    n_physics = 0

    for i, inpath in enumerate(input_files):
        if (i + 1) % 100 == 1:
            print("  [{}/{}] {}".format(i+1, len(input_files), os.path.basename(inpath)))
        try:
            infile = dataio.I3File(inpath, "r")
        except Exception as e:
            print("  WARNING: could not open {}: {}".format(inpath, e))
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
        if n_physics > 0:
            written_gcd = True

    outfile.close()
    size_mb = os.path.getsize(output_path) / 1e6
    print("Done. {} physics frames, {:.1f} MB".format(n_physics, size_mb))
    return n_physics


def append_to_master(new_file, master_file):
    """Append frames from new_file to master_file (rewrite master in place)."""
    import tempfile, shutil

    tmp = master_file + ".tmp"
    print("\nAppending {} to master file...".format(new_file))
    print("  Master: {}".format(master_file))

    # Collect new physics/DAQ frames (skip GCD — master already has one)
    new_frames = []
    infile = dataio.I3File(new_file, "r")
    for frame in infile:
        stop = frame.Stop
        if stop in (icetray.I3Frame.DAQ, icetray.I3Frame.Physics):
            new_frames.append(frame)
    infile.close()
    print("  New frames to append: {}".format(len(new_frames)))

    # Rewrite master + new frames to tmp, then replace
    outfile = dataio.I3File(tmp, "w")
    infile  = dataio.I3File(master_file, "r")
    for frame in infile:
        outfile.push(frame)
    infile.close()
    for frame in new_frames:
        outfile.push(frame)
    outfile.close()

    shutil.move(tmp, master_file)
    size_mb = os.path.getsize(master_file) / 1e6
    print("  Master updated: {:.1f} MB".format(size_mb))


def main():
    output_file = DEFAULT_OUT
    do_append   = "--append" in sys.argv
    for i, arg in enumerate(sys.argv[1:], 1):
        if arg == "--output" and i < len(sys.argv):
            output_file = sys.argv[i + 1]

    # Find all 2020-2021 coincidence files
    coinc_files = []
    for year in ["2020", "2021"]:
        coinc_files += sorted(glob.glob(
            os.path.join(STEP3_DIR, year, "*", "*", "*_coinc.i3.zst")))
    coinc_files = [f for f in coinc_files if os.path.getsize(f) > 5000]

    if not coinc_files:
        print("ERROR: No coincidence files found for 2020-2021 in {}".format(STEP3_DIR))
        print("Has step3 finished? Check: condor_q bcharett (on NPX)")
        sys.exit(1)

    print("Found {} coincidence files for 2020-2021".format(len(coinc_files)))
    for year in ["2020", "2021"]:
        n = sum(1 for f in coinc_files if "/{}/".format(year) in f)
        print("  {}: {} files".format(year, n))

    n_events = merge_files(coinc_files, output_file)

    if do_append:
        if not os.path.isfile(MASTER):
            print("WARNING: Master file not found at {}.".format(MASTER))
            print("  Run merge_output.py first to build the master, then re-run with --append.")
        else:
            append_to_master(output_file, MASTER)
            print("\nMaster file now contains 2012-2021 data.")

    print("\nSummary:")
    print("  2020-2021 merged file: {}".format(output_file))
    print("  Events: {}".format(n_events))
    if not do_append:
        print("\n  To append to master file, re-run with --append:")
        print("    python3 ~/dmice/merge_2020_2021.py --append")


if __name__ == "__main__":
    main()
