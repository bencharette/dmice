#!/usr/bin/env python3
"""
step3_rerun_all_years.py

Re-run step3 coincidence finding for years 2014-2019, outputting to a new
directory (step3_coincidences_v2/) so existing results aren't overwritten.

Usage:
    python3 step3_rerun_all_years.py [--years 2014,2015,2016,2017,2018,2019] [--submit]

Run on NPX after ensuring step1 and step2 are complete.
"""

import os, sys, glob

YEARS        = [int(y) for y in
                next((a.split("=")[1] for a in sys.argv if a.startswith("--years=")),
                     "2014,2015,2016,2017,2018,2019").split(",")]
SUBMIT       = "--submit" in sys.argv

USER         = os.environ.get("USER", "bcharett")
BASEDIR      = f"/data/user/{USER}/dmice_coincidences_2011_2022"
STEP1_DIR    = os.path.join(BASEDIR, "step1_muons")
STEP2_DIR    = os.path.join(BASEDIR, "step2_subruns")
STEP3_DIR    = os.path.join(BASEDIR, "step3_coincidences_v2")   # new output dir
CONDOR_DIR   = f"/scratch/{USER}/dmice_condor/step3_rerun"
LOG_DIR      = os.path.join(BASEDIR, "step3_logs_v2")
WRAPPER_DIR  = os.path.join(BASEDIR, "step3_wrappers_v2")

COINC_SCRIPT = f"/home/{USER}/dmice_work/find_dmice_coincidences.py"
CVMFS_BASE   = "/cvmfs/icecube.opensciencegrid.org/py3-v4.3.0"
ICETRAY_VER  = "v1.12.1"
TIME_WINDOW  = "-10,60"


def make_dirs(*paths):
    for p in paths:
        os.makedirs(p, exist_ok=True)


def create_wrapper(muon_dir, i3file, detector, out_file, wrapper_path):
    script = f"""#!/bin/bash
set -e
ARCH=$({CVMFS_BASE}/os_arch.sh)
exec {CVMFS_BASE}/${{ARCH}}/metaprojects/icetray/{ICETRAY_VER}/env-shell.sh \\
    python3 {COINC_SCRIPT} {muon_dir} {i3file} {out_file} \\
    --detector {detector} --time-window {TIME_WINDOW}
"""
    with open(wrapper_path, "w") as f:
        f.write(script)
    os.chmod(wrapper_path, 0o755)


def main():
    total = 0
    skipped_done = 0
    skipped_missing = 0
    batches = {}

    for year in YEARS:
        subrun_files = sorted(glob.glob(
            os.path.join(STEP2_DIR, str(year), "*", "*", "*_subrun.txt")))
        print(f"{year}: {len(subrun_files)} subrun files")

        for subrun_file in subrun_files:
            rel   = subrun_file[len(STEP2_DIR)+1:]
            parts = rel.split(os.sep)
            if len(parts) < 4:
                continue
            yr, month, detector = parts[0], parts[1], parts[2]

            with open(subrun_file) as f:
                content = f.read().strip()
            if not content:
                continue
            fields = content.split("\t")
            if len(fields) < 2:
                continue
            i3file = fields[1].strip()
            if not os.path.isfile(i3file):
                skipped_missing += 1
                continue

            stem     = os.path.basename(subrun_file).replace("_subrun.txt", "")
            out_file = os.path.join(STEP3_DIR, yr, month, detector, stem + "_coinc.i3.zst")

            if os.path.isfile(out_file):
                skipped_done += 1
                continue

            wrap_dir = os.path.join(WRAPPER_DIR, yr, month, detector)
            out_dir  = os.path.join(STEP3_DIR, yr, month, detector)
            make_dirs(wrap_dir, out_dir)

            muon_dir = os.path.join(STEP1_DIR, yr, month, detector)
            wrapper  = os.path.join(wrap_dir, stem + ".sh")
            create_wrapper(muon_dir, i3file, detector, out_file, wrapper)

            key = f"{yr}/{month}/{detector}"
            batches.setdefault(key, []).append((stem, wrapper))
            total += 1

    print(f"\nJobs to submit: {total}")
    print(f"Already done:   {skipped_done}")
    print(f"Missing i3:     {skipped_missing}")

    if total == 0:
        print("Nothing to do.")
        return

    # Write condor submit files
    for key in sorted(batches):
        yr, month, detector = key.split("/")
        log_sub  = os.path.join(CONDOR_DIR, yr, month, detector)
        log_out  = os.path.join(LOG_DIR,    yr, month, detector)
        make_dirs(log_sub, log_out)

        sub_path = os.path.join(log_sub, "submit.sub")
        with open(sub_path, "w") as f:
            f.write("Universe       = vanilla\n")
            f.write("GetEnv         = False\n")
            f.write("request_cpus   = 1\n")
            f.write("request_memory = 4GB\n")
            f.write("request_disk   = 2GB\n")
            f.write("Executable     = $(WRAPPER)\n")
            f.write(f"Output         = {log_out}/$(STEM).out\n")
            f.write(f"Error          = {log_out}/$(STEM).err\n")
            f.write(f"Log            = {log_sub}/condor.log\n\n")
            f.write("Queue STEM, WRAPPER from (\n")
            for stem, wrapper in batches[key]:
                f.write(f"  {stem}, {wrapper}\n")
            f.write(")\n")

        if SUBMIT:
            ret = os.system(f"condor_submit {sub_path}")
            if ret != 0:
                print(f"  FAILED: {sub_path}")

    if not SUBMIT:
        print("\nDry run. Re-run with --submit to submit jobs.")
    else:
        print("\nAll jobs submitted. Monitor with: condor_q")


if __name__ == "__main__":
    main()
