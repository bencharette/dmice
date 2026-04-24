#!/usr/bin/env python3
"""
Resubmit step1 (vetoRootMaster) for 2020 and 2021 only.

Fixes the original failure: wrapper used `exec python` which didn't have ROOT
in its path on worker nodes. This version uses `exec python3` which correctly
picks up $SROOT/lib/root from the py3-v4.3.0 environment.

Run ON NPX (submit node):
    python3 ~/dmice/step1_resubmit_2020_2021.py [--submit]
"""

import os
import sys
import glob

YEARS      = [2020, 2021]
MONTHS     = ["{:02d}".format(m) for m in range(1, 13)]
DETECTORS  = ["det1", "det2"]

DMICE_BASE  = "/data/exp/DM-Ice/{year}/filtered/pole/data/tree/{month}/std_processing"
VETO_SCRIPT = "/home/bcharett/dmice_work/vetoRootMaster.py"
CVMFS_SETUP = "/cvmfs/icecube.opensciencegrid.org/py3-v4.3.0/setup.sh"

USER        = os.environ.get("USER", "bcharett")
BASEDIR     = "/data/user/{}/dmice_coincidences_2011_2022".format(USER)
STEP1_OUT   = os.path.join(BASEDIR, "step1_muons")
CONDOR_DIR  = "/scratch/{}/dmice_condor/step1_2020_2021".format(USER)
STEP1_LOGS     = os.path.join(BASEDIR, "step1_logs")
STEP1_WRAPPERS = os.path.join(BASEDIR, "step1_wrappers_2020_2021")  # separate dir to avoid collision


def make_dirs(*paths):
    for p in paths:
        os.makedirs(p, exist_ok=True)


def create_wrapper(root_file, month, year, out_file, detector, wrapper_path):
    script = """#!/bin/bash
set -e

# Load IceCube py3-v4.3.0 environment
eval $({cvmfs})

# Add ROOT Python bindings explicitly
export PYTHONPATH=$SROOT/lib/root:$PYTHONPATH
export LD_LIBRARY_PATH=$SROOT/lib/root:$LD_LIBRARY_PATH

# Run vetoRootMaster with python3
exec python3 {veto} {root_file} {month} {year} {out_file} {det}
""".format(
        cvmfs=CVMFS_SETUP,
        veto=VETO_SCRIPT,
        root_file=root_file,
        month=month,
        year=year,
        out_file=out_file,
        det=detector,
    )
    with open(wrapper_path, "w") as f:
        f.write(script)
    os.chmod(wrapper_path, 0o755)


def create_submit_file(year, month, detector, root_files):
    out_dir     = os.path.join(STEP1_OUT,      str(year), month, detector)
    log_dir     = os.path.join(CONDOR_DIR,     str(year), month, detector)
    stdout_dir  = os.path.join(STEP1_LOGS,     str(year), month, detector)
    wrapper_dir = os.path.join(STEP1_WRAPPERS, str(year), month, detector)
    make_dirs(out_dir, log_dir, stdout_dir, wrapper_dir)

    submit_path = os.path.join(log_dir, "submit.sub")
    wrappers = []

    for root_file in root_files:
        basename = os.path.basename(root_file)
        stem     = basename.replace(".root", "")
        out_file = os.path.join(out_dir, stem + ".txt")
        wrapper  = os.path.join(wrapper_dir, stem + ".sh")

        # Force regeneration — old wrappers used `python`, not `python3`
        create_wrapper(root_file, month, year, out_file, detector, wrapper)
        wrappers.append((stem, wrapper))

    if not wrappers:
        return None

    with open(submit_path, "w") as f:
        f.write("Universe       = vanilla\n")
        f.write("GetEnv         = False\n")
        f.write("request_cpus   = 1\n")
        f.write("request_memory = 2GB\n")
        f.write("request_disk   = 1GB\n")
        f.write("Executable     = $(WRAPPER)\n")
        f.write("Output         = {}/$(STEM).out\n".format(stdout_dir))
        f.write("Error          = {}/$(STEM).err\n".format(stdout_dir))
        f.write("Log            = {}/condor.log\n".format(log_dir))
        f.write("\n")
        f.write("Queue STEM, WRAPPER from (\n")
        for stem, wrapper in wrappers:
            f.write("  {}, {}\n".format(stem, wrapper))
        f.write(")\n")

    return submit_path, len(wrappers)


def main():
    submit = "--submit" in sys.argv
    submit_files = []
    total_jobs = 0

    for year in YEARS:
        for month in MONTHS:
            dmice_dir = DMICE_BASE.format(year=year, month=month)
            if not os.path.isdir(dmice_dir):
                continue
            root_files = sorted(glob.glob(os.path.join(dmice_dir, "*.root")))
            root_files = [f for f in root_files if "darknoise" not in f.lower()]
            if not root_files:
                continue
            for detector in DETECTORS:
                result = create_submit_file(year, month, detector, root_files)
                if result:
                    sub, n = result
                    submit_files.append(sub)
                    total_jobs += n
                    print("  {}/{} {}: {} files".format(year, month, detector, n))

    print("\nTotal jobs to submit: {}".format(total_jobs))

    if not submit:
        print("\nDry run — no jobs submitted. Re-run with --submit to submit.")
        return

    print("\nSubmitting to condor...")
    failed = 0
    for sub in submit_files:
        ret = os.system("condor_submit {}".format(sub))
        if ret != 0:
            print("  FAILED: {}".format(sub))
            failed += 1
    print("\nDone. {} / {} submit files failed.".format(failed, len(submit_files)))
    print("Monitor with: condor_q")


if __name__ == "__main__":
    main()
