#!/usr/bin/env python3
"""
Step 3 condor submission for 2020-2021 only.
Run on NPX after step2 (subrun finding) has completed on cobalt.

Usage:
    python3 step3_submit_2020_2021.py [--submit]
"""

import os
import sys
import glob

YEARS         = ["2020", "2021"]
USER          = os.environ.get("USER", "bcharett")
BASEDIR       = "/data/user/{}/dmice_coincidences_2011_2022".format(USER)
STEP1_DIR     = os.path.join(BASEDIR, "step1_muons")
STEP2_DIR     = os.path.join(BASEDIR, "step2_subruns")
STEP3_DIR     = os.path.join(BASEDIR, "step3_coincidences")
CONDOR_DIR    = "/scratch/{}/dmice_condor/step3_2020_2021".format(USER)
STEP3_LOGS    = os.path.join(BASEDIR, "step3_logs")
STEP3_WRAPPERS = os.path.join(BASEDIR, "step3_wrappers_2020_2021")

COINC_SCRIPT  = "/home/bcharett/dmice_work/find_dmice_coincidences.py"
CVMFS_BASE    = "/cvmfs/icecube.opensciencegrid.org/py3-v4.3.0"
ICETRAY_VER   = "v1.12.1"
TIME_WINDOW   = "-10,60"


def make_dirs(*paths):
    for p in paths:
        os.makedirs(p, exist_ok=True)


def create_wrapper(muon_dir, i3file, detector, out_file, wrapper_path):
    script = """#!/bin/bash
set -e
ARCH=$({cvmfs}/os_arch.sh)
exec {cvmfs}/${{ARCH}}/metaprojects/icetray/{ver}/env-shell.sh \\
    python3 {coinc} {muon_dir} {i3file} {out_file} \\
    --detector {det} --time-window {tw}
""".format(
        cvmfs=CVMFS_BASE, ver=ICETRAY_VER, coinc=COINC_SCRIPT,
        muon_dir=muon_dir, i3file=i3file, out_file=out_file,
        det=detector, tw=TIME_WINDOW,
    )
    with open(wrapper_path, "w") as f:
        f.write(script)
    os.chmod(wrapper_path, 0o755)


def main():
    submit = "--submit" in sys.argv
    all_submits = []
    total = 0
    skipped = 0

    subrun_files = []
    for year in YEARS:
        subrun_files += sorted(glob.glob(
            os.path.join(STEP2_DIR, year, "*", "*", "*_subrun.txt")))

    print("Found {} subrun files for 2020-2021".format(len(subrun_files)))

    batches = {}
    for subrun_file in subrun_files:
        rel   = subrun_file[len(STEP2_DIR)+1:]
        parts = rel.split(os.sep)
        if len(parts) < 4:
            continue
        year, month, detector = parts[0], parts[1], parts[2]

        with open(subrun_file) as f:
            content = f.read().strip()
        if not content:
            continue
        fields = content.split("\t")
        if len(fields) < 2:
            continue
        i3file = fields[1].strip()

        if not os.path.isfile(i3file):
            print("  WARNING: i3 not found: {}".format(i3file))
            continue

        stem     = os.path.basename(subrun_file).replace("_subrun.txt", "")
        out_file = os.path.join(STEP3_DIR, year, month, detector, stem + "_coinc.i3.zst")

        if os.path.isfile(out_file):
            skipped += 1
            continue

        wrapper_dir = os.path.join(STEP3_WRAPPERS, year, month, detector)
        make_dirs(wrapper_dir, os.path.join(STEP3_DIR, year, month, detector))
        wrapper = os.path.join(wrapper_dir, stem + ".sh")
        muon_dir = os.path.join(STEP1_DIR, year, month, detector)
        create_wrapper(muon_dir, i3file, detector, out_file, wrapper)

        key = "{}/{}/{}".format(year, month, detector)
        batches.setdefault(key, []).append((stem, wrapper))
        total += 1

    for key in sorted(batches):
        year, month, detector = key.split("/")
        log_dir    = os.path.join(CONDOR_DIR, year, month, detector)
        stdout_dir = os.path.join(STEP3_LOGS, year, month, detector)
        make_dirs(log_dir, stdout_dir)

        sub_path = os.path.join(log_dir, "submit.sub")
        with open(sub_path, "w") as f:
            f.write("Universe       = vanilla\n")
            f.write("GetEnv         = False\n")
            f.write("request_cpus   = 1\n")
            f.write("request_memory = 4GB\n")
            f.write("request_disk   = 2GB\n")
            f.write("Executable     = $(WRAPPER)\n")
            f.write("Output         = {}/$(STEM).out\n".format(stdout_dir))
            f.write("Error          = {}/$(STEM).err\n".format(stdout_dir))
            f.write("Log            = {}/condor.log\n".format(log_dir))
            f.write("\n")
            f.write("Queue STEM, WRAPPER from (\n")
            for stem, wrapper in batches[key]:
                f.write("  {}, {}\n".format(stem, wrapper))
            f.write(")\n")
        all_submits.append(sub_path)

    print("New jobs:    {}".format(total))
    print("Skipped:     {}".format(skipped))

    if not submit:
        print("\nDry run. Re-run with --submit to submit.")
        return

    print("\nSubmitting...")
    failed = 0
    for sub in all_submits:
        ret = os.system("condor_submit {}".format(sub))
        if ret != 0:
            failed += 1
    print("Done. {} failures.".format(failed))


if __name__ == "__main__":
    main()
