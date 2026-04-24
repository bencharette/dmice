#!/usr/bin/env python3
"""
submit_reco_condor.py

Submit run_sim_all_recos.py as parallel Condor jobs, one per chunk.

Usage (on NPX):
    python3 ~/dmice/submit_reco_condor.py --npz PATH --out-dir DIR [--det DETECTOR] [--n-chunks 10] [--submit]

Example:
    python3 ~/dmice/submit_reco_condor.py \\
        --npz ~/dmice_work/output/muons_binned_5bins_200pbin_det1_repacked.npz \\
        --out-dir ~/dmice_work/output/comparison/det1_chunks \\
        --det det1 --n-chunks 10 --submit
"""

import os, sys, argparse, glob

CVMFS     = "/cvmfs/icecube.opensciencegrid.org/py3-v4.3.0"
ITRAY_VER = "v1.12.1"
RECO      = os.path.expanduser("~/dmice/reco/run_sim_all_recos.py")
USER      = os.environ.get("USER", "bcharett")
SCRATCH   = f"/scratch/{USER}/dmice_condor/reco"   # Condor logs — must be /scratch
DATA_BASE = f"/data/user/{USER}/dmice_condor/reco"  # wrappers + output

parser = argparse.ArgumentParser()
parser.add_argument("--npz",      required=True)
parser.add_argument("--out-dir",  required=True)
parser.add_argument("--det",      default=None, help="Detector override e.g. det_center")
parser.add_argument("--n-chunks", type=int, default=10)
parser.add_argument("--submit",   action="store_true")
args = parser.parse_args()

npz     = os.path.expanduser(args.npz)
out_dir = os.path.expanduser(args.out_dir)
stem    = os.path.basename(npz).replace("_repacked.npz", "").replace(".npz", "")

os.makedirs(out_dir, exist_ok=True)
wrap_dir    = os.path.join(DATA_BASE, stem, "wrappers")  # wrappers on /data
stdout_dir  = os.path.join(DATA_BASE, stem, "logs")      # stdout/stderr on /data
condor_log  = os.path.join(SCRATCH,   stem)              # condor log on /scratch
os.makedirs(wrap_dir,   exist_ok=True)
os.makedirs(stdout_dir, exist_ok=True)
os.makedirs(condor_log, exist_ok=True)

wrappers = []
for chunk_id in range(args.n_chunks):
    csv_out  = os.path.join(out_dir, f"chunk_{chunk_id:03d}.csv")
    det_arg  = args.det if args.det else "none"
    wrap_path = os.path.join(wrap_dir, f"chunk_{chunk_id:03d}.sh")

    script = f"""#!/bin/bash
set -e
ARCH=$({CVMFS}/os_arch.sh)
exec {CVMFS}/${{ARCH}}/metaprojects/icetray/{ITRAY_VER}/env-shell.sh \\
    python3 -u {RECO} {npz} {csv_out} {det_arg} {chunk_id} {args.n_chunks}
"""
    with open(wrap_path, "w") as f:
        f.write(script)
    os.chmod(wrap_path, 0o755)
    wrappers.append((f"chunk_{chunk_id:03d}", wrap_path))

# Write submit file
sub_path = os.path.join(DATA_BASE, stem, "submit.sub")
with open(sub_path, "w") as f:
    f.write("Universe       = vanilla\n")
    f.write("GetEnv         = False\n")
    f.write("request_cpus   = 1\n")
    f.write("request_memory = 4GB\n")
    f.write("request_disk   = 1GB\n")
    f.write("Executable     = $(WRAPPER)\n")
    f.write(f"Output         = {stdout_dir}/$(STEM).out\n")
    f.write(f"Error          = {stdout_dir}/$(STEM).err\n")
    f.write(f"Log            = {condor_log}/condor.log\n\n")
    f.write("Queue STEM, WRAPPER from (\n")
    for stem_name, wrapper in wrappers:
        f.write(f"  {stem_name}, {wrapper}\n")
    f.write(")\n")

print(f"NPZ:     {npz}")
print(f"Chunks:  {args.n_chunks}")
print(f"Out dir: {out_dir}")
print(f"Submit:  {sub_path}")

if args.submit:
    ret = os.system(f"condor_submit {sub_path}")
    if ret == 0:
        print(f"\nSubmitted {args.n_chunks} jobs. Monitor: condor_q")
        print(f"Merge when done: python3 ~/dmice/merge_reco_chunks.py --out-dir {out_dir}")
    else:
        print("condor_submit failed")
else:
    print("\nDry run. Add --submit to submit jobs.")
