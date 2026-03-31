#!/bin/bash
# Condor wrapper for DMice targeted muon simulation.
# Arguments: PROCESS_ID NEVENTS  (run number = 2000 + PROCESS_ID)
set -e

RUN=$((2000 + $1))
NEVENTS=${2:-500}

echo "[wrapper] Starting run=${RUN}, nevents=${NEVENTS} on $(hostname)"
echo "[wrapper] Date: $(date)"

# Prometheus lives in the user's home dir on all execute nodes (via AFS)
export PYTHONPATH="/home/bcharett/prometheus:${PYTHONPATH}"

# h5py and pyarrow are available in the system python3
python3 /home/bcharett/dmice/simulate_muons.py \
    --run "${RUN}" \
    --nevents "${NEVENTS}" \
    --det both

echo "[wrapper] Done. $(date)"
