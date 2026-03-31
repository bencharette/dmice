#!/bin/bash
# run_prometheus_to_i3.sh
# Converts Prometheus parquet output to .i3 (parquet -> npz -> i3).
# Works on NPX: step 1 uses system python3, step 2 uses CVMFS IceTray.
#
# Usage: bash run_prometheus_to_i3.sh <input.parquet> <output.i3>

set -e

INPUT_PARQUET="${1}"
OUTPUT_I3="${2}"

if [[ -z "$INPUT_PARQUET" || -z "$OUTPUT_I3" ]]; then
    echo "Usage: bash run_prometheus_to_i3.sh <input.parquet> <output.i3>"
    exit 1
fi

# Auto-detect local vs NPX environment
if [ -d "${HOME}/dmice_work" ]; then
    # NPX: scripts live in ~/dmice_work, use CVMFS IceTray
    SCRIPTS_DIR="${HOME}/dmice_work"
    ENV_SHELL="/cvmfs/icecube.opensciencegrid.org/py3-v4.3.0/RHEL_9_x86_64/metaprojects/icetray/v1.12.1/env-shell.sh"
else
    # Local: scripts live in ~/.icevenv/i3/scripts, use local IceTray build
    SCRIPTS_DIR="${HOME}/.icevenv/i3/scripts"
    ENV_SHELL="${HOME}/.icevenv/i3/icetray/build/env-shell.sh"
fi

NPZ_TMP="${OUTPUT_I3%.i3}.npz"

echo "[1/2] parquet -> npz..."
python3 "${SCRIPTS_DIR}/parquet_to_npz.py" "${INPUT_PARQUET}" "${NPZ_TMP}"

echo "[2/2] npz -> i3 (IceTray env)..."
"${ENV_SHELL}" python3 "${SCRIPTS_DIR}/prometheus_to_i3.py" "${NPZ_TMP}" "${OUTPUT_I3}"

rm -f "${NPZ_TMP}"
echo "[DONE] ${OUTPUT_I3}"
