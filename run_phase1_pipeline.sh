#!/bin/bash
# run_phase1_pipeline.sh
# Phase 1 validation: parquet -> i3 -> sim_linefit_comparison for all 20 runs.
# Run inside a screen session on NPX.
#
# Usage: bash run_phase1_pipeline.sh

set -e
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Auto-detect local vs NPX environment
if [ -d "${HOME}/dmice_work" ]; then
    ENV_SHELL="/cvmfs/icecube.opensciencegrid.org/py3-v4.3.0/RHEL_9_x86_64/metaprojects/icetray/v1.12.1/env-shell.sh"
    SIM_BASE="/data/user/bcharett/dmice_sim_output"
else
    ENV_SHELL="${HOME}/.icevenv/i3/icetray/build/env-shell.sh"
    SIM_BASE="${HOME}/dmice_sim_output"
fi
OUTDIR="${SCRIPT_DIR}/phase1_output"
mkdir -p "${OUTDIR}"

echo "========================================================"
echo "Phase 1 pipeline: $(date)"
echo "Output dir: ${OUTDIR}"
echo "========================================================"

for RUN in $(seq 2000 2019); do
    PARQUET="${SIM_BASE}/run_0${RUN}/${RUN}_photons.parquet"
    I3FILE="${OUTDIR}/run_${RUN}_sim.i3"
    CSV="${OUTDIR}/run_${RUN}_results.csv"
    PLOT="${OUTDIR}/run_${RUN}_plot.png"

    echo ""
    echo "=== Run ${RUN} ==="

    if [ ! -f "${PARQUET}" ]; then
        echo "  ERROR: parquet not found: ${PARQUET}"
        continue
    fi

    # Step 1+2: parquet -> npz -> i3
    bash "${SCRIPT_DIR}/run_prometheus_to_i3.sh" "${PARQUET}" "${I3FILE}"

    # Step 3: linefit comparison (i3 file contains its own geometry frame)
    "${ENV_SHELL}" python "${SCRIPT_DIR}/sim_linefit_comparison.py" \
        -i "${I3FILE}" -g "${I3FILE}" \
        --output "${CSV}" \
        --plot "${PLOT}"

    echo "  Run ${RUN} done."
done

echo ""
echo "========================================================"
echo "All runs processed. Merging results..."
echo "========================================================"
python3 "${SCRIPT_DIR}/merge_phase1_results.py" \
    "${OUTDIR}" \
    "${SCRIPT_DIR}/phase1_validation.png"

echo ""
echo "======================================================== "
echo "Pipeline complete: $(date)"
echo "Plot: ${SCRIPT_DIR}/phase1_validation.png"
echo "========================================================"
