#!/bin/bash
# run_2020_2021_pipeline.sh
# Orchestrates the DM-Ice coincidence pipeline for 2020 and 2021.
# Run ON NPX in a screen session — condor commands run locally, SSHes to cobalt only for
# step2 (subrun finding) and step4 (IceTray merge).
#
#   screen -S dmice_2020_2021
#   bash ~/dmice/run_2020_2021_pipeline.sh 2>&1 | tee ~/dmice/logs/pipeline_2020_2021.log
#
# Steps:
#   1. Submit step1 condor jobs on NPX (vetoRootMaster.py per ROOT file)
#   2. Wait for step1 to finish
#   3. Run step2 (find IceCube subruns) on cobalt  [SSH to cobalt]
#   4. Submit step3 condor jobs on NPX (find_dmice_coincidences.py per subrun)
#   5. Wait for step3 to finish
#   6. Merge all coincidence files on cobalt         [SSH to cobalt]

set -euo pipefail

BASEDIR="/data/user/bcharett/dmice_coincidences_2011_2022"
CONDOR_STEP1="${BASEDIR}/condor/step1"
STEP1_MUONS="${BASEDIR}/step1_muons"
SCRIPT_DIR="/home/bcharett/dmice_work"
ICETRAY_ENV="/cvmfs/icecube.opensciencegrid.org/py3-v4.3.0/RHEL_9_x86_64/metaprojects/icetray/v1.12.1/env-shell.sh"
YEARS="2020 2021"
CLUSTER_IDS_FILE="${BASEDIR}/step1_cluster_ids_2020_2021.txt"
STEP3_CLUSTERS_FILE="${BASEDIR}/step3_cluster_ids_2020_2021.txt"

mkdir -p ~/dmice/logs

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

# ─────────────────────────────────────────────────────────────────────────────
# STEP 1: Submit step1 condor jobs on NPX
# ─────────────────────────────────────────────────────────────────────────────
log "=== STEP 1: Submitting step1 condor jobs on NPX ==="

> "${CLUSTER_IDS_FILE}"

for year in ${YEARS}; do
    for det in det1 det2; do
        for submit_file in $(find "${CONDOR_STEP1}/${year}" -name 'submit.sub' -path "*/${det}/*" | sort); do
            # Check how many output files already exist
            month=$(echo "${submit_file}" | grep -oP '\d{2}(?=/'"${det}"')' | tail -1)
            expected_count=$(grep -c '^Queue' "${submit_file}" 2>/dev/null || echo 0)
            actual_count=$(find "${STEP1_MUONS}/${year}/${month}/${det}" -name '*.txt' 2>/dev/null | wc -l)
            if [ "${actual_count}" -ge "${expected_count}" ] && [ "${expected_count}" -gt 0 ]; then
                log "  SKIP (already done): ${year}/${month}/${det} (${actual_count}/${expected_count} files)"
                continue
            fi

            log "  Submitting: ${year}/${month}/${det}"
            cluster_out=$(condor_submit "${submit_file}" 2>&1)
            cluster_id=$(echo "${cluster_out}" | grep -oP 'cluster \K[0-9]+' || true)
            if [ -n "${cluster_id}" ]; then
                echo "${cluster_id}" >> "${CLUSTER_IDS_FILE}"
                log "    -> cluster ${cluster_id}"
            else
                log "    WARNING: no cluster ID from output: ${cluster_out}"
            fi
        done
    done
done

CLUSTER_COUNT=$(wc -l < "${CLUSTER_IDS_FILE}" 2>/dev/null || echo 0)
log "Submitted ${CLUSTER_COUNT} clusters."

# ─────────────────────────────────────────────────────────────────────────────
# STEP 1 WAIT: Poll NPX until all step1 clusters finish
# ─────────────────────────────────────────────────────────────────────────────
log "=== Waiting for step1 condor jobs to finish ==="

while true; do
    cluster_ids=$(cat "${CLUSTER_IDS_FILE}" 2>/dev/null || true)
    if [ -z "${cluster_ids}" ]; then
        log "No clusters to wait for — continuing."
        break
    fi

    still_running=0
    for cid in ${cluster_ids}; do
        count=$(condor_q "${cid}" 2>/dev/null | grep -cP "^\s*${cid}\." || true)
        still_running=$((still_running + count))
    done

    if [ "${still_running}" -eq 0 ]; then
        log "All step1 condor jobs finished."
        break
    fi
    log "  ${still_running} step1 jobs still queued/running. Sleeping 5 min..."
    sleep 300
done

# Count step1 output
for year in ${YEARS}; do
    count=$(find "${STEP1_MUONS}/${year}" -name '*.txt' 2>/dev/null | wc -l)
    log "  Step1 output ${year}: ${count} muon txt files"
done

# ─────────────────────────────────────────────────────────────────────────────
# STEP 2: Find IceCube subruns (runs on cobalt)
# ─────────────────────────────────────────────────────────────────────────────
log "=== STEP 2: Finding IceCube subruns ==="
ssh cobalt-14 "bash '${SCRIPT_DIR}/step2_run.sh'" 2>&1 | while IFS= read -r line; do log "  ${line}"; done
log "Step 2 complete."

# ─────────────────────────────────────────────────────────────────────────────
# STEP 3: Submit step3 condor jobs on NPX
# ─────────────────────────────────────────────────────────────────────────────
log "=== STEP 3: Submitting step3 condor jobs on NPX ==="

step3_out=$(cd "${SCRIPT_DIR}" && python3 step3_submit.py --submit 2>&1)
log "${step3_out}"
echo "${step3_out}" | grep -oP 'cluster \K[0-9]+' > "${STEP3_CLUSTERS_FILE}" || true
STEP3_COUNT=$(wc -l < "${STEP3_CLUSTERS_FILE}" 2>/dev/null || echo 0)
log "Step3 submitted ${STEP3_COUNT} clusters."

# ─────────────────────────────────────────────────────────────────────────────
# STEP 3 WAIT: Poll until all step3 clusters finish
# ─────────────────────────────────────────────────────────────────────────────
log "=== Waiting for step3 condor jobs to finish ==="

while true; do
    cluster_ids=$(cat "${STEP3_CLUSTERS_FILE}" 2>/dev/null || true)
    if [ -z "${cluster_ids}" ]; then
        log "No step3 clusters to wait for — continuing."
        break
    fi

    still_running=0
    for cid in ${cluster_ids}; do
        [ -z "${cid}" ] && continue
        count=$(condor_q "${cid}" 2>/dev/null | grep -cP "^\s*${cid}\." || true)
        still_running=$((still_running + count))
    done

    if [ "${still_running}" -eq 0 ]; then
        log "All step3 condor jobs finished."
        break
    fi
    log "  ${still_running} step3 jobs still queued/running. Sleeping 5 min..."
    sleep 300
done

# ─────────────────────────────────────────────────────────────────────────────
# STEP 4: Merge on cobalt (inside IceTray env)
# ─────────────────────────────────────────────────────────────────────────────
log "=== STEP 4: Merging all coincidence files ==="
ssh cobalt-14 "'${ICETRAY_ENV}' python3 '${SCRIPT_DIR}/merge_output.py'" 2>&1 | while IFS= read -r line; do log "  ${line}"; done

log "=== PIPELINE COMPLETE ==="
log "Merged output: ${BASEDIR}/all_dmice_coincidences_2011_2022.i3.zst"
log "Copy locally:  scp cobalt-14:${BASEDIR}/all_dmice_coincidences_2011_2022.i3.zst ~/dmice_results/"
