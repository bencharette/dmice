#!/bin/bash
# pipeline_2020_2021_followup.sh
#
# Runs step2 + step3 + merge for 2020-2021 DM-Ice coincidences,
# AFTER step1 condor jobs (step1_resubmit_2020_2021.py) have finished.
#
# Run this ON COBALT in a screen session:
#   screen -S pipeline_2020_2021
#   bash ~/dmice/pipeline_2020_2021_followup.sh 2>&1 | tee ~/dmice_work/pipeline_2020_2021_followup.log
#
# Prerequisites:
#   - step1_resubmit_2020_2021.py --submit has been run on NPX and completed
#   - Verify: find /data/user/bcharett/dmice_coincidences_2011_2022/step1_muons/{2020,2021} -name "*.txt" | wc -l
#
# Steps:
#   1. Step 2 (cobalt): find IceCube subruns for each 2020-2021 DM-Ice muon
#   2. Step 3 (NPX via ssh): submit condor jobs to find coincident IceCube events
#   3. Wait for step3 condor jobs
#   4. Merge 2020-2021 coincidences into a new i3 file (cobalt)
#   5. Append to existing all_dmice_coincidences_2011_2022.i3.zst

set -euo pipefail

BASEDIR="/data/user/bcharett/dmice_coincidences_2011_2022"
STEP1_DIR="${BASEDIR}/step1_muons"
STEP2_DIR="${BASEDIR}/step2_subruns"
STEP3_DIR="${BASEDIR}/step3_coincidences"
SCRIPT_DIR="/home/bcharett/dmice_work"
DMICE_DIR="/home/bcharett/dmice"
ICETRAY_ENV="/cvmfs/icecube.opensciencegrid.org/py3-v4.3.0/RHEL_9_x86_64/metaprojects/icetray/v1.12.1/env-shell.sh"
YEARS="2020 2021"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

# ── Check step1 produced output ───────────────────────────────────────────────
log "=== Checking step1 output for 2020-2021 ==="
N_MUONS=$(find "${STEP1_DIR}/2020" "${STEP1_DIR}/2021" -name "*.txt" 2>/dev/null | wc -l)
log "  Muon txt files found: ${N_MUONS}"
if [ "${N_MUONS}" -eq 0 ]; then
    log "ERROR: No step1 output found for 2020-2021."
    log "  Wait for condor jobs to finish, then re-run."
    log "  Check: condor_q bcharett (on NPX)"
    exit 1
fi

# ── Step 2: Find IceCube subruns (cobalt, serial) ─────────────────────────────
log "=== STEP 2: Finding IceCube subruns for 2020-2021 ==="
PROCESSED=0; FOUND=0; NOT_FOUND=0
SUBRUN_SCRIPT="${SCRIPT_DIR}/subrunDouble_fixed.py"

for year in 2020 2021; do
    for muon_file in $(find "${STEP1_DIR}/${year}" -name "*.txt" | sort); do
        PROCESSED=$((PROCESSED + 1))
        rel="${muon_file#${STEP1_DIR}/}"
        month=$(echo "${rel}" | cut -d'/' -f2)
        detector=$(echo "${rel}" | cut -d'/' -f3)
        stem=$(basename "${muon_file}" .txt)

        out_dir="${STEP2_DIR}/${year}/${month}/${detector}"
        mkdir -p "${out_dir}"
        out_file="${out_dir}/${stem}_subrun.txt"

        if [ -f "${out_file}" ] && [ -s "${out_file}" ]; then
            FOUND=$((FOUND + 1))
            continue
        fi

        python3 "${SUBRUN_SCRIPT}" "${muon_file}" "${out_file}" "${year}" 2>/dev/null

        if [ -f "${out_file}" ] && [ -s "${out_file}" ]; then
            FOUND=$((FOUND + 1))
        else
            NOT_FOUND=$((NOT_FOUND + 1))
            echo "${muon_file}" >> "${STEP2_DIR}/unmatched_muons_2020_2021.txt"
        fi

        if [ $((PROCESSED % 200)) -eq 0 ]; then
            log "  Progress: ${PROCESSED} | Found: ${FOUND} | Not found: ${NOT_FOUND}"
        fi
    done
done

log "Step 2 complete — found: ${FOUND}, not matched: ${NOT_FOUND}"

N_SUBRUNS=$(find "${STEP2_DIR}/2020" "${STEP2_DIR}/2021" -name "*.txt" 2>/dev/null | wc -l)
if [ "${N_SUBRUNS}" -eq 0 ]; then
    log "ERROR: Step 2 found no subruns. Cannot continue."
    exit 1
fi

# ── Step 3: Submit condor jobs via NPX ────────────────────────────────────────
log "=== STEP 3: Submitting coincidence finder jobs (NPX) ==="

# Copy the targeted step3 submit script to NPX and run it
scp "${DMICE_DIR}/step3_submit_2020_2021.py" npx:~/dmice/ 2>/dev/null || true
CLUSTER_IDS=$(ssh npx "python3 ~/dmice/step3_submit_2020_2021.py --submit 2>/dev/null" | grep -oP 'cluster \K[0-9]+' || true)
log "  Submitted clusters: ${CLUSTER_IDS:-none printed, check condor_q}"

log "=== Waiting for step3 condor jobs to finish ==="
log "  Polling every 5 minutes..."
while true; do
    RUNNING=$(ssh npx "condor_q bcharett 2>/dev/null | grep -c 'bcharett'" || echo 0)
    log "  Jobs remaining: ${RUNNING}"
    if [ "${RUNNING}" -eq 0 ]; then
        break
    fi
    sleep 300
done
log "All step3 jobs finished."

# ── Step 4: Merge 2020-2021 coincidences ──────────────────────────────────────
log "=== STEP 4: Merging 2020-2021 coincidence files ==="

N_COINC=$(find "${STEP3_DIR}/2020" "${STEP3_DIR}/2021" -name "*_coinc.i3.zst" 2>/dev/null | wc -l)
log "  Coincidence files found: ${N_COINC}"

MERGED_2020_2021="${BASEDIR}/dmice_coincidences_2020_2021.i3.zst"
${ICETRAY_ENV} python3 "${DMICE_DIR}/merge_2020_2021.py" --output "${MERGED_2020_2021}"

log "  Merged file: ${MERGED_2020_2021}"
ls -lh "${MERGED_2020_2021}"

log "=== DONE ==="
log "Next: append 2020-2021 to the master coincidence file."
log "  See merge_2020_2021.py --append for instructions."
