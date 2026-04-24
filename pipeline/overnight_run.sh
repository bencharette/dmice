#!/bin/bash
# overnight_run.sh
#
# Run on Cobalt in a screen session:
#   screen -S overnight
#   bash ~/dmice/overnight_run.sh 2>&1 | tee ~/dmice_work/overnight.log
#
# Phases:
#   1. Run all recos on current fixed file (2012-2019, 8366 events)
#   2. Poll NPX until step3 condor jobs finish
#   3. Merge 2020-2021 coincidences and append to fixed file
#   4. Re-run all recos on full dataset (2012-2021)

set -euo pipefail

BASEDIR="/data/user/bcharett/dmice_coincidences_2011_2022"
FIXED_FILE="$BASEDIR/all_dmice_coincidences_2011_2022_fixed.i3.zst"
MERGED_2021="$BASEDIR/dmice_coincidences_2020_2021.i3.zst"
OUT_DIR="$HOME/dmice_work/output"
DMICE_DIR="$HOME/dmice"
ICETRAY="/cvmfs/icecube.opensciencegrid.org/py3-v4.3.0/RHEL_9_x86_64/metaprojects/icetray/v1.12.1/env-shell.sh"

mkdir -p "$OUT_DIR"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

# ── Phase 1: recos on fixed file (2012-2019) ──────────────────────────────────
log "=== PHASE 1: Recos on fixed file (2012-2019, 8366 events) ==="
log "  Output: $OUT_DIR/real_all_recos.csv"

$ICETRAY python3 "$DMICE_DIR/reco/run_all_recos_real.py" \
    > "$OUT_DIR/real_all_recos_phase1.log" 2>&1
log "Phase 1 complete."

# ── Phase 2: Poll NPX condor until step3 jobs finish ─────────────────────────
log "=== PHASE 2: Waiting for step3 condor jobs on NPX ==="
log "  Sleeping 3 minutes for jobs to appear in queue..."
sleep 180

while true; do
    TOTAL=$(ssh npx "condor_q bcharett 2>/dev/null | grep 'Total for bcharett:'" \
            | awk '{print $4}' || echo "1")
    log "  Condor jobs remaining for bcharett: $TOTAL"
    if [ "${TOTAL:-1}" -eq 0 ]; then
        break
    fi
    sleep 300
done
log "Step3 condor jobs finished."

# ── Phase 3: Merge 2020-2021 and append to fixed file ────────────────────────
log "=== PHASE 3: Merging 2020-2021 coincidences ==="

# Count how many step3 output files exist
N_COINC=$(find "$BASEDIR/step3_coincidences/2020" "$BASEDIR/step3_coincidences/2021" \
          -name "*_coinc.i3.zst" 2>/dev/null | wc -l)
log "  Step3 coincidence files found: $N_COINC"

if [ "$N_COINC" -eq 0 ]; then
    log "WARNING: No step3 coincidence files found for 2020-2021. Skipping merge."
    log "  Check step3 logs in $BASEDIR/step3_logs/2020 and /2021"
    log "  Skipping Phase 4."
    exit 0
fi

# Merge 2020-2021 into a standalone file
$ICETRAY python3 "$DMICE_DIR/analysis/merge_2020_2021.py" --output "$MERGED_2021" \
    > "$OUT_DIR/merge_2020_2021.log" 2>&1
log "  Merged file: $MERGED_2021"

# Append 2020-2021 frames into the fixed master file (in-place)
log "  Appending 2020-2021 to fixed master file..."
$ICETRAY python3 - <<'PYEOF'
import os, shutil
from icecube import icetray, dataio

BASEDIR    = "/data/user/bcharett/dmice_coincidences_2011_2022"
FIXED_FILE = os.path.join(BASEDIR, "all_dmice_coincidences_2011_2022_fixed.i3.zst")
NEW_FILE   = os.path.join(BASEDIR, "dmice_coincidences_2020_2021.i3.zst")
TMP        = FIXED_FILE + ".tmp"

print("Appending {} to {}".format(NEW_FILE, FIXED_FILE))

# Collect new frames (physics/DAQ only — master already has GCD)
new_frames = []
inf = dataio.I3File(NEW_FILE, "r")
for frame in inf:
    if frame.Stop in (icetray.I3Frame.DAQ, icetray.I3Frame.Physics):
        new_frames.append(frame)
inf.close()
print("  New frames: {}".format(len(new_frames)))

# Rewrite fixed file + new frames to tmp
outf = dataio.I3File(TMP, "w")
inf  = dataio.I3File(FIXED_FILE, "r")
n_orig = 0
for frame in inf:
    outf.push(frame)
    if frame.Stop == icetray.I3Frame.Physics:
        n_orig += 1
inf.close()
for frame in new_frames:
    outf.push(frame)
outf.close()

shutil.move(TMP, FIXED_FILE)
size_mb = os.path.getsize(FIXED_FILE) / 1e6
print("  Fixed file updated: {:.1f} MB ({} + {} physics frames)".format(
      size_mb, n_orig, len(new_frames)))
PYEOF

log "Phase 3 complete."

# ── Phase 4: Re-run recos on full dataset (2012-2021) ────────────────────────
log "=== PHASE 4: Recos on full dataset (2012-2021) ==="
log "  Output: $OUT_DIR/real_all_recos.csv (overwriting phase 1)"

$ICETRAY python3 "$DMICE_DIR/reco/run_all_recos_real.py" \
    > "$OUT_DIR/real_all_recos_phase4.log" 2>&1
log "Phase 4 complete."

log "=== ALL DONE ==="
log "  CSV:  $OUT_DIR/real_all_recos.csv"
log "  Logs: $OUT_DIR/real_all_recos_phase1.log"
log "        $OUT_DIR/merge_2020_2021.log"
log "        $OUT_DIR/real_all_recos_phase4.log"
