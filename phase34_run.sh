#!/bin/bash
# phase34_run.sh — merge 2020-2021 and re-run recos on full dataset
# Run on Cobalt:
#   screen -S phase34
#   bash ~/dmice/phase34_run.sh 2>&1 | tee ~/dmice_work/phase34.log

set -euo pipefail

BASEDIR="/data/user/bcharett/dmice_coincidences_2011_2022"
FIXED_FILE="$BASEDIR/all_dmice_coincidences_2011_2022_fixed.i3.zst"
MERGED_2021="$BASEDIR/dmice_coincidences_2020_2021.i3.zst"
OUT_DIR="$HOME/dmice_work/output"
DMICE_DIR="$HOME/dmice"
ICETRAY="/cvmfs/icecube.opensciencegrid.org/py3-v4.3.0/RHEL_9_x86_64/metaprojects/icetray/v1.12.1/env-shell.sh"

mkdir -p "$OUT_DIR"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

# ── Phase 3: Merge 2020-2021 coincidences ────────────────────────────────────
log "=== PHASE 3: Merging 2020-2021 coincidences ==="

N_COINC=$(find "$BASEDIR/step3_coincidences/2020" "$BASEDIR/step3_coincidences/2021" \
          -name "*_coinc.i3.zst" 2>/dev/null | wc -l)
log "  Step3 coincidence files: $N_COINC"

# Merge 2020-2021 into standalone file
$ICETRAY python3 "$DMICE_DIR/merge_2020_2021.py" --output "$MERGED_2021" \
    > "$OUT_DIR/merge_2020_2021.log" 2>&1
log "  Merged file written: $MERGED_2021"
ls -lh "$MERGED_2021"

# Append 2020-2021 frames to the fixed master file (in-place)
log "  Appending to fixed master file..."
$ICETRAY python3 - <<'PYEOF'
import os, shutil
from icecube import icetray, dataio

BASEDIR    = "/data/user/bcharett/dmice_coincidences_2011_2022"
FIXED_FILE = os.path.join(BASEDIR, "all_dmice_coincidences_2011_2022_fixed.i3.zst")
NEW_FILE   = os.path.join(BASEDIR, "dmice_coincidences_2020_2021.i3.zst")
TMP        = FIXED_FILE + ".tmp"

print("Appending {} to {}".format(NEW_FILE, FIXED_FILE))

new_frames = []
inf = dataio.I3File(NEW_FILE, "r")
for frame in inf:
    if frame.Stop in (icetray.I3Frame.DAQ, icetray.I3Frame.Physics):
        new_frames.append(frame)
inf.close()
print("  New frames to append: {}".format(len(new_frames)))

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
log "  Output: $OUT_DIR/real_all_recos.csv"

$ICETRAY python3 "$DMICE_DIR/run_all_recos_real.py" \
    > "$OUT_DIR/real_all_recos_phase4.log" 2>&1

log "Phase 4 complete."
log "=== ALL DONE ==="
log "  CSV:      $OUT_DIR/real_all_recos.csv"
log "  Phase4 log: $OUT_DIR/real_all_recos_phase4.log"
