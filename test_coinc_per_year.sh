#!/bin/bash
# Quick per-year coincidence test: run find_dmice_coincidences.py on one
# representative subrun from each year and compare to the existing step3 file.
#
# Run on Cobalt inside IceTray environment:
#   /cvmfs/icecube.opensciencegrid.org/py3-v4.3.0/RHEL_9_x86_64/metaprojects/icetray/v1.12.1/env-shell.sh bash ~/dmice/test_coinc_per_year.sh

BASEDIR=/data/user/bcharett/dmice_coincidences_2011_2022
COINC_SCRIPT=~/dmice_work/find_dmice_coincidences.py
OUTDIR=~/dmice_work/output/coinc_test
mkdir -p "$OUTDIR"

for YEAR in 2014 2015 2016 2017 2018 2019; do
    echo "=== $YEAR ==="

    # Pick the first det1 subrun file that points to an existing i3 file
    SUBRUN_FILE=""
    I3FILE=""
    for f in "$BASEDIR/step2_subruns/$YEAR"/**/det1/*.txt; do
        candidate=$(awk '{print $2}' "$f" 2>/dev/null)
        if [ -f "$candidate" ]; then
            SUBRUN_FILE="$f"
            I3FILE="$candidate"
            break
        fi
    done

    if [ -z "$I3FILE" ]; then
        echo "  No valid subrun found for $YEAR, skipping"
        continue
    fi

    # Derive paths
    STEM=$(basename "$SUBRUN_FILE" _subrun.txt)
    REL=${SUBRUN_FILE#$BASEDIR/step2_subruns/}   # e.g. 2017/06/det1/dmice_run...
    MONTH=$(echo "$REL" | cut -d/ -f2)
    MUON_DIR="$BASEDIR/step1_muons/$YEAR/$MONTH/det1"
    EXISTING="$BASEDIR/step3_coincidences/$YEAR/$MONTH/det1/${STEM}_coinc.i3.zst"
    OUTFILE="$OUTDIR/${YEAR}_${STEM}_test.i3.zst"

    echo "  Subrun: $STEM"
    echo "  IceCube file: $(basename $I3FILE)"
    echo "  Muon dir: $MUON_DIR"

    # Run coincidence finder
    python3 "$COINC_SCRIPT" "$MUON_DIR" "$I3FILE" "$OUTFILE" \
        --detector det1 --time-window -10,60

    # Count P-frames in new vs existing output
    NEW_COUNT=$(python3 -c "
from icecube import icetray, dataio
n=0
try:
    f=dataio.I3File('$OUTFILE')
    while f.more():
        fr=f.pop_frame()
        if fr.Stop==icetray.I3Frame.Physics: n+=1
    f.close()
except: pass
print(n)
" 2>/dev/null)

    OLD_COUNT=$(python3 -c "
from icecube import icetray, dataio
n=0
try:
    f=dataio.I3File('$EXISTING')
    while f.more():
        fr=f.pop_frame()
        if fr.Stop==icetray.I3Frame.Physics: n+=1
    f.close()
except: pass
print(n)
" 2>/dev/null)

    echo "  Result: new=$NEW_COUNT  existing=$OLD_COUNT"
    echo ""
done

echo "Done. Test outputs in $OUTDIR"
