#!/bin/bash
# launch_offset_sims.sh
# Launches simulate_muons_offset.py for each d_perp offset in separate screen sessions on WARD.
# Run this script ON WARD:
#   bash ~/dmice/launch_offset_sims.sh

BLO_RESOURCE="$HOME/.icevenv/BLO/resources"
PPC_EXE="$BLO_RESOURCE/PPC_executables/PPC_CUDA/ppc"
GEO_FILE="$BLO_RESOURCE/geofiles/icecube_with_dmice.geo"
PPC_TABLES="$BLO_RESOURCE/PPC_tables/south_pole"
SCRIPT="$HOME/dmice/sim/simulate_muons_offset.py"
LOGDIR="$HOME/dmice_work/output/logs"
mkdir -p "$LOGDIR"

declare -A RUNS
RUNS[0]=1000      # on-axis high-stats
RUNS[50]=500      # 50 m offset
RUNS[100]=500     # 100 m offset
RUNS[200]=500     # 200 m offset

for OFFSET in 0 50 100 200; do
    N=${RUNS[$OFFSET]}
    SESSION="dmice_off${OFFSET}m"
    LOGFILE="$LOGDIR/offset_${OFFSET}m.log"

    CMD="BLO_PPC_EXE=$PPC_EXE BLO_GEO_FILE=$GEO_FILE BLO_PPC_TABLES=$PPC_TABLES python3 $SCRIPT --offset $OFFSET --n $N"
    nohup bash -c "$CMD" > "$LOGFILE" 2>&1 &
    PID=$!
    echo "Launched PID $PID  (offset=${OFFSET}m, n=${N})  → $LOGFILE"
done

echo ""
echo "Monitor with:"
echo "  screen -ls"
echo "  tail -f ~/dmice_work/output/logs/offset_50m.log"
echo "  tail -f ~/dmice_work/output/logs/offset_0m.log"
