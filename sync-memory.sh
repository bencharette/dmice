#!/bin/bash
# Sync local Claude memory to shared /data/user path on Cobalt/NPX.
# Run this after updating memory on local machine.

LOCAL_MEM="$HOME/.claude/projects/-home-ben/memory"
REMOTE="cobalt-14:/data/user/bcharett/dmice_claude_memory/"

echo "Syncing memory to $REMOTE ..."
scp "$LOCAL_MEM"/*.md "$REMOTE" && echo "Done."
