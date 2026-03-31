#!/bin/bash

# DMice Tmux Workspace Launcher
# Creates a tmux session with Neovim and Claude Code in side-by-side panes

MACHINE="${1:-local}"
SESSION="dmice-$MACHINE"

# Kill existing session if running
tmux kill-session -t "$SESSION" 2>/dev/null

# Create new session with large window (240 cols x 60 rows)
tmux new-session -d -s "$SESSION" -x 240 -y 60

# Start Neovim in left pane
tmux send-keys -t "$SESSION" "nvim-dmice-$MACHINE" Enter

# Split window vertically (creates right pane)
tmux split-window -h -t "$SESSION"

# Start Claude Code in right pane
tmux send-keys -t "$SESSION" "claude-code ~/dmice" Enter

# Select left pane (Neovim)
tmux select-pane -t "$SESSION:0.0"

# Attach to session
tmux attach-session -t "$SESSION"
