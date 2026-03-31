# Claude Code Config for DMice

Everything needed to restore my exact Claude Code setup on a new machine.

## File layout

```
claude-config/
  memory/           → Claude persistent memory (project context, preferences, feedback)
  claude-viz/       → Live TUI for monitoring Claude tool calls
  settings.json     → Global Claude Code settings + hooks
  settings.local.json → Machine-specific permissions (edit paths before using)
```

## Setup steps

### 1. Memory files

Find (or create) the project memory directory for this repo. Claude Code stores
project memories at:

```
~/.claude/projects/<url-encoded-path>/memory/
```

For `/home/<you>/dmice` the path will be something like:
`~/.claude/projects/-home-<you>-dmice/memory/`

Copy all files from `claude-config/memory/` into that directory.

To onboard a fresh Claude instance immediately, paste the contents of
`memory/project_dmice.md` at the start of your conversation, prefixed with:

> "Here is context about the project I'm working on:"

### 2. claude-viz (optional — live tool-call monitor)

```bash
cp -r claude-config/claude-viz ~/claude-viz
pip install rich        # only dependency
python3 ~/claude-viz/viz.py   # run in a separate terminal while using Claude
```

### 3. Claude settings

The hooks in `settings.json` wire up claude-viz automatically.
Merge its `hooks` block into your `~/.claude/settings.json`.

`settings.local.json` contains machine-specific allowed commands (SSH paths,
IceTray env paths, etc.). Review and adapt paths before copying to
`~/.claude/settings.local.json`.

## Memory contents

| File | What it stores |
|------|----------------|
| `project_dmice.md` | Full physics/workflow context: machines, scripts, coordinate conventions, bugs fixed, current phase |
| `user_role.md` | User background (IceCube physicist, IceTray/HTCondor experience) |
| `feedback_comparison_script.md` | Always use `sim_linefit_comparison.py` from `/home/bench/dmice/` |
| `reference_onboarding.md` | How to paste context into a new Claude conversation |
| `MEMORY.md` | Index of all memory files (loaded automatically by Claude) |
