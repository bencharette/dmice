# DMice Commands Reference

Quick reference for all DMice commands, keybindings, and workflows.

---

## Launching Neovim Profiles

```bash
nvim-dmice-local    # Launch Neovim with local machine profile
nvim-dmice-npx      # Launch Neovim with NPX machine profile
nvim-dmice-cobalt   # Launch Neovim with Cobalt machine profile
```

---

## Vim Commands (inside Neovim)

All commands are available in any DMice Neovim profile.

### Memory & Inbox

| Command | Action |
|---------|--------|
| `:DMiceMemory` | Open memory file for current machine |
| `:DMiceInbox` | Open inbox file for current machine |
| `:DMiceSend <machine>` | Append timestamped message to inbox-<machine>.md |

Examples:
```vim
:DMiceMemory           " Open ~/dmice/memory/local.md (if on local machine)
:DMiceInbox            " Open ~/dmice/handoff/inbox-local.md
:DMiceSend npx         " Write message to inbox-npx.md
:DMiceSend cobalt      " Write message to inbox-cobalt.md
```

### Synchronization

| Command | Action |
|---------|--------|
| `:DMiceSyncPull` | Pull handoff/ from all remote machines |
| `:DMiceSyncPush` | Push handoff/ to all remote machines |
| `:DMiceSyncAll` | Pull → Push → Open inbox |

### Dashboard

| Command | Action |
|---------|--------|
| `:DMiceDashboard` | Open reference dashboard with all commands and paths |

---

## Keybindings (Vim normal mode)

| Key | Action |
|-----|--------|
| `<leader>dm` | Open memory file (`:DMiceMemory`) |
| `<leader>di` | Open inbox file (`:DMiceInbox`) |
| `<leader>da` | Sync all then open inbox (`:DMiceSyncAll`) |
| `<leader>dd` | Open dashboard (`:DMiceDashboard`) |

### Examples

```vim
<leader>dm          " Quickly access machine memory
<leader>di          " Check inbox for messages
<leader>da          " Pull updates, push changes, view inbox
<leader>dd          " See all available commands
```

---

## Common Workflows

### Sending a Message to Another Machine

```vim
:DMiceSend npx           " Open npx inbox with timestamp prompt
" Type your message and save
:w
```

### Syncing Before Work

```vim
<leader>da               " Pull from remotes, push local changes, show inbox
```

Or step by step:

```vim
:DMiceSyncPull           " Get latest from npx and cobalt-14
:DMiceInbox              " Check messages for you
:DMiceSendLocal          " Reply to any pending items
:DMiceSyncPush           " Send your changes back
```

### Checking Machine Status

```vim
<leader>dd               " View dashboard with all machine info
<leader>dm               " Read machine-specific memory notes
```

---

## File Locations

### Memory Files (per-machine, NOT synced)

```
~/dmice/memory/local.md     # Local machine knowledge
~/dmice/memory/npx.md       # NPX machine knowledge
~/dmice/memory/cobalt.md    # Cobalt machine knowledge
```

### Handoff Files (synced across machines)

```
~/dmice/handoff/inbox-local.md      # Messages for LOCAL machine
~/dmice/handoff/inbox-npx.md        # Messages for NPX machine
~/dmice/handoff/inbox-cobalt.md     # Messages for Cobalt machine
~/dmice/handoff/done.md             # Archive of completed items
```

### Neovim Profiles

```
~/.config/nvim-dmice-local/init.lua     # Local machine Neovim config
~/.config/nvim-dmice-npx/init.lua       # NPX machine Neovim config
~/.config/nvim-dmice-cobalt/init.lua    # Cobalt machine Neovim config
```

---

## Machine-Specific Context

### LOCAL Machine

**Role:** Development, editing, orchestration

**Quick commands:**
```bash
ssh npx                # Connect to NPX submit node
ssh cobalt-14          # Connect to Cobalt compute node
nvim-dmice-local       # Open local Neovim profile
```

**Primary tasks:**
- Edit scripts and configs
- Submit jobs to NPX via Condor
- Orchestrate cross-machine workflows
- Sync handoff files to NPX and Cobalt

### NPX Machine

**Role:** Condor submission and job monitoring

**Quick commands:**
```bash
ssh cobalt-14          # Connect to Cobalt
nvim-dmice-npx        # Open NPX Neovim profile
condor_submit job.sub  # Submit Condor job (requires interactive Kerberos)
condor_q              # Check job queue
```

**Primary tasks:**
- Submit HTCondor jobs
- Monitor job progress
- Retrieve output files

### COBALT Machine

**Role:** Interactive compute and debugging

**Quick commands:**
```bash
nvim-dmice-cobalt     # Open Cobalt Neovim profile
screen -S <name>      # Start long-running session
python script.py      # Run scripts interactively
steamshovel <file.i3> # Open IceCube event viewer (source IceTray env first)
```

**Primary tasks:**
- Run scripts interactively
- Debug pipelines
- Long-running computations via screen/tmux

---

## SSH Configuration (Already Set Up)

```
ssh npx                    # Connect via proxy to NPX
ssh cobalt-14              # Connect via proxy to Cobalt-14
ssh pub                    # Connect to pub bastion
```

All connections route through `pub.icecube.wisc.edu` via ProxyJump.

---

## Steamshovel (IceCube Event Viewer)

Runs on Cobalt only. Source the IceTray environment first, then launch:

```bash
source /cvmfs/icecube.opensciencegrid.org/py3-v4.3.0/RHEL_9_x86_64/metaprojects/icetray/v1.12.1/env-shell.sh
steamshovel <file.i3>
```

Note: use `screen` for the session — do NOT use `nohup` with env-shell.sh.

---

## Tips & Tricks

1. **Use `<leader>da` often** — syncs before every work session
2. **Check `<leader>dd` when stuck** — dashboard lists all available commands
3. **Keep memory files updated** — add insights for future reference
4. **Archive to done.md** — move completed items for record-keeping
5. **Use screen on Cobalt** — never use nohup with env-shell.sh
6. **One message per change** — use `:DMiceSend <machine>` for each update

---

## Quick Reference Card

```
LAUNCH PROFILES:
  nvim-dmice-local | nvim-dmice-npx | nvim-dmice-cobalt

VIM COMMANDS:
  :DMiceMemory          :DMiceInbox          :DMiceSend <machine>
  :DMiceSyncPull        :DMiceSyncPush       :DMiceSyncAll
  :DMiceDashboard

KEYBINDINGS:
  <leader>dm  (memory)   <leader>di  (inbox)   <leader>da  (sync+inbox)
  <leader>dd  (dashboard)

MACHINES:
  local (dev)            npx (submit)         cobalt-14 (compute)
```

