---
name: DMice GitHub repo and BLO repo
description: GitHub URLs for the dmice project repo and the BLO simulation package
type: reference
---

**DMice repo:** https://github.com/bencharette/dmice
- SSH: `git@github.com:bencharette/dmice.git`
- Contains all simulation scripts, phase 1 pipeline, claude-config/, and BLO setup folder
- SSH auth confirmed working via `git@github.com` (bencharette key)

**BLO (BlueLightOrchestra):** https://github.com/kcarloni/BlueLightOrchestra.jl
- Julia package for Cherenkov simulation (PROPOSAL + PPC)
- Clone to `~/.icevenv/BLO/` on each machine

**Repo folder structure (as of 2026-03-31):**
- `BLO/` — WARD setup: `batch_dm_ice_sim.py`, `icecube_with_dmice.geo`, README with install steps
- `claude-config/` — Claude Code memory files, settings.json, claude-viz scripts
- Root — simulation/analysis scripts (simulate_muons.py, sim_linefit_comparison.py, etc.)
