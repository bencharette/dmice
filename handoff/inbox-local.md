# Inbox for LOCAL Machine

Messages from NPX and Cobalt for the LOCAL machine.

---
Starting muon simulations with Prometheus for 2018 data

---
## Session summary 2026-04-08

### Completed this session
- **Angular error plots** (`output/ang_err_combined.png`): both IC-only LineFit + DM-Ice Pivot side-by-side, 0–20° y-axis, 1000-event BLO data. Script: `plot_ang_err_vs_ndoms_energy.py` (fixed path typo + column rename).
- **blo_python.py**: Python port of BlueLightOrchestra.jl. Working on Cobalt (CPU PPC). See memory for full details.
- **PPC compiled on Cobalt** at `~/dmice_work/ppc_cpu/ppc` with NSTR=100 (patched from 94).
- **BLO resources synced** to cobalt: `~/dmice/BlueLightOrchestra.jl/resources/`

### COMPLETED 2026-04-09: WARD BLO runs
- 2-event test run ✓ (440 + 641 DOMs)
- 200-event binned downgoing sim ✓ (`muons_binned_200ev.npz`)
- Angular error plots, event displays, dist plots all generated
- See project_ward_blo_run.md and project_dmice.md for full details
