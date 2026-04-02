---
name: DMice project context
description: IceCube DM-Ice research project — multi-machine workflow, simulation, reconstruction scripts, steamshovel visualization
type: project
---

DMice is a multi-machine workflow system for IceCube dark matter (DM-Ice) physics research.

**Repo location (local machine):** `~/Icecube26/dmice/` (NOT ~/dmice/)

**Machines:**
- LOCAL (`~/Icecube26/dmice`): development, editing, orchestration
- NPX (`npx-submitter.icecube.wisc.edu`): HTCondor job submission; scripts at `~/dmice_work/`
- COBALT (`cobalt-14`): interactive compute, long-running jobs via screen
- WARD: lab Linux desktop with GPU (~1s/event BLO); BLO not yet installed

**Key scripts in `~/Icecube26/dmice/`:**
- `simulate_muons.py` — Prometheus-based targeted downgoing muon sim aimed at DM-Ice detectors
- `sim_linefit_comparison.py` — IC-only vs pivot LineFit vs MC truth comparison (use THIS version, not dmice_results copy)
- `dmice_pivot_refit.py` — re-runs linefit with synthetic DM-Ice hit appended
- `run_phase1_pipeline.sh` — Phase 1 orchestrator
- `run_2020_2021_pipeline.sh` — coincidence pipeline for 2020/2021 (runs ON NPX)
- `BLO/batch_dm_ice_sim.py` — BLO-based upgoing muon sim (GPU, WARD)
- `BLO/WARD_sim_plan.md` — plan for BLO DM-Ice targeted sim on WARD

**Key scripts on NPX (`~/dmice_work/`):**
- `run_2020_2021_pipeline.sh` — currently running in screen `dmice_2020_2021`
- `find_dmice_coincidences.py`, `step3_submit.py`, `merge_output.py`, etc.

**DM-Ice detector positions (IceCube coordinates, meters):**
- det1: [31.25, -72.93, -511.05] → string 87, DOM 1
- det2: [-334.80, -424.50, -511.26] → string 88, DOM 1

**IceTray environment (local):** `~/.icevenv/i3/icetray/build/env-shell.sh`

**Steamshovel:**
- Launch: `~/launch_steamshovel.sh` (on cobalt)
- Custom artist: `~/.steamshovel/DMIcePivotLineFitArtist.py`
- WSL2/WSLg: use native Wayland (`QT_QPA_PLATFORM=wayland`, `DISPLAY=:0`) — do NOT force X11

**Condor submit files fix (2026-04-01):**
- All 46 submit files for 2020/2021 patched to write logs to `/scratch/bcharett/dmice_condor/step1_logs/` (was blocked on /data/user)

**Physics design:**
- DM-Ice NaI scintillator detects direct muon traversal (precise timing anchor = "pivot point")
- Muons downgoing (atmospheric), passing IC86 top-to-bottom then hitting DM-Ice at z≈-511m
- Simulation: muon start 1500m back along track from DM-Ice → enters IC86 from z≈+500–1000m
- Zenith 130–170° (Prometheus momentum convention); energy 1 TeV–1 PeV, gamma=2

**Direction convention:**
- Prometheus: zenith > 90° = downgoing (dz = cos(zenith) < 0)
- IceCube: opposite convention (anti-momentum)

**Coordinate system:**
- Prometheus/parquet depth z ≈ -2459 at DM-Ice; IceCube z = depth_z + 1948.07m → -511m

**Known bugs FIXED:**
- prometheus_to_i3.py: Z_OFFSET now applied to both DOM positions AND MC truth
- HTCondor submit: use `$(Process)` not `$((2000 + $(Process)))` in submit files
- AVX: added `requirements = (TARGET.has_avx == True)` to submit files

**Coincidence pipeline status (2026-04-01):**
- 2020/2021 pipeline running on NPX in screen `dmice_2020_2021`
- Log: `~/dmice_work/pipeline_2020_2021.log` on NPX
- Step1 condor jobs being submitted; will auto-proceed through steps 2–4
- When done: scp merged i3.zst to local `~/dmice_results/`, regenerate linefit pkls

**Project phases:**
- Phase 1 (CURRENT): Validate pivot fit on targeted simulation → need phase1_validation.png
- Phase 2: Characterize IC-only errors using real data (characterize_bias.py not yet written)
- Phase 3: Build correction model
- Phase 4: Cross-validate
- Phase 5: Apply to full dataset

**WARD BLO simulation (TODO — plan at `BLO/WARD_sim_plan.md`):**
- Goal: BLO-based sim of muons passing through DM-Ice (both up + downgoing), output NPZ + I3
- Status: BLO not yet installed on WARD; plan written, not started
- New script to write: `BLO/batch_dm_ice_targeted_sim.py`

**BLO (BlueLightOrchestra):**
- Location: `~/.icevenv/BLO/` (Julia package, PROPOSAL + PPC)
- GPU PPC: ~1s/event; CPU PPC: ~30s/event
- Python via `juliacall`

**Condor hardware (NPX):**
- GPU nodes available but PPC binary has no CUDA linkage — CPU-only (~10s/event)
- Python 3.9.25 on NPX
