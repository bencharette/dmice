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
- WARD: lab Linux desktop with GPU; BLO installed and working (CPU PPC only — no CUDA yet)

**Key scripts in `~/Icecube26/dmice/`:**
- `simulate_muons.py` — Prometheus-based targeted downgoing muon sim aimed at DM-Ice detectors
- `sim_linefit_comparison.py` — IC-only vs pivot LineFit vs MC truth comparison (use THIS version, not dmice_results copy)
- `dmice_pivot_refit.py` — re-runs linefit with synthetic DM-Ice hit appended
- `run_phase1_pipeline.sh` — Phase 1 orchestrator
- `run_2020_2021_pipeline.sh` — coincidence pipeline for 2020/2021 (runs ON NPX)
- `BLO/batch_dm_ice_sim.py` — BLO-based upgoing muon sim (GPU, WARD)
- `BLO/batch_dm_ice_targeted_sim.py` — BLO targeted sim: both directions, both detectors, NPZ output
- `BLO/blo_npz_to_i3.py` — converts NPZ output to IceTray I3 format
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
- All 46 submit files for 2020/2021 patched: Output/Error on /data/user/, Log on /scratch/bcharett/dmice_condor/step1_logs/

**Physics design:**
- DM-Ice NaI scintillator detects direct muon traversal (precise timing anchor = "pivot point")
- Muons downgoing (atmospheric), passing IC86 top-to-bottom then hitting DM-Ice at z≈-511m
- Simulation: muon start 1500m back along track from DM-Ice → enters IC86 from z≈+500–1000m
- Zenith 130–170° (Prometheus momentum convention); energy 1 TeV–1 PeV, gamma=2

**Direction convention:**
- Prometheus/BLO: zenith > 90° = downgoing (dz = cos(zenith) < 0)
- IceCube: opposite convention (anti-momentum)

**Coordinate system:**
- BLO uses depth coords: depth_z = IceCube_z - 1948.07m
- DM-Ice depth z ≈ -2459m → IceCube z = -511m
- IC86 spans IceCube z ≈ -500 to +500m

**Known bugs FIXED:**
- prometheus_to_i3.py: Z_OFFSET now applied to both DOM positions AND MC truth
- HTCondor submit: use `$(Process)` not `$((2000 + $(Process)))` in submit files
- AVX: added `requirements = (TARGET.has_avx == True)` to submit files

**Coincidence pipeline status (2026-04-01):**
- 2020/2021 pipeline running on NPX in screen `dmice_2020_2021`
- Log: `~/dmice_work/pipeline_2020_2021.log` on NPX
- Step1 condor jobs submitted (~31k jobs); will auto-proceed through steps 2–4
- When done: scp merged i3.zst to local `~/dmice_results/`, regenerate linefit pkls

**WARD BLO simulation status (2026-04-01):**
- BLO fully installed on WARD (Julia, juliacall, CPU PPC compiled; no CUDA yet)
- `batch_dm_ice_targeted_sim.py` written and working
- 1000-event simulation currently running on WARD
- NPZ→I3 converter: `BLO/blo_npz_to_i3.py` (run on machine with IceTray)

**Project phases:**
- Phase 1 (CURRENT): Validate pivot fit on targeted simulation → need phase1_validation.png
- Phase 2: Characterize IC-only errors using real data (characterize_bias.py not yet written)
- Phase 3: Build correction model
- Phase 4: Cross-validate
- Phase 5: Apply to full dataset

**BLO (BlueLightOrchestra):**
- Location: `~/.icevenv/BLO/` (Julia package, PROPOSAL + PPC)
- CPU PPC: ~30s/event; GPU PPC: ~1s/event (needs CUDA)
- Python via `juliacall`
- WARD: pip3 install needs `--break-system-packages` flag (Python 3.14, PEP 668)

**Condor hardware (NPX):**
- GPU nodes available but PPC binary has no CUDA linkage — CPU-only (~10s/event)
- Python 3.9.25 on NPX
