---
name: DMice project context
description: IceCube DM-Ice research project — multi-machine workflow, simulation, reconstruction scripts, steamshovel visualization
type: project
---

DMice is a multi-machine workflow system for IceCube dark matter (DM-Ice) physics research.

**Machines:**
- LOCAL (`/home/bench/dmice`): development, editing, orchestration
- NPX (`npx-submitter.icecube.wisc.edu`): HTCondor job submission
- COBALT (`cobalt-14`): interactive compute, long-running jobs via screen

**Key scripts in `/home/bench/dmice/`:**
- `simulate_muons.py` — TARGETED downgoing muon simulation via Prometheus. Generates a custom LI injection .h5 file placing muons AT the DM-Ice detector position (back-projected 1500m along track so they enter IC86 from the top), then runs Prometheus PPC. Accepts `--run`, `--nevents`, `--det` CLI args.
- `condor_sim/submit_dmice_sim.sub` — HTCondor submit file for NPX (20 jobs × 500 events, run numbers 2000–2019)
- `condor_sim/run_sim.sh` — Condor wrapper script
- `sim_linefit_comparison.py` — Compares IC-only LineFit vs IC+DM-Ice pivot LineFit vs MC truth; runs inside IceTray env on NPX/Cobalt
- `dmice_pivot_refit.py` — Re-runs linefit.simple with synthetic DM-Ice hit appended to IC pulses; IceTray pipeline

**Key scripts on NPX (`~/dmice_work/`):**
- `parquet_to_npz.py` — converts Prometheus parquet output to .npz (no IceTray needed)
- `prometheus_to_i3.py` — converts .npz to .i3 (needs IceTray; use local env-shell)
- `find_dmice_coincidences.py` — finds IceCube events coincident with DM-Ice muon detections
- `vetoRootMaster.py` — step1: reads DM-Ice ROOT files, extracts muon detections to TXT
- `subrunDouble_fixed.py` — step2: finds IceCube Level2 subrun file for each DM-Ice muon TXT
- `step3_submit.py` — step3: generates + submits condor jobs for coincidence finding
- `merge_output.py` / `merge_coincidences.py` — merges step3 output i3.zst files

**DM-Ice detector positions (IceCube coordinates, meters):**
- det1: [31.25, -72.93, -511.05]  → IceCube string 87, DOM 1 in icecube_with_dmice.geo
- det2: [-334.80, -424.50, -511.26] → IceCube string 88, DOM 1
- In Prometheus depth coordinates: det1 z = -2459.12m, det2 z = -2459.33m

**IceTray environment (local build):**
`/home/bench/.icevenv/i3/icetray/build/env-shell.sh`

**Steamshovel setup:**
- Launch script: `~/launch_steamshovel.sh`
- Custom artist: `~/.steamshovel/DMIcePivotLineFitArtist.py` — loaded via `~/.steamshovel/startup.py`
- Visualizes: DM-Ice pivot LineFit (orange), PoleMuonLinefit (blue), MC truth (green), DM-Ice detector spheres (violet)
- WSL2/WSLg: use `DISPLAY=:0`, `WAYLAND_DISPLAY=wayland-0`, `XDG_RUNTIME_DIR=/run/user/1000`, `QT_QPA_PLATFORM=wayland`. Do NOT use `ip route`-based DISPLAY, `LIBGL_ALWAYS_INDIRECT=1`, or `unset WAYLAND_DISPLAY` — forcing X11/Xwayland breaks cursor forwarding. Native Wayland mode fixes the invisible cursor issue.

**Physics design (DM-Ice targeting):**
- DM-Ice NaI scintillator detects DIRECT muon traversal only (not Cherenkov photons)
- This gives a precise timing anchor (no Cherenkov propagation uncertainty) = "pivot point"
- Muons are DOWNGOING (atmospheric), passing through IC86 top-to-bottom, then hitting DM-Ice
- Simulation places muon start 1500m back along track from DM-Ice → enters IC86 from z≈+500–1000m
- Zenith range: 130–170° (Prometheus momentum convention = downgoing); |dz| > 0.5
- Energy: 1 TeV – 1 PeV, gamma=2
- Condor: 20 jobs × 100 events = 2,000 events total, run numbers 2000–2019
- Output: `/data/user/bcharett/dmice_sim_output/run_NNNNN/NNNNN_photons.parquet`

**Known bugs FIXED:**
- `prometheus_to_i3.py`: MC truth z was in depth coords, not IceCube coords → now applies Z_OFFSET (+1948.07m) to MC truth position, matching DOM positions and DMICE_POS in analysis scripts
- Previous simulation (run 1337) used wrong (lower) energy range and was volume injection, not targeted

**Direction convention:**
- Prometheus stores directions as [zenith, azimuth] in MOMENTUM direction convention
- zenith > 90° = downgoing momentum (dz = cos(zenith) < 0)
- Our simulated muons: zenith 130–170° → dz = -0.64 to -0.98 (downgoing) ✓
- LI output file also uses Prometheus momentum convention (NOT IceCube anti-momentum convention)
- Previous simulate_muons.py comments claiming "IceCube convention" were WRONG

**Coordinate system:**
- Prometheus / parquet `sensor_pos_z` and MC truth z: depth coordinates (z ≈ -2459 at DM-Ice)
- IceCube-centred z = depth_z + Z_OFFSET (Z_OFFSET = 1948.07m)
- DM-Ice at depth z ≈ -2459 → IceCube z = -511m (below IC86 main array)
- IC86 spans IceCube z ≈ -500 to +500m

**Condor hardware (NPX):**
- GPU nodes: SuperMicro 4027GR-TR with 8×GTX980 or 4×GTX1080 AVAILABLE
- BUT: PPC binary (`~/prometheus/resources/PPC_executables/PPC_CUDA/ppc`) has NO CUDA linkage — CPU-only
- Jobs run in CPU mode (`device=-1`), ~10s/event for PROPOSAL propagation; PPC is slower
- Python version on NPX: 3.9.25

**IceTray environment (local build):**
`/home/bench/.icevenv/i3/icetray/build/env-shell.sh`
Run scripts with: `/home/bench/.icevenv/i3/icetray/build/env-shell.sh python3 script.py`

**prometheus_to_i3.py:** `/home/bench/.icevenv/i3/scripts/prometheus_to_i3.py`
- Takes .npz (from parquet_to_npz.py), builds I3MCPESeriesMap and InIcePulses
- Uses `icecube_with_dmice.geo` as default geo, Z_OFFSET=1948.07m
- Z_OFFSET now applied to BOTH DOM positions AND MC truth vertex position (bug fixed 2026-03-23)
- Must be run inside IceTray env

**HTCondor submit file bug (FIXED 2026-03-23):**
- `$((2000 + $(Process)))` is bash arithmetic, NOT valid HTCondor syntax — passes literal string to script
- Fix: pass `$(Process)` in submit file, compute `RUN=$((2000 + $1))` in run_sim.sh
- First batch (cluster 10323162) failed with "invalid int value: '$((2000'" — all 20 jobs errored
- Second batch (cluster 10323163) submitted with fix — jobs running correctly

**AVX issue on Condor nodes (FIXED 2026-03-23):**
- jaxlib (used by Prometheus) is built with AVX instructions — crashes on older CPU nodes
- Error: `RuntimeError: This version of jaxlib was built using AVX instructions...`
- 7/20 jobs in cluster 10323163 landed on non-AVX nodes and failed immediately
- Fix: added `requirements = (TARGET.has_avx == True)` to submit files
- Resubmit cluster 10327036 submitted for the 7 failed jobs (processes 1,5,6,8,10,11,14)
- Both submit files now have the AVX requirement baked in

**Test run results (run 9998, 2 events, det1):**
- 708 and 208 photon hits (avg 458) — vs 13.7 avg on old untargeted run
- Confirms targeted injection pipeline works end-to-end

**Coincidence data — key paths on cobalt:**
- Proper pipeline dir: `/data/user/bcharett/dmice_coincidences_2011_2022/`
  - Merged file (2012–2019): `all_dmice_coincidences_2011_2022.i3.zst` (will be updated to 2012–2021 when pipeline finishes)
  - Step3 per-subrun outputs: `step3_coincidences/YEAR/MONTH/DET/*_coinc.i3.zst`
  - ROOT source files: `/data/exp/DM-Ice/YEAR/filtered/pole/data/tree/MONTH/std_processing/*.root`
- Old ad-hoc outputs (ignore, only 2017 data): `~/dmice_work/full_coincidence_outputs/`, `~/dmice_work/parallel_coincidence_outputs/`
- Local copy: `/home/bench/dmice_results/all_dmice_coincidences_2011_2022.i3.zst` (25MB, md5: 31fba8a046b7566ccf810f4374e38e72) — currently 2012–2019 only, stale after pipeline finishes
- ~14,949 subrun files processed (2012–2019), 1,789 physics frames in merged file
- Local linefit pkls already generated from this data: `linefit_all_years.pkl`, `linefit_fixed_speed.pkl`, `linefit_pivot.pkl` (all in `/home/bench/.icevenv/DMIce_data/`)
- Linefit comparison plots already made: `linefit_all_years.png`, `angle_comparison.png`, `v_comparison.png`, `angular_resolution.png`, `angular_shift.png`

**Coincidence pipeline — 2020/2021 extension (launched 2026-04-01):**
- Screen `dmice_2020_2021` running on cobalt, log: `~/dmice_work/pipeline_2020_2021.log`
- Script: `~/dmice_work/run_2020_2021_pipeline.sh` (also at `/home/bench/dmice/run_2020_2021_pipeline.sh`)
- Step1 wrappers and condor submit files already existed for 2020–2021 at `condor/step1/2020` and `condor/step1/2021` — just needed submitting
- 2020: ~14,646 step1 jobs (LC23-highHV2-v9 format); 2021: ~16,860 jobs (DM34-SLC_no-DM2_v1 format)
- 2021 files say "no-DM2" — det2 absent, so det2 step1 jobs for 2021 will produce empty output (harmless)
- Pipeline chains: step1 condor (NPX) → step2 on cobalt → step3 condor (NPX) → merge
- When done: re-merge will update `all_dmice_coincidences_2011_2022.i3.zst` to cover 2012–2021
- After merge: scp to local `~/dmice_results/` and regenerate linefit pkls

**Project phases (from local notes):**
- Phase 1 (CURRENT): Validate pivot fit on targeted simulation — need Phase 1 plot showing pivot reduces angular error for DM-Ice-passing muons
- Phase 2: Characterize IC-only errors using real data — bin IC-only vs pivot shifts by zenith/n_doms/speed. `characterize_bias.py` not yet written.
- Phase 3: Build correction model (lookup table, linear regression, or ML)
- Phase 4: Cross-validate correction model
- Phase 5: Apply to full dataset

**Phase 1 pipeline (built 2026-03-25, in `/home/bench/dmice/`):**
- `run_phase1_pipeline.sh` — main orchestrator: loops runs 2000–2019, calls steps below
- `run_prometheus_to_i3.sh` — parquet→npz→i3 wrapper; auto-detects LOCAL vs NPX (checks for `~/dmice_work/`)
  - LOCAL: uses `/home/bench/.icevenv/i3/scripts/` + local IceTray env-shell
  - NPX: uses `~/dmice_work/` scripts + CVMFS IceTray env-shell
- `sim_linefit_comparison.py` — per-run linefit comparison, outputs CSV + plot
- `merge_phase1_results.py` — merges all run CSVs → `phase1_all_runs.csv` + `phase1_validation.png`
- GCD design: pass the sim i3 file as its own GCD (`-g I3FILE -i I3FILE`) since prometheus_to_i3.py embeds geometry

**Running the Phase 1 pipeline:**
- On NPX (parquets at `/data/user/bcharett/dmice_sim_output/`): `cd ~/dmice && bash run_phase1_pipeline.sh` in a screen session
- Locally: sync parquets from NPX to `~/dmice_sim_output/run_0XXXX/XXXX_photons.parquet` first

**Current status (2026-03-25):**
- Condor clusters 10323163 + 10327036 submitted 2026-03-23; both should be complete by now (~6hr max ETA)
- Need to verify all 20 parquets exist on NPX, then run Phase 1 pipeline
- Next: run pipeline → get `phase1_validation.png` showing pivot reduces angular error
- After Phase 1: start `characterize_bias.py` for Phase 2 (can use existing pkl files in `/home/bench/.icevenv/DMIce_data/`)

**BLO (BlueLightOrchestra) — Julia-based Cherenkov simulation:**
- Location: `~/.icevenv/BLO/` (Julia package)
- Uses PROPOSAL (particle propagation, via PyCall) + PPC (photon propagation)
- CPU PPC: `resources/PPC_executables/PPC/ppc` (~30s/event)
- GPU/CUDA PPC: `resources/PPC_executables/PPC_CUDA/ppc` (~1.5s/event)
- Python interface via `juliacall` — see `examples/example.py`
- Key calls: `BLO.propagate(p_init, dist)` → losses; `BLO.run_ppc(p_init, losses; use_gpu=True)` → hits; `BLO.process_hits(hits)` → Table with pos/time/nhits/string_id/sensor_id per DOM
- `process_hits` returns a TypedTable; `sum(uhits.nhits)` = total DOM hits

**WARD machine:**
- Lab computer with a good GPU
- BLO simulation speed: ~1s/event (GPU)
- Use for interactive GPU simulations instead of NPX (which is CPU-only submit node)
- BLO setup on WARD: clone dmice repo → follow `BLO/README.md` (install Julia, clone BLO, Pkg.instantiate(), pip install juliacall, build PPC_CUDA binary)
- No Prometheus-specific WARD script yet; `simulate_muons.py` is NPX/CPU only

**Batch upgoing muon sim script (2026-03-27):**
- `/home/bench/dmice/batch_dm_ice_sim.py` (also at `BLO/batch_dm_ice_sim.py` in GitHub repo)
- Collects 100 upgoing through-going muon events with >200 DOM hits
- Uses BLO with GPU, injection parameters matching `simulate_dm_ice_through.py`
- Output: `~/dmice_work/output/blo_muons_200hits.npz`
- Run on WARD: `screen -S dmice_sim && python ~/dmice/BLO/batch_dm_ice_sim.py`
