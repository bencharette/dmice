# DMice Scripts Reference

Status: ✅ active/validated | ⚠️ partial/superseded | ❌ broken/tabled | 🗄️ archive

---

## Core Library

| Script | Purpose | Status | Notes |
|--------|---------|--------|-------|
| `blo_python.py` | BLO Python port: PROPOSAL propagation + GPU/CPU PPC photon sim. Core library imported by all sim scripts. Includes `smt8_trigger()`, `_hlc_mask()`, `process_hits()`, `run_ppc()`. | ✅ | Run on WARD. `BLO_PPC_EXE` + `BLO_PPC_TABLES` env vars required. |
| `dmice_likelihood.py` | Combined IC Pandel + DM-Ice Gaussian likelihood functions. | ❌ | Approximate uniform-ice Pandel diverges vs SpiceMie. See TABLED-01 in RESULTS.md. |

---

## Simulation Scripts

| Script | Purpose | Status | Machine | Output |
|--------|---------|--------|---------|--------|
| `simulate_muons_binned.py` | **Primary sim**: N events per log-spaced energy bin (100 GeV–100 TeV). Args: `--n-per-bin`, `--output`. Includes SMT8 trigger flag per event. | ✅ | WARD | `muons_binned_{N}ev.npz` |
| `simulate_muons_test_bin0.py` | 200-event quick test, Bin 0 only (100–398 GeV). | ✅ | WARD | `muons_test_bin0_200ev.npz` |
| `simulate_muons_offset.py` | Muons aimed at DM-Ice with configurable miss distance. Used for offset sensitivity studies. | ✅ | WARD | `muons_offset_{D}m_{N}ev.npz` |
| `simulate_muons.py` | Early targeted sim (DM-Ice aimed, Prometheus-based). | ⚠️ | WARD | Superseded by `simulate_muons_binned.py` |
| `batch_dm_ice_sim.py` | Batch BLO sim collecting events with >200 DOM hits. Upgoing only. | ⚠️ | WARD | Superseded |
| `simulate_dm_ice_through.py` | Prometheus dark matter ice-through simulation (exotic DM particles). | 🗄️ | — | Prometheus-based, separate from BLO pipeline |
| `run_blo_2events.py` | 2-event smoke test of BLO+PPC pipeline. | ✅ | WARD | `blo_dmice_targeted_det1det2_both_2events.npz` |

---

## Reconstruction Scripts

| Script | Purpose | Status | Machine | Output |
|--------|---------|--------|---------|--------|
| `run_sim_all_recos.py` | **Primary reco**: LineFit + Pivot LineFit + SPEFit + MPEFit (all seeded from Pivot) on BLO npz. Requires IceTray. | ✅ | Cobalt | `comparison/sim_all_recos.csv` |
| `run_splinempe_pivot.py` | SplineMPE + Pivot SPEFit/MPEFit on BLO npz. Heavier than `run_sim_all_recos.py`. | ✅ | Cobalt | `splinempe_pivot_comparison.csv` |
| `run_linefit.py` | Standalone LineFit only on BLO npz. | ✅ | Cobalt | — |
| `run_pivot_mpefit.py` | Pivot MPEFit only (fast single-fit version). | ✅ | Cobalt | — |
| `run_all_recos_real.py` | All recos on real IceCube coincidence data. | ✅ | Cobalt | `real_all_recos.csv` |
| `run_truncated_energy.py` | TruncatedEnergy estimator on reco'd events. | ⚠️ | Cobalt | — |
| `dmice_pivot_refit.py` | Appends synthetic DM-Ice pivot hit to IC pulses, reruns LineFit. Early prototype. | ⚠️ | Cobalt | Superseded by `run_sim_all_recos.py` |
| `inject_dmice_times.py` | Injects `dm_t_injected_ns` from npz as I3Double into i3 frames. Helper used by reco pipeline. | ✅ | Cobalt | — |

---

## Plotting Scripts

| Script | Purpose | Status | Output |
|--------|---------|--------|--------|
| `plot_smt8_efficiency.py` | SMT8 trigger efficiency vs energy with binomial errors. | ✅ | `smt8_efficiency_vs_energy.png` |
| `plot_sim_all_recos.py` | Angular error comparison across all fit types from `sim_all_recos.csv`. | ✅ | `ang_err_vs_energy.png` |
| `plot_ang_err_vs_energy_3panels.py` | 3-panel angular error vs energy (LF / SPE / MPE). | ✅ | — |
| `plot_energy_vs_ang_err.py` | Scatter: energy vs angular error per fit. | ✅ | `energy_vs_ang_err_*.png` |
| `plot_pivot_mpefit_comparison.py` | Side-by-side standard vs pivot MPEFit. | ✅ | `mpe_pivot_comparison.png` |
| `plot_sim_distributions.py` | DOM hit distributions from BLO simulation. | ✅ | — |
| `plot_blo_truth_distributions.py` | MC truth energy/zenith/azimuth distributions. | ✅ | — |
| `plot_hits_vs_energy_comparison.py` | DOM hits vs energy across sim batches. | ✅ | — |
| `plot_median_comparison.py` | Median angular error table plots across analyses. | ✅ | — |
| `plot_dmice_timing_diagram.py` | NaI scintillation timing model diagrams (slides). | ✅ | `dmice_timing_slide*.png` |
| `plot_real_recos.py` | Angular error distributions from real data recos. | ✅ | `real_recos_*.png` |
| `plot_event_display.py` | Simple event display (DOM positions + hit times). | ✅ | `event_display_ev*.png` |
| `replot_benchmark.py` | Regenerate benchmark plots from existing CSV (no IceTray). | ✅ | — |
| `npz_linefit_sim_comparison.py` | LineFit vs Pivot LineFit on BLO npz (no IceTray). | ✅ | — |
| `npz_linefit_comparison.py` | Same but on Prometheus parquet output. | ⚠️ | Prometheus-era |
| `compare_splinempe_seeds.py` | Compare SplineMPE results from different seeds. | ⚠️ | — |

---

## Data Pipeline Scripts

| Script | Purpose | Status | Machine |
|--------|---------|--------|---------|
| `merge_2020_2021.py` | Merge 2020+2021 IceCube coincidence i3 files. | ✅ | Cobalt |
| `rebuild_master_i3.py` | Rebuild master coincidence i3 from all 2012–2021 step3 files. | ✅ | Cobalt |
| `score_real_coincidences.py` | Score/filter real DM-Ice coincidences. | ✅ | Cobalt |
| `build_dmice_timing_model.py` | Build NaI scintillation timing model from real 2012 data. Output: `dmice_timing_model.npz`. | ✅ | Cobalt |
| `parquet_to_npz.py` | Prometheus parquet → npz (Step 1, system python3). | ✅ | NPX |
| `prometheus_to_i3.py` | npz → i3 (Step 2, IceTray). | ✅ | NPX/Cobalt |
| `merge_phase1_results.py` | Merge phase 1 pipeline output plots. | ⚠️ | — |
| `step1_resubmit_2020_2021.py` | Resubmit failed Condor jobs for 2020–2021 pipeline. | 🗄️ | NPX |
| `step3_submit_2020_2021.py` | Submit step3 Condor jobs for 2020–2021. | 🗄️ | NPX |
| `check_track_through_dmice.py` | Verify BLO tracks geometrically pass through DM-Ice crystals. Diagnostic. | ✅ | — |
| `sim_linefit_comparison.py` | LineFit vs Pivot on BLO i3 files (early version). | ⚠️ | Cobalt |
| `npz_to_linefit_csv.py` | Export LineFit results from npz to CSV. | ✅ | — |

---

## Shell / Submission Scripts

| Script | Purpose | Status | Machine |
|--------|---------|--------|---------|
| `run_2020_2021_pipeline.sh` | Orchestrate full DM-Ice coincidence pipeline (parquet → i3 → recos). | ✅ | NPX (screen) |
| `pipeline_2020_2021_followup.sh` | Follow-up jobs for 2020–2021 pipeline. | 🗄️ | NPX |
| `phase34_run.sh` | Run phase 3+4 of pipeline. | 🗄️ | NPX |
| `launch_offset_sims.sh` | Launch offset sensitivity sim batch on WARD. | ✅ | WARD |
| `overnight_run.sh` | Wrapper to run long jobs overnight. | 🗄️ | — |
| `run_phase1_pipeline.sh` | Phase 1 validation: parquet → i3 → sim comparison for 20 runs. | 🗄️ | NPX |
| `run_prometheus_condor.sh` | Submit Prometheus sim via Condor. | 🗄️ | NPX |
| `run_prometheus_to_i3.sh` | Run parquet→npz→i3 pipeline. | 🗄️ | NPX |
| `run_benchmark.sub` | Condor submit file for benchmark jobs. | 🗄️ | NPX |
| `run_all_recos_real.sub` | Condor submit file for real data reco jobs. | ✅ | NPX |
| `simulate_muons.sub` | Condor submit file for muon sim jobs. | ⚠️ | NPX |
| `start-work.sh` | Session startup script (env setup). | ✅ | Local |
| `sync-memory.sh` | Sync memory/handoff files across machines via git. | ✅ | Local |

---

## Visualisation / IceTray Helpers

| Script | Purpose | Status |
|--------|---------|--------|
| `steamshovel_artists.py` | ICLineFitArtist + PivotLineFitArtist for Steamshovel. Auto-registered on Cobalt via startup.py. | ✅ |
| `load_dmice_artists.py` | Loads custom artists into a running Steamshovel session. | ✅ |

---

## Repacking Utility

| Script | Purpose | Status | Notes |
|--------|---------|--------|-------|
| `BLO/repack_npz.py` | Convert dtype=object ragged arrays → flat+offsets format for NumPy 1.x compat on Cobalt. **Run on WARD after every sim before scp to Cobalt.** | ✅ | Required step — Cobalt IceTray uses NumPy 1.x |
