# reco/ — Reconstruction Scripts

Run track reconstructions on simulated and real DM-Ice coincidence events.

## Scripts

| File | Description |
|------|-------------|
| `blo_python.py` | Python port of Blue Light Orchestra reconstruction |
| `build_dmice_timing_model.py` | Fit the NaI timing model (Gaussian mean/sigma from sim) |
| `check_track_through_dmice.py` | Check whether a reconstructed track passes through DM-Ice volume |
| `dmice_likelihood.py` | NaI likelihood term: Gaussian timing + geometric d_perp |
| `dmice_pivot_refit.py` | Refit tracks using the pivot-point SplineMPE method |
| `run_2d_vtx_mpe.py` | 2D vertex + MPE joint reconstruction |
| `run_all_recos_real.py` | Run all reconstructions on real coincidence i3 file |
| `run_all_recos_real.sub` | Condor submit file for `run_all_recos_real.py` |
| `run_benchmark.sub` | Condor submit file for benchmarking `run_splinempe_pivot.py` |
| `run_blo_2events.py` | Test BLO on 2 hand-picked events |
| `run_itermpe_events.py` | Run iterative MPE on selected events |
| `run_linefit.py` | Run standard LineFit reconstruction |
| `run_pivot_mpefit.py` | Run pivot-point MPEfit |
| `run_sim_all_recos.py` | Run all reconstructions on simulated events |
| `run_splinempe_nai_lambda.py` | λ scan: SplineMPE + NaI likelihood (scipy refinement, Pandel IC term) |
| `run_nai_lambda_pandel.py` | λ scan variant: vertex-anchored Pandel + NaI Gaussian (newer approach) |
| `run_splinempe_nai_gulliver.py` | λ scan: SplineMPE + NaI via proper Gulliver I3EventLogLikelihood (correct spline IC term) |
| `run_splinempe_pivot.py` | SplineMPE with pivot-point seed |
| `run_truncated_energy.py` | Run truncated energy estimator |
| `score_real_coincidences.py` | Score quality of real coincidence events |
| `step1_resubmit_2020_2021.py` | Resubmit failed step1 Condor jobs for 2020-2021 |
| `step3_rerun_all_years.py` | Rerun step3 coincidence finder for all years |
| `step3_submit_2020_2021.py` | Submit step3 Condor jobs for 2020-2021 |
| `submit_reco_condor.py` | Generic Condor submitter for `run_sim_all_recos.py` chunks |
| `merge_reco_chunks.py` | Merge chunked Condor reco output into single file |
| `rebuild_master_i3.py` | Rebuild the master coincidence i3 file from scratch |
