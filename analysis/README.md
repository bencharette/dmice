# analysis/ — Analysis Scripts

Compare reconstructions, study coincidence rates, and extract physics results.

## Scripts

| File | Description |
|------|-------------|
| `compare_coinc_cuts.py` | Compare Gaussian timing cut vs geometric d_perp cut for coincidence selection |
| `compare_coinc_cuts_v2.py` | Updated version with additional cut variations |
| `compare_splinempe_seeds.py` | Compare SplineMPE seeds (LineFit vs pivot) on angular error |
| `count_coinc_per_year.py` | Count coincidence events per year for rate analysis |
| `extract_dm_amplitude.py` | Extract NaI pulse amplitude from ROOT files |
| `merge_2020_2021.py` | Merge 2020-2021 step3 coincidence i3 files into a single file |
| `merge_phase1_results.py` | Merge phase 1 CSV results from multiple runs |
| `ml_coinc_classifier.py` | ML classifier (BDT/RF) for coincidence vs background |
| `npz_linefit_comparison.py` | Compare LineFit results between npz and i3 pipelines |
| `npz_linefit_sim_comparison.py` | Compare sim LineFit vs real LineFit angular distributions |
| `npz_to_linefit_csv.py` | Convert npz event arrays to LineFit CSV format |
| `profile_era_keys.py` | Profile which i3 frame keys exist across different data eras |
| `sim_linefit_comparison.py` | Compare simulation LineFit output to truth directions |
