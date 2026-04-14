# DMice Simulation Results Ledger

One row per simulation run or analysis result. Status: ✅ validated | ⚠️ partial | ❌ broken/tabled.

---

## Simulation Runs

### SIM-01 — BLO 2-event test
- **Script**: `run_blo_2events.py`
- **Machine**: WARD (GPU PPC)
- **Output**: `~/dmice_work/output/blo_dmice_targeted_det1det2_both_2events.npz`
- **Events**: 2 (one per DM-Ice detector)
- **Purpose**: Smoke-test BLO Python port + PPC pipeline
- **Status**: ✅ PPC runs, DM-Ice strings 87/88 produce hits

---

### SIM-02 — BLO 200-event initial run
- **Script**: `simulate_muons_offset.py` (offset=0)
- **Machine**: WARD (GPU PPC)
- **Output**: `~/dmice_work/output/blo_muons_200hits.npz` + `_repacked.npz`
- **Events**: 200, uniform energy, tracks aimed at DM-Ice
- **Purpose**: First full sim batch, used for LineFit/SPEFit/MPEFit baseline
- **Status**: ✅

---

### SIM-03 — BLO 1000-event targeted (det1+det2)
- **Script**: `simulate_muons_offset.py` (targeted)
- **Machine**: WARD (GPU PPC)
- **Output**: `~/dmice_work/output/blo_dmice_targeted_det1det2_both_1000events.npz` + `_repacked.npz`
- **Events**: 1000, alternating det1/det2 targeting
- **Purpose**: Larger batch for SplineMPE + Pivot comparison
- **Status**: ✅

---

### SIM-04 — 5000-event binned (5 energy bins)
- **Script**: `simulate_muons_binned.py --n-per-bin 1000`
- **Machine**: WARD (GPU PPC)
- **Output**: `~/dmice_work/output/muons_binned_5000ev.npz` + `_repacked.npz` (36–46 MB)
- **Events**: 5000 (1000/bin × 5 bins, 100 GeV – 100 TeV)
- **Seed**: 42
- **Purpose**: Main SplineMPE pivot comparison dataset
- **Status**: ✅ Used in ANALYSIS-03 (splinempe_pivot_comparison.csv)

---

### SIM-05 — 200-event Bin 0 test
- **Script**: `simulate_muons_test_bin0.py`
- **Machine**: WARD (GPU PPC)
- **Output**: `~/dmice_work/output/muons_test_bin0_200ev.npz` + `_repacked.npz`
- **Events**: 200, Bin 0 only (100–398 GeV)
- **Seed**: 99
- **Purpose**: Quick test of SMT8 trigger + Pivot seed angular resolution
- **Key result**: SMT8 efficiency 77/200 = **38.5%** at lowest energy bin
- **Status**: ✅

---

### SIM-06 — 1000-event binned (200/bin trigger efficiency)
- **Script**: `simulate_muons_binned.py --n-per-bin 200`
- **Machine**: WARD (GPU PPC)
- **Output**: `~/dmice_work/output/muons_binned_1000ev.npz` + `_repacked.npz` (9.1 MB)
- **Events**: 1000 (200/bin × 5 bins)
- **Seed**: 42
- **Purpose**: SMT8 trigger efficiency vs energy; full reco comparison
- **SMT8 results**:
  | Bin | Energy | Efficiency |
  |-----|--------|-----------|
  | 0 | 100–398 GeV | 37.5% |
  | 1 | 398 GeV–1.58 TeV | 81.5% |
  | 2 | 1.58–6.31 TeV | 96.0% |
  | 3 | 6.31–25.1 TeV | 99.5% |
  | 4 | 25.1–100 TeV | 98.5% |
  | All | — | 82.6% |
- **Plot**: `~/dmice_work/output/smt8_efficiency_vs_energy.png`
- **Status**: ✅

---

## Reconstruction Analyses

### ANALYSIS-01 — LineFit/SPEFit/MPEFit on SIM-02
- **Script**: `run_linefit.py`
- **Input**: SIM-02 repacked npz
- **Machine**: Cobalt (IceTray v1.12.1)
- **Purpose**: Baseline angular resolution without DM-Ice constraint
- **Status**: ✅

---

### ANALYSIS-02 — SplineMPE + Pivot on SIM-03
- **Script**: `run_splinempe_pivot.py`
- **Input**: SIM-03 repacked npz
- **Machine**: Cobalt
- **Purpose**: SplineMPE with and without Pivot LineFit seed; combined DM likelihood
- **Status**: ⚠️ SplineMPE + Pivot validated; Combined DM likelihood diverges (see TABLED-01)

---

### ANALYSIS-03 — Full reco comparison on SIM-04 (5000 events)
- **Script**: `run_splinempe_pivot.py`
- **Input**: `muons_binned_5000ev_repacked.npz`
- **Machine**: Cobalt
- **Output**: `~/dmice_work/output/splinempe_pivot_comparison.csv` (5000 rows)
- **Plots**: `splinempe_pivot_comparison*.png` (4 plots)
- **Key results** (median angular error):
  | Fit | Median err |
  |-----|-----------|
  | LineFit | ~10° |
  | Pivot LineFit | ~1–2° |
  | SPEFit | ~3° |
  | MPEFit (std seed) | ~3° |
  | MPEFit (pivot seed) | ~2.9° |
- **Note**: `mpe_dm_ang_err` column absent — combined DM fit was never run on this batch
- **Status**: ✅ (SplineMPE results validated)

---

### ANALYSIS-04 — LineFit/SPE/MPE + Pivot on SIM-06 (1000 events)
- **Script**: `run_sim_all_recos.py`
- **Input**: `muons_binned_1000ev_repacked.npz`
- **Machine**: Cobalt
- **Output**: `~/dmice_work/output/comparison/sim_all_recos.csv`
- **Plot**: `~/dmice_work/output/ang_err_vs_energy.png`
- **Key results** (median angular error, all 5 bins combined):
  | Fit | Valid | Median err |
  |-----|-------|-----------|
  | LineFit | 878 | 4.48° |
  | Pivot LineFit | 869 | 1.71° |
  | SPEFit | 859 | 1.01° |
  | Pivot SPEFit | 859 | 1.01° |
  | MPEFit | 842 | **0.55°** |
  | Pivot MPEFit | 842 | **0.55°** |
  | Combined SPE/MPE | — | ❌ ~60° diverged |
- **Note**: SPEFit = Pivot SPEFit and MPEFit = Pivot MPEFit because both are seeded from the same Pivot LineFit — expected behaviour
- **Status**: ✅

---

## Tabled / Dead Ends

### TABLED-01 — Combined DM-Ice + IceCube likelihood fit
- **Scripts**: `run_splinempe_pivot.py` (`DMCombinedFitModule`), `run_sim_all_recos.py`
- **Approach**: Pandel PDF (IC) + Gaussian(280 ns, 81 ns) DM-Ice term, Nelder-Mead optimizer
- **Result**: ~60° median error — worse than LineFit
- **Root cause**: `_PANDEL_LA=98 m`, `_PANDEL_LS=30 m` (uniform ice approximation) does not match SpiceMie layered model used in PPC simulation. Optimizer diverges.
- **Fix required**: Replace Pandel with SplineMPE tables or MilliPede likelihood. Non-trivial.
- **Status**: ❌ Tabled — use SPEFit/MPEFit + Pivot seed instead

---

## Real Data Analyses

### REAL-01 — Real coincidence reconstruction
- **Script**: `run_all_recos_real.py`
- **Input**: Real IceCube data (2020–2021 coincidences)
- **Output**: `~/dmice_work/output/real_all_recos.csv`
- **Plots**: `real_recos_*.png`
- **Status**: ✅ (used as cross-check for sim pipeline)
