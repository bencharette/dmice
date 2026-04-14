# BLO Simulation Pipeline — Command Reference

## Machine usernames
- **WARD**: `bencharette@WARD`
- **Cobalt**: `bcharett@cobalt-14`
- **Local**: `/home/ben/`

---

## Step 1: Run simulation on WARD (GPU)

```bash
ssh ward
nohup bash -c '
  BLO_PPC_EXE=~/.icevenv/BLO/resources/PPC_executables/PPC_CUDA/ppc \
  BLO_PPC_TABLES=~/.icevenv/BLO/resources/PPC_tables/south_pole \
  python3 ~/dmice/simulate_muons_binned.py \
  > ~/dmice_work/output/sim5000.log 2>&1
' &
```

**Check progress:**
```bash
ssh ward "tail -3 ~/dmice_work/output/sim5000.log"
```

**Output:** `~/dmice_work/output/muons_binned_5000ev.npz`

---

## Step 2: Repack NPZ on WARD

Required for numpy 1.x compatibility on cobalt. Run after sim finishes.

```bash
ssh ward "python3 ~/dmice/BLO/repack_npz.py \
  ~/dmice_work/output/muons_binned_5000ev.npz \
  ~/dmice_work/output/muons_binned_5000ev_repacked.npz"
```

**Output:** `~/dmice_work/output/muons_binned_5000ev_repacked.npz`

---

## Step 3: Copy repacked NPZ from WARD to cobalt

Run from local machine or cobalt.

```bash
# From local:
scp ward:~/dmice_work/output/muons_binned_5000ev_repacked.npz \
    bcharett@cobalt-14:~/dmice_work/output/

# From cobalt directly:
ssh cobalt-14 "scp bencharette@WARD:~/dmice_work/output/muons_binned_5000ev_repacked.npz \
    ~/dmice_work/output/"
```

---

## Step 4: Copy analysis scripts to cobalt

```bash
scp ~/dmice/run_splinempe_pivot.py bcharett@cobalt-14:~/dmice/
```

---

## Step 5: Update NPZ path in analysis script

In `run_splinempe_pivot.py`, line ~27:
```python
NPZ_FILE = os.path.expanduser(
    "~/dmice_work/output/muons_binned_5000ev_repacked.npz")
```

---

## Step 6: Run reconstruction benchmark on cobalt

```bash
ssh cobalt-14
/cvmfs/icecube.opensciencegrid.org/py3-v4.3.0/RHEL_9_x86_64/metaprojects/icetray/v1.12.1/env-shell.sh \
  python3 ~/dmice/run_splinempe_pivot.py \
  > ~/dmice_work/output/splinempe_pivot_comparison.log 2>&1
```

**Or in a screen session (recommended for long runs):**
```bash
ssh cobalt-14 "screen -dmS benchmark bash -c '
/cvmfs/icecube.opensciencegrid.org/py3-v4.3.0/RHEL_9_x86_64/metaprojects/icetray/v1.12.1/env-shell.sh \
  python3 ~/dmice/run_splinempe_pivot.py \
  > ~/dmice_work/output/splinempe_pivot_comparison.log 2>&1
'"
```

**Check progress:**
```bash
ssh cobalt-14 "tail -5 ~/dmice_work/output/splinempe_pivot_comparison.log"
```

**Outputs:**
- `~/dmice_work/output/splinempe_pivot_comparison.csv`
- `~/dmice_work/output/splinempe_pivot_comparison.png`
- `~/dmice_work/output/splinempe_pivot_comparison_per_bin.png`

---

## Step 7: Copy results back to local

```bash
scp bcharett@cobalt-14:~/dmice_work/output/splinempe_pivot_comparison*.{png,csv} \
    ~/dmice_work/output/
```

---

## Full pipeline one-liner reference

```
WARD (GPU sim) → WARD (repack) → cobalt (reco benchmark) → local (plots)
```

---

## Environment variables for BLO on WARD

| Variable | Value |
|----------|-------|
| `BLO_PPC_EXE` | `~/.icevenv/BLO/resources/PPC_executables/PPC_CUDA/ppc` |
| `BLO_PPC_TABLES` | `~/.icevenv/BLO/resources/PPC_tables/south_pole` |
| `BLO_GEO_FILE` | `~/.icevenv/BLO/resources/geofiles/icecube_with_dmice.geo` |

---

## IceTray environment on cobalt

```bash
/cvmfs/icecube.opensciencegrid.org/py3-v4.3.0/RHEL_9_x86_64/metaprojects/icetray/v1.12.1/env-shell.sh
```

Call as a wrapper: `env-shell.sh python3 script.py` — do NOT source it over SSH.

---

## Key file locations

| File | Machine | Path |
|------|---------|------|
| Sim script | WARD + local | `~/dmice/simulate_muons_binned.py` |
| Repack script | WARD + local | `~/dmice/BLO/repack_npz.py` |
| Benchmark script | cobalt + local | `~/dmice/run_splinempe_pivot.py` |
| 200-event repacked NPZ | cobalt | `~/dmice_work/output/muons_binned_200ev_repacked.npz` |
| 5000-event repacked NPZ | cobalt (after scp) | `~/dmice_work/output/muons_binned_5000ev_repacked.npz` |
| Timing model (BLO) | cobalt | `~/dmice_work/output/dmice_timing_model.npz` |
| Timing model (real data) | cobalt | `~/dmice_work/output/dmice_timing_model_calibrated.npz` |
| Benchmark CSV | cobalt | `~/dmice_work/output/splinempe_pivot_comparison.csv` |
| Benchmark plots | cobalt + local | `~/dmice_work/output/splinempe_pivot_comparison*.png` |
