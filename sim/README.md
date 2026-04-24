# sim/ — Simulation Scripts

Generate muon events through the DM-Ice detector using Prometheus + BLO/PPC.

## Scripts

| File | Description |
|------|-------------|
| `simulate_muons.py` | Main muon simulation — runs Prometheus + PPC photon propagation, outputs `.npz` |
| `simulate_muons_binned.py` | Binned simulation for systematic studies (energy/zenith bins) |
| `simulate_muons_offset.py` | Simulation with lateral offset (`d_perp`) from detector axis |
| `simulate_muons_test_bin0.py` | Quick test run for a single energy bin |
| `simulate_dm_ice_through.py` | Simulate muon tracks that pass through the DM-Ice volume |
| `batch_dm_ice_sim.py` | Batch launcher for multiple simulation configurations |
| `inject_dmice_times.py` | Inject NaI hit times into existing i3 frames for timing model tests |
| `prometheus_to_i3.py` | Convert Prometheus `.npz` photon output to IceTray `.i3` format |
| `parquet_to_npz.py` | Convert Prometheus `.parquet` output to `.npz` |
| `run_prometheus_condor.sh` | Condor job wrapper — calls `simulate_muons.py` (keep with `.sub`) |
| `simulate_muons.sub` | HTCondor submit file for `run_prometheus_condor.sh` |
| `launch_offset_sims.sh` | Launch `simulate_muons_offset.py` for multiple offsets on WARD |

## Usage

```bash
# Run simulation locally (WARD GPU env required for PPC)
python3 ~/dmice/sim/simulate_muons.py

# Submit to Condor (NPX)
condor_submit ~/dmice/sim/simulate_muons.sub

# Launch offset scans on WARD
bash ~/dmice/sim/launch_offset_sims.sh
```
