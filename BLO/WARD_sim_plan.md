# Plan: BLO DM-Ice Through-Going Muon Simulation on WARD

## Context
The lab's Linux simulation desktop (WARD) has a good GPU (~1s/event with BLO). We want to simulate muons that physically pass through the DM-Ice NaI scintillator detectors using BLO (BlueLightOrchestra), which is faster than Prometheus. This covers both downgoing (atmospheric, like simulate_muons.py) and upgoing muons. Output should be NPZ (for analysis) and I3 (for IceTray/steamshovel).

BLO is not yet installed on WARD.

---

## Phase 1: Install BLO on WARD

Follow `BLO/README.md` exactly. Steps:

```bash
# 1. Install Julia
curl -fsSL https://install.julialang.org | sh

# 2. Clone BLO
mkdir -p ~/.icevenv && git clone https://github.com/kcarloni/BlueLightOrchestra.jl ~/.icevenv/BLO

# 3. Instantiate Julia packages
cd ~/.icevenv/BLO && julia --project=. -e 'using Pkg; Pkg.instantiate()'

# 4. Python bridge
pip install juliacall

# 5. Copy DM-Ice geometry
cp ~/dmice/BLO/icecube_with_dmice.geo ~/.icevenv/BLO/resources/geofiles/

# 6. Compile GPU PPC binary
cd ~/.icevenv/BLO/resources/PPC_executables/PPC_CUDA && make
```

Verify with: `python3 -c "import juliacall; print('OK')"`

---

## Phase 2: Write the Simulation Script

**New file:** `~/dmice/BLO/batch_dm_ice_targeted_sim.py`

Based on `batch_dm_ice_sim.py` but with targeted DM-Ice injection geometry.

### Key design differences from batch_dm_ice_sim.py

| | batch_dm_ice_sim.py | new script |
|---|---|---|
| Direction | Upgoing only (0–90°) | Both (up + down) |
| Injection | Through IceCube center | Back-projected from DM-Ice position |
| Energy (downgoing) | — | 1 TeV – 1 PeV, gamma=2 |
| Energy (upgoing) | 100 GeV – 1 PeV, gamma=1 | same |
| Hit threshold | >200 total DOM hits | >50 (downgoing hits fewer DOMs) |
| Output | NPZ only | NPZ + flag for I3 conversion |

### Injection geometry

For both directions, back-project 2000m from DM-Ice position along anti-momentum direction:

```python
# DM-Ice positions (IceCube coords, meters)
DMICE_POSITIONS = {
    "det1": np.array([31.25, -72.93, -511.05]),
    "det2": np.array([-334.80, -424.50, -511.26]),
}

dx = np.sin(zenith) * np.cos(azimuth)
dy = np.sin(zenith) * np.sin(azimuth)
dz = np.cos(zenith)
direction = np.array([dx, dy, dz])

# Inject 2000m back along track from DM-Ice
BACKPROJECT_DIST = 2000.0  # meters
start_pos = dmice_pos - BACKPROJECT_DIST * direction
```

**Downgoing:** zenith 130–170° (same as simulate_muons.py), dz < 0
**Upgoing:** zenith 10–50°, dz > 0

### CLI args
```
--nevents N         Number of accepted events (default: 100)
--det {1,2,both}    Target DM-Ice detector (default: both)
--direction {up,down,both}  Muon direction (default: both)
--outdir PATH       Output directory (default: ~/dmice_work/output/)
--no-gpu            Disable GPU (fall back to CPU)
```

### Output

**NPZ** (`blo_dmice_targeted_{det}_{direction}_{N}events.npz`):
Same structure as `batch_dm_ice_sim.py` — per-event arrays of energy, zenith, azimuth, n_hits, n_doms, plus per-DOM ragged arrays (x, y, z, t, nhits, string, sensor). Also stores `det_id` and `direction_type` per event.

**I3 conversion:** Add `--convert-i3` flag that calls `prometheus_to_i3.py` logic post-hoc. Since IceTray may not be on WARD, this step is designed to run separately on cobalt or local machine by passing the NPZ through the existing `prometheus_to_i3.py` script.

---

## Phase 3: Test Run

```bash
# Quick test on WARD (10 events, det1, downgoing only)
python ~/dmice/BLO/batch_dm_ice_targeted_sim.py \
    --nevents 10 --det 1 --direction down
```

Check: output NPZ exists, events have reasonable hit counts (expect 50–500 DOMs for 1–3 TeV downgoing muon through IC86+DM-Ice).

---

## Phase 4: I3 Conversion (on cobalt or local)

```bash
# On machine with IceTray:
source /path/to/env-shell.sh python3 prometheus_to_i3.py \
    --npz blo_dmice_targeted_det1_down_100events.npz \
    --geo ~/dmice/BLO/icecube_with_dmice.geo \
    --output blo_dmice_targeted_det1_down_100events.i3.zst
```

Uses existing `prometheus_to_i3.py` (at `~/dmice/prometheus_to_i3.py` locally or `~/dmice_work/prometheus_to_i3.py` on NPX/cobalt).

---

## Critical Files

- `BLO/batch_dm_ice_sim.py` — reference/template for new script
- `BLO/README.md` — WARD install instructions
- `simulate_muons.py` — reference for downgoing injection geometry and DM-Ice positions
- `prometheus_to_i3.py` — I3 conversion (run on cobalt/local, not WARD)
- `BLO/icecube_with_dmice.geo` — geometry file, must be on WARD

---

## Verification

1. BLO install: `python3 -c "import juliacall; print('OK')"` on WARD
2. Test sim produces NPZ with correct structure
3. At least one event has DOM hits near DM-Ice z position (~-511m)
4. I3 file opens in steamshovel showing track passing through DM-Ice sphere
