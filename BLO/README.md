# BLO Setup for WARD

Setup instructions for running [BlueLightOrchestra](https://github.com/kcarloni/BlueLightOrchestra.jl) (BLO) on WARD.

BLO uses PROPOSAL (particle propagation) + PPC (photon propagation via GPU). At ~1s/event on WARD's GPU.

---

## Files in this folder

| File | Purpose |
|------|---------|
| `batch_dm_ice_sim.py` | Batch simulation: collects 100 upgoing muon events with >200 DOM hits |
| `icecube_with_dmice.geo` | Custom IceCube+DM-Ice geometry file (not in BLO repo by default) |

---

## One-time setup on WARD

### 1. Install Julia

```bash
curl -fsSL https://install.julialang.org | sh
# restart shell or: source ~/.bashrc
julia --version   # verify
```

### 2. Clone BLO

```bash
mkdir -p ~/.icevenv
git clone https://github.com/kcarloni/BlueLightOrchestra.jl ~/.icevenv/BLO
```

### 3. Instantiate the Julia environment

```bash
cd ~/.icevenv/BLO
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

This downloads all Julia dependencies (PROPOSAL, PPC, etc.). Takes a few minutes first time.

### 4. Install the Python bridge

```bash
pip install juliacall
```

### 5. Copy the DM-Ice geofile into BLO resources

```bash
cp ~/dmice/BLO/icecube_with_dmice.geo ~/.icevenv/BLO/resources/geofiles/
```

### 6. Build the GPU PPC binary

The PPC CUDA binary must be compiled on WARD (GPU-specific):

```bash
cd ~/.icevenv/BLO/resources/PPC_executables/PPC_CUDA
make
ls ppc   # verify
```

If `make` fails, check that CUDA toolkit is installed: `nvcc --version`

---

## Running the simulation

```bash
screen -S dmice_sim
python ~/dmice/BLO/batch_dm_ice_sim.py
```

Output: `~/dmice_work/output/blo_muons_200hits.npz` (100 events, >200 hits each)

To reload:
```python
import numpy as np
d = np.load('~/dmice_work/output/blo_muons_200hits.npz', allow_pickle=True)
d['energy_GeV']   # shape (100,)
d['dom_x'][0]     # x positions of hit DOMs in event 0 [m]
```

---

## Quick API reference

```python
import juliacall
from juliacall import Main as jl

jl.seval("""
using Pkg; Pkg.activate("/path/to/.icevenv/BLO")
using BlueLightOrchestra
using BlueLightOrchestra.AstroParticleUnits
""")
BLO = jl.BlueLightOrchestra
jlx = getattr(jl, "*")

# Create a particle
p = BLO.ParticleState(
    jlx(1000, jl.GeV),          # energy
    jlx([0, 0, -1.3], jl.km),   # position (depth coords)
    [0, 0, -1],                  # direction unit vector
    jl.PDGID(13),                # PDG: 13 = muon
    jlx(0.0, jl.ns),             # time
)

losses = BLO.propagate(p, jlx(1.5, jl.km))
hits   = BLO.run_ppc(p, losses, suppress_error=True, use_gpu=True)
uhits  = BLO.process_hits(hits)   # TypedTable: pos, time, nhits, string_id, sensor_id

total_hits = int(sum(uhits.nhits))
```
