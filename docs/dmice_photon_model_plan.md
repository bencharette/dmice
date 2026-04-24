# Plan: DM-Ice NaI Photon Detection Model for Direct Likelihood

## Context

The DM-Ice pivot seeding approach (Pivot LineFit → Pivot MPEFit) doesn't improve reconstruction
accuracy beyond what MPEFit already achieves (~0.46° on clean BLO sim). The real gain comes
from treating DM-Ice timing as a **direct term in the reconstruction likelihood**:

```
log L_total = log L_IceCube(pulses | track) + log L_DM-Ice(t_DM | track)
```

For this to work, we need a photon detection model for the DM-Ice NaI crystal that maps
a track hypothesis → expected hit time + probability distribution. This plan describes how
to build that model using the BLO simulation we already have.

---

## What Physical Quantities Need Modeling

### 1. Expected hit time (geometric term)

Given a muon track with:
- Vertex (x₀, y₀, z₀) at time t₀
- Direction unit vector d̂ = (dx, dy, dz)
- DM-Ice detector position r_DM

The muon passes the point of closest approach (PCA) to r_DM at:
```
t_PCA = t₀ + [(r_DM - r₀) · d̂] / c
d_perp = ||(r_DM - r₀) - [(r_DM - r₀) · d̂] d̂||   # perpendicular distance
```

Cherenkov photons emitted at the PCA reach r_DM at geometric time:
```
t_geo = t_PCA + d_perp / (c · sin θ_c)  ≈ t_PCA + d_perp / 0.2998 · (1/sin41°)
```
(This is the same formula used in Pandel/SPE fits for IceCube DOMs.)

### 2. Timing residual distribution

Observed hit time = t_geo + Δt, where Δt includes:
- **Scattering delay**: photons scattered in ice arrive later (Pandel tail)
- **NaI scintillation delay**: ~250ns decay time, but Cherenkov photons are prompt and the
  PMT fires on first photon — this is a ~1ns smearing, negligible vs scattering
- **Electronic jitter**: few ns

The residual PDF `p(Δt | d_perp, zenith)` is the key output of the model.

### 3. Detection probability (amplitude term)

Not all tracks produce hits — probability that the NaI fires at all:
```
P(hit | d_perp, zenith) = efficiency × (photon yield integrated over visible track)
```
This enters the likelihood as a Poisson amplitude term: for expected n_exp hits,
L_amp ∝ Poisson(n_obs | n_exp).

---

## Data Sources Already Available

**BLO simulation output** (the 200-event npz + any future runs) contains:
- `dom_x, dom_y, dom_z, dom_t, dom_nhits, dom_string, dom_sensor` per event
- DM-Ice DOMs are in the geofile at z_BLO ≈ −2459 m (IceCube z ≈ −511 m)
- MC truth: `energy_GeV, zenith_rad, azimuth_rad` → full track direction
- Injection geometry: track passes through DM-Ice det1 or det2

**Key insight**: BLO simulation already propagates Cherenkov photons through the ice model
(PPC/SPICEMie) and records which DOMs fired and when. The DM-Ice entries in the hit arrays
*are* the photon detection signal. No new simulation is needed for the model — just extract
and fit the residual distribution from existing events.

---

## Step-by-Step Model Construction

### Step 1: Extract DM-Ice hits from BLO sim

Script: `~/dmice/build_dmice_timing_model.py`

```python
# For each event in the npz:
# 1. Identify which DOM entries correspond to DM-Ice detectors
#    (match dom_string/dom_sensor to known DM-Ice OMKeys, or by position)
# 2. For each DM-Ice hit, compute:
#    - MC truth PCA time: t_geo from truth track + DM-Ice position
#    - Timing residual: Δt = dom_t[j] - t_geo
#    - Perpendicular distance: d_perp from truth track to DM-Ice position
#    - cos(zenith) of the track
# 3. Record (d_perp, zenith, Δt, n_photons) for each event with a DM-Ice hit
```

DM-Ice DOM identification (from geofile `icecube_with_dmice.geo`):
- det1: string 87, DOM 1 (or whatever the geofile assigns)
- det2: string 88, DOM 1
- BLO positions: det1=[31.25, -72.93, -2459.12], det2=[-334.80, -424.50, -2459.33] (BLO z)
- IceCube z: add Z_OFFSET=1948.07 → det1_z ≈ -511.05 m

### Step 2: Compute geometric expected times

For each event with a DM-Ice hit, compute t_geo from MC truth:
```python
def t_geometric(track_vertex, track_dir, t0, dm_pos, c=0.2998, n_ice=1.3195):
    """
    Returns expected first-photon arrival time at dm_pos for a muon track.
    track_dir: unit vector (travel direction)
    n_ice: ice refractive index at Cherenkov frequency
    """
    theta_c = np.arccos(1.0 / n_ice)   # ≈ 41°
    r = dm_pos - track_vertex
    s = np.dot(r, track_dir)           # signed distance along track to PCA
    d_perp = np.sqrt(np.dot(r, r) - s**2)
    t_pca = t0 + s / c
    t_geo = t_pca + d_perp / (c * np.sin(theta_c))
    return t_geo, d_perp
```

Residual: `Δt = t_observed - t_geo`

### Step 3: Fit the residual distribution

With ~200 events (more if we run larger sims), fit `p(Δt | d_perp)`:

**Option A — Pandel-like (recommended):**
Use the same functional form as IceCube's Pandel PDF:
```
p(Δt | d, λ_a, λ_s) ∝ Pandel(Δt; d, λ_a, λ_s)
```
where λ_a = absorption length, λ_s = scattering length. These are SPICEMie ice model parameters.
Since BLO already uses SPICEMie, the photons we extract already encode these properties — fit
the effective λ_a, λ_s by matching the Pandel function to the observed Δt histogram.

**Option B — Lookup table:**
Bin events by d_perp (e.g., 0–50m, 50–150m, 150–300m) and store empirical CDFs.
Simple, model-free, but requires many events per bin.

**Option C — Gaussian + exponential:**
```
p(Δt | d) = (1-f) · N(μ(d), σ(d)) + f · Exp(τ(d)) for Δt > 0
```
Fit μ, σ, τ, f as smooth functions of d_perp (polynomial or spline).

**Recommended starting point:** Option A (Pandel) — physically motivated and matches existing
IceCube infrastructure, so the DM-Ice term will be compatible with MPEFit/SplineMPE.

### Step 4: Detection efficiency (Poisson amplitude)

For each event, record whether DM-Ice fired at all vs. track geometry:
- Bin by d_perp and zenith
- Efficiency = (events with DM-Ice hit) / (total events in bin)
- Fit efficiency function: `ε(d_perp) = ε₀ · exp(-d_perp / λ_eff)`

This will also tell us the useful range of the model (d_perp < ~100–200 m probably).

### Step 5: Implement as IceTray module

```python
# ~/dmice/dmice_likelihood.py

class DMIceLikelihoodTerm(icetray.I3Module):
    """
    Adds log L_DM-Ice to a particle hypothesis's likelihood.
    
    Reads: frame[pulse_key] for DM-Ice OMKeys
           frame[track_key] for current track hypothesis
    Writes: frame["DMIce_LogL"] = I3Double(log_likelihood)
    """
    def __init__(self, context):
        super().__init__(context)
        self.AddParameter("TrackKey", "I3Particle with track hypothesis", "MPEFit")
        self.AddParameter("PulseKey", "Pulse map key", "InIcePulses")
        self.AddParameter("ModelFile", "Path to timing model (npz or pickle)", "")

    def Configure(self):
        self.track_key = self.GetParameter("TrackKey")
        self.pulse_key = self.GetParameter("PulseKey")
        model = np.load(self.GetParameter("ModelFile"))
        # Load fitted Pandel params: lambda_a, lambda_s as function of d_perp
        self.model = model

    def Physics(self, frame):
        track = frame[self.track_key]
        # 1. Get DM-Ice hit times from pulse map
        # 2. Compute t_geo from track hypothesis
        # 3. Evaluate Pandel PDF at Δt = t_obs - t_geo
        # 4. Add log likelihood
        frame["DMIce_LogL"] = dataclasses.I3Double(log_l)
        self.PushFrame(frame)
```

For true integration into the minimizer (not just post-hoc scoring), this needs to plug into
gulliver's likelihood service interface. This is the harder step — see below.

### Step 6: Gulliver likelihood service integration

For SplineMPE to use DM-Ice, implement `I3EventLogLikelihood` interface:
```python
# Inherits from: lilliput.I3EventLogLikelihood (C++ binding)
# Method: GetLogLikelihood(hypothesis) → float
# Called by: gulliver minimizer during fit
```

This allows MPEFit to be replaced by a combined fit that minimizes over (θ,φ,t₀,x₀,y₀,z₀)
using both IceCube DOM pulses AND DM-Ice timing simultaneously.

**Simpler alternative:** Run MPEFit first (standard IceCube fit), then apply DM-Ice correction
as a Bayesian update — valid if DM-Ice adds independent information orthogonal to IceCube.

---

## Implementation Roadmap

### Phase A: Model extraction (1–2 days)
**Script:** `~/dmice/build_dmice_timing_model.py`
- Load existing 200-event npz
- Extract DM-Ice hit times vs MC truth track geometry
- Compute Δt residuals and d_perp for each event
- Fit Pandel params (or Gaussian+exp)
- Output: `~/dmice_work/output/dmice_timing_model.npz` with fitted parameters + lookup tables
- Plot: residual distributions binned by d_perp

### Phase B: Scoring module (1 day)
**Script:** `~/dmice/dmice_likelihood.py`
- Load timing model
- Score existing MPEFit tracks against DM-Ice hits on real coincidence data
- Plot log L_DM-Ice distribution — does it discriminate real DM-Ice tracks?

### Phase C: Dedicated sim run (2–3 days on WARD)
- Run simulate_muons_binned.py with larger statistics (1000+ events)
- Wider d_perp coverage (inject tracks at various offset distances, not just through detector)
- Need off-axis events: inject tracks aimed 10m, 50m, 100m, 200m from DM-Ice
- This gives the efficiency vs. d_perp curve

### Phase D: Gulliver integration (1 week, ambitious)
- Implement I3EventLogLikelihood subclass
- Test against MPEFit-only baseline on BLO sim (ground truth available)
- Metric: median angular error on 200-event sim

---

## Critical Files

| File | Purpose |
|------|---------|
| `~/dmice/run_sim_all_recos.py` | Reference: IceTray NPZInjector, MPEFit setup |
| `~/dmice/run_pivot_mpefit.py` | Reference: pivot_linefit(), DOM pulse extraction |
| `~/dmice/BlueLightOrchestra.jl/resources/geofiles/icecube_with_dmice.geo` | DM-Ice DOM positions |
| `~/dmice_work/output/muons_binned_200ev_repacked.npz` | Input sim data for Phase A |
| `~/dmice_work/output/comparison/sim_all_recos.csv` | Current performance baseline |

## Reusable Functions

- `pivot_linefit_ic()` in `run_sim_all_recos.py:95` — DM-Ice transit time calculation (same geometry)
- `load_ragged()` in `run_sim_all_recos.py:74` — flat+offsets npz format reader
- `ang_err_deg()` in `run_sim_all_recos.py:114` — angular difference calculation
- `load_geo()` in `run_sim_all_recos.py:41` — parse geofile to get DOM positions

---

## Verification

1. **Residual plots**: Δt distributions should be Pandel-shaped (skewed, positive tail). If they
   look Gaussian-only, something is wrong with t_geo calculation.

2. **Sanity check**: For tracks that pass directly through DM-Ice (d_perp ≈ 0), Δt should peak
   near 0 with minimal scatter. For d_perp = 100m, peak should shift by ~100/(0.2998·sin41°) ≈ 508 ns.

3. **Benchmark**: After building scoring module, apply to 200-event BLO sim:
   - Score MPEFit (ground truth available) vs DM-Ice log L
   - Events with high DM-Ice L should have better MPEFit accuracy (already have DM-Ice info)
   - This validates the model captures real signal

4. **Real data sanity**: Apply scorer to real coincidence events (2012/2013/2018).
   Events with d_perp < 50m should have much higher L_DM-Ice than those with d_perp > 200m.
