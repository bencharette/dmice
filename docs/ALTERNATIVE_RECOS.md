# Alternative Reconstructions for DM-Ice Pivot Seeding

## The Core Idea

The DM-Ice pivot works by anchoring the track's time coordinate at the NaI crystal.
The accuracy of the pivot depends on how well the **seed direction** estimates the
track's true direction — a better seed means a more accurate d_perp and t_geo
calculation at DM-Ice, which means a tighter timing anchor.

Current pipeline (LineFit seed):
```
LineFit (5.75° median) → Pivot anchor → MPEFit pivot (0.39°)
```

The question: what if we used a better direction as the pivot seed?

---

## Currently Implemented

| Method | Median error (5000-ev sim) | Notes |
|--------|---------------------------|-------|
| LineFit | 5.75° | Standard, fastest |
| Pivot LineFit | 2.09° | DM-Ice anchored LineFit |
| MPEFit (std) | 0.51° | Pandel MPE, LineFit seed |
| **MPEFit (pivot)** | **0.39°** | Pandel MPE, Pivot LF seed |
| SPEFit (std) | 1.05° | Pandel SPE1st, LineFit seed |
| SPEFit (pivot) | 0.93° | Pandel SPE1st, Pivot LF seed |

---

## Most Promising: MPEFit-Seeded Pivot (now implemented)

**What it does:**
```
LineFit → MPEFit(std) → Pivot anchor using MPEFit direction → MPEFit(pivot2)
```

**Why it's better:**
The LineFit-seeded pivot uses a 5.75° direction to project DM-Ice onto the track.
MPEFit direction is ~10× more accurate (0.51°), so the computed t_geo at DM-Ice
is much closer to the true photon arrival time. The pivot anchor is therefore
tighter, and MPEFit(pivot2) starts from a better position.

**Expected gain:** The LineFit→pivot improvement was ~0.12° (0.51°→0.39°).
MPEFit-seeded pivot should be larger since the seed quality jump is much bigger.

**Implemented as:** `RealMPE_Piv2` in `run_all_recos_real.py`

---

## Other Options to Consider

### 1. I3ImprovLineFit (Improved LineFit)

**IceTray name:** `I3ImprovLineFit`

**What it does:** Standard LineFit with two improvements:
- Amplitude weighting: DOMs with more charge get higher weight
- Iterative outlier rejection: removes hits that are far from the fit, then refits

**Expected accuracy:** ~3–4° (better than 5.75° from vanilla LineFit)

**Why it helps the pivot:** Better seed direction → better t_geo at DM-Ice → better anchor.
Cheaper than MPEFit so useful if compute time is a concern.

**How to add:**
```python
from icecube import linefit
tray.AddModule("I3LineFit",
    Name         = "ImprovLineFit",
    InputRecoPulses = IC_PULSES,
    StoragePolicy = "OnlyBestFit",
    MinHits      = 4,
    AmpWeightPower = 1.0,   # amplitude weighting
)
```

---

### 2. IterativePandelFit (Iterative MPE)

**IceTray name:** `I3IterativePandelFitter` via `lilliput.segments`

**What it does:** Runs MPEFit multiple times. After each iteration, removes pulses
with large timing residuals (outliers), then refits. Converges to a cleaner solution.

**Expected accuracy:** ~0.3° (better than single-pass MPEFit 0.51°)

**Why it helps:** Even cleaner direction estimate for the pivot seed. Also more
robust against scattering tails in the residual distribution.

**How to add:**
```python
tray.Add(icecube.lilliput.segments.I3IterativePandelFitter,
    fitname     = "IterMPE",
    domllh      = "MPE",
    pulses      = IC_PULSES,
    seeds       = ["LineFit"],
    n_iterations = 3,
)
```

Then use `IterMPE` as the pivot seed direction.

---

### 3. DipoleFit

**IceTray name:** `I3DipoleFit`

**What it does:** Computes the direction from the charge-weighted dipole moment
of the hit DOM positions. No timing used — purely geometric from charge distribution.

**Expected accuracy:** ~5–8° (similar to LineFit for throughgoing muons)

**Why it might help:** Independent of timing, so not correlated with LineFit errors.
Could be useful as a cross-check or for events where timing is unreliable.

**Note:** DipoleFit has a 180° ambiguity (gives track axis not direction sense).
Needs to be disambiguated using LineFit or another directional fit.

**How to add:**
```python
tray.AddModule("I3DipoleFit",
    Name            = "DipoleFit",
    InputRecoPulses = IC_PULSES,
    MinHits         = 5,
    AmpWeightPower  = 1,
)
```

---

### 4. Tensor of Inertia (ToI)

**IceTray name:** `I3TensorOfInertia`

**What it does:** Finds the principal axis of the hit charge distribution using
the inertia tensor. The smallest eigenvalue axis = track direction.

**Expected accuracy:** ~5–10° (coarser than LineFit for muons)

**Why it might help:** Very robust — works even when timing is noisy or missing.
Good for understanding the event topology before timing-based fits.

**Note:** Also has 180° ambiguity. Primarily useful as a cross-check.

**How to add:**
```python
tray.AddModule("I3TensorOfInertia",
    Name            = "ToI",
    InputRecoPulses = IC_PULSES,
    MinHits         = 5,
    AmplitudeOption = 1,
)
```

---

### 5. DirectWalk + LineFit (direct photon seed)

**What it does:** First runs `I3DirectHitsValues` to identify "direct" photons
(small timing residuals, most likely unscattered). Then runs LineFit on only
those direct-hit DOMs. Cleaner direction estimate than using all pulses.

**Expected accuracy:** ~2–3° (better than full-pulse LineFit, closer to MPEFit)

**Why it helps:** Direct photons are geometrically closer to the true Cherenkov cone,
so the LineFit on direct hits is a better proxy for the true track direction.
Faster than MPEFit, could be used as an intermediate pivot seed.

**How to add:**
```python
from icecube import direct_hits
tray.AddModule("I3DirectHitsCalculator",
    DirectHitSeriesMapName  = "DirectHits",
    PulseSeriesMapName      = IC_PULSES,
    ParticleName            = "LineFit",
    DirectHitDefinitionSeries = [
        direct_hits.I3DirectHitsDefinition("A", -15, 75),   # tight
    ],
)
# Then LineFit on DirectHits["A"] only
```

---

## Recommended Testing Order

1. **MPEFit-seeded pivot** ← already implemented, results from current run
2. **IterativePandelFit seed** — likely the next biggest gain
3. **I3ImprovLineFit seed** — faster alternative to IterMPE if compute is limited
4. DirectWalk + LineFit — useful if IterMPE isn't available on the cluster

DipoleFit and ToI are lower priority for this use case — LineFit and MPEFit
already give a better directional seed for throughgoing muons.

---

## Key Files

| File | Purpose |
|------|---------|
| `run_all_recos_real.py` | Real data pipeline (MPEFit-pivot2 now included) |
| `run_splinempe_pivot.py` | BLO sim benchmark (add IterMPE here next) |
| `DMICE_TIMING_IMPLEMENTATION.md` | Timing model and known issues |
