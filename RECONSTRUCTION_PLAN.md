# DM-Ice Muon Track Reconstruction — Implementation Plan

**Project:** IceCube / DM-Ice NaI coincidence analysis
**Date:** 2026-04-15
**Status:** Active development

Related: [[INDEX]] | [[RESULTS]] | [[ALTERNATIVE_RECOS]] | [[READING_LIST]] | [[DMICE_TIMING_IMPLEMENTATION]]

---

## Executive Summary

The current pipeline achieves MPEFit 0.55° median angular error (simulation,
1000-event binned sample, ANALYSIS-04) and Pivot LineFit 1.71°. The combined
DM-Ice + IceCube likelihood fit (TABLED-01) diverged to ~60° because the
uniform-ice Pandel approximation does not match the PPC/SpiceMie photon
propagation model used in simulation.

This document lays out a phased plan: Phase A fixes the current pivot seeding
properly; Phase B builds a correct custom joint likelihood using SplineMPE tables
or a full Gulliver service; Phase C explores machine-learning track reconstruction
as a longer-term investment. All phases use the existing DM-Ice timing anchor:
Gaussian(μ=+280 ns, σ=81 ns) for real data; Gaussian(μ=−42 ns, σ=47 ns) for BLO sim.

The most important single fact about the physics: **DM-Ice fires on direct
ionisation/scintillation (d⊥ ≈ 0), not on Cherenkov photons.** The hit time is
therefore a near-exact position+time anchor on the muon track, smeared only by
NaI scintillation physics, not by ice propagation.

---

## 1. Existing LineFit Variants in IceCube

### 1.1 I3LineFit (standard PoleMuonLinefit)

The analytic weighted least-squares solution minimising

    χ² = Σᵢ wᵢ (rᵢ − r₀ − v tᵢ)²

in the time component, with wᵢ = qᵢ^AmpWeightPower (usually 1). The direction is
the unit vector along v, and the speed |v| is unconstrained (but physically
expected near c).

**IceTray call:**
```python
tray.Add("I3LineFit",
    Name            = "LineFit",
    InputRecoPulses = "InIcePulses",
    AmpWeightPower  = 1.0,
)
```

Current performance (ANALYSIS-04, 1000 events): 4.48° median.

### 1.2 Amplitude-Weighted LineFit

Already the default (AmpWeightPower=1.0). Increasing to 2.0 weights bright DOMs
more heavily. This is suboptimal for muons (bright DOMs are often near hadronic
showers, not on the Cherenkov cone), but worth a quick test. Expected gain: small.

### 1.3 Iterative Outlier Rejection LineFit

Runs LineFit, computes timing residuals for every hit, removes hits with |residual|
> N·σ_residual, reruns LineFit on the cleaned hit set. Repeats for a fixed number
of iterations. This is the "Improved LineFit" mentioned in ALTERNATIVE_RECOS.md.

In practice on cobalt, it is faster to implement iterative rejection in Python
directly on the npz event data (no IceTray overhead) for benchmarking.

**Expected gain as pivot seed:** reduces LineFit direction error from ~4.5° to
~2–3°, tightening the pivot anchor projected time.

### 1.4 Direct-Hit LineFit

Run LineFit only on "direct" hits — pulses with timing residuals within a tight
window (typically −15 to +75 ns) relative to the Cherenkov cone of a seed track.
These are predominantly unscattered photons and give a cleaner geometric direction.

**IceTray call:**
```python
from icecube import direct_hits
tray.AddModule("I3DirectHitsCalculator",
    DirectHitSeriesMapName    = "DirectHits",
    PulseSeriesMapName        = "InIcePulses",
    ParticleName              = "LineFit",
    DirectHitDefinitionSeries = [
        direct_hits.I3DirectHitsDefinition("A", -15, 75),
    ],
)
tray.Add("I3LineFit",
    Name            = "DirectLineFit",
    InputRecoPulses = "DirectHitsA",
    AmpWeightPower  = 1.0,
)
```

### 1.5 SPEFit as Intermediate Seed

SPEFit achieves ~1.01° median (ANALYSIS-04) and is computationally cheap.
Using SPEFit as the pivot seed — instead of LineFit — is the logical middle step:

```
LineFit → SPEFit (std) → Pivot anchor using SPEFit direction → MPEFit (pivot)
```

### 1.6 IterativePandelFitter (IterMPE)

Runs MPEFit, removes outlier pulses with large residuals, reruns. Can converge
to ~0.3° median for high-energy events. Recommended in ALTERNATIVE_RECOS.md.

**IceTray call:**
```python
tray.Add(icecube.lilliput.segments.I3IterativePandelFitter,
    fitname      = "IterMPE",
    domllh       = "MPE",
    pulses       = "InIcePulses",
    seeds        = ["LineFit"],
    n_iterations = 3,
)
```

Then pass "IterMPE" as the pivot seed direction. **This is the highest-priority
addition in Phase A.**

---

## 2. Custom Likelihood-Based Reconstruction Using the DM-Ice Constraint

### 2.1 Why the Current Combined Fit Failed (TABLED-01)

The `DMCombinedFitModule` in `run_sim_all_recos.py` uses analytic Pandel parameters
(`_PANDEL_LA=98 m`, `_PANDEL_LS=30 m`) for uniform ice. PPC simulation uses the
SpiceMie layered ice model, which has depth-dependent absorption and scattering
lengths differing by factors of 2–5 from the bulk averages. The optimizer
(Nelder-Mead on the combined objective) finds a direction where the Pandel misfit
is minimized at the expense of the DM-Ice Gaussian term, yielding ~60° median errors.

**The fix is not to tune the Pandel constants — it is to replace the IC likelihood
term with a table-based evaluation that matches the simulation.**

### 2.2 Correct Approach: SplineMPE Tables as the IC Likelihood

SplineMPE replaces the analytic Pandel PDF with bicubic spline interpolation of
pre-computed photon arrival time CDFs as a function of (d_perp, t_residual),
evaluated on the SpiceMie ice model. These tables exist on cobalt at:

```
/cvmfs/icecube.opensciencegrid.org/data/photon-tables/splines/
    ems_mie_z20_a10.abs.fits    # absorption spline
    ems_mie_z20_a10.prob.fits   # probability spline
```

**Architecture for a correct custom fit:**

```
log L_total = log L_SplineMPE(IC pulses | track) + log L_DM-Ice(t_obs | track, DM position)
```

where:
- `log L_SplineMPE` = SplineMPE multi-photon likelihood evaluated via IceTray's photospline service
- `log L_DM-Ice` = -0.5 * ((t_obs - t_geo(track, dm_pos) - μ) / σ)²

The minimiser optimises over (zenith, azimuth, x₀, y₀, z₀, t₀) — all 6 track
parameters — simultaneously.

### 2.3 Implementation via Gulliver Framework

The cleanest IceCube-native approach is to implement the DM-Ice term as a Gulliver
`I3EventLogLikelihood` service, then add it to the existing MPEFit minimization.

**DM-Ice Gulliver service sketch:**

```python
# dmice_gulliver_llh.py
from icecube import gulliver, icetray, dataclasses
import math

MU_NS    =  280.0   # real data (use -41.9 for BLO sim)
SIGMA_NS =   81.0
C_M_NS   = 0.2998
THETA_C  = math.acos(1.0 / 1.3195)

class DMIceLogLikelihood(icetray.I3Module):
    """
    Gulliver I3EventLogLikelihood service: DM-Ice NaI timing term.
    Adds log L = -0.5*((t_obs - t_geo - μ)/σ)² to any Gulliver fit.
    """
    def Configure(self):
        self.dm_t_key = self.GetParameter("DMIceTime")
        self.dm_pos   = self.GetParameter("DMIcePos")

    def SetEvent(self, frame):
        if self.dm_t_key not in frame:
            return False
        self.t_obs = frame[self.dm_t_key].value   # already μ-corrected
        return True

    def GetLogLikelihood(self, hypothesis):
        p = hypothesis.particle
        # project DM-Ice pos onto track, compute geometric time
        dx = math.sin(p.dir.zenith) * math.cos(p.dir.azimuth)
        dy = math.sin(p.dir.zenith) * math.sin(p.dir.azimuth)
        dz = -math.cos(p.dir.zenith)
        rx = self.dm_pos[0] - p.pos.x
        ry = self.dm_pos[1] - p.pos.y
        rz = self.dm_pos[2] - p.pos.z
        s  = rx*dx + ry*dy + rz*dz
        t_geo = p.time + s / C_M_NS  # d_perp ≈ 0 for DM-Ice (direct ionisation)
        return -0.5 * ((self.t_obs - t_geo) / SIGMA_NS)**2
```

**Combined service in tray:**
```python
tray.AddService("I3CombinedLogLikelihood", "CombinedLLH",
    LogLikelihoods = ["SplineMPELLH", "DMIceLLH"],
    Multipliers    = [1.0, 1.0],
)
```

Use Minuit (MIGRAD) as the minimizer, not Nelder-Mead. Minuit is faster and
more robust for smooth likelihood surfaces.

### 2.4 Multi-Seed Timing Uncertainty

The NaI scintillation delay has σ = 81 ns. Run three seeds:

```
Seeds: {dm_t_corrected − σ, dm_t_corrected, dm_t_corrected + σ}
→ Run Pivot MPEFit from each
→ Keep the result with highest IC log L
```

Adds ~3× compute but makes the result robust to single-photon fluctuations
(81 ns → ~24 m positional uncertainty along track at c).

---

## 3. Machine Learning Approaches

### 3.1 Context: ML Track Reconstruction in IceCube

IceCube has invested heavily in graph neural network (GNN) reconstruction since
~2019. The key published work is IceCube's DynEdge / GNN-Reco (graphnet) framework.
These methods treat detector hits as a point cloud and train end-to-end to predict
track direction from raw pulses.

**Published performance:**
- DeepCore (6–100 GeV): GNNs outperform SPEFit by ~30% at 10 GeV
- Muon tracks (>1 TeV): GNNs match or slightly beat SplineMPE at much lower inference cost

**DM-Ice novelty:** No existing IceCube ML model includes DM-Ice as an input feature.
Adding the DM-Ice hit time as a node in the hit graph is a genuinely novel contribution.

### 3.2 Graph Neural Network Architecture (DynEdge)

```
Input: per-hit features (x, y, z, t, charge, is_dmice) for each pulse
Graph: k-nearest neighbours in (x,y,z,t) space
Layers: 4–8 EdgeConv message-passing layers
Pooling: global mean/sum
Output: (sin θ cos φ, sin θ sin φ, cos θ) + uncertainty estimate (κ)
```

**Adding DM-Ice:** inject DM-Ice hit as extra node with `is_dmice=1.0` flag.
All IC DOM nodes get `is_dmice=0.0`. The network learns that this node carries
a special position+time anchor.

### 3.3 Training Data Requirements

Current dataset: 5000 events (SIM-04). GNN training needs ~10⁵–10⁶ events.

```bash
# On WARD:
python3 ~/dmice/simulate_muons_binned.py --n-per-bin 20000
# → 100k events, ~10 GPU-hours
```

### 3.4 Available Frameworks

**graphnet (recommended starting point):**
- Repository: `github.com/graphnet-team/graphnet`
- Contains: DynEdge model, IceCube dataset loaders, training scripts
- Install: `pip install --user graphnet`

**Key limitation:** BLO sim models DM-Ice as a pseudo-Cherenkov DOM (physically
wrong — real DM-Ice fires only at d_perp ≈ 0). GNN trained on BLO sim will have
systematic biases on real data. Domain adaptation or fine-tuning on the 586 real
pivot events is required for deployment.

### 3.5 Simpler Alternative: BDT on Reco Features

A 2-week (vs 2-month) ML approach: train a BDT/random forest on the outputs of
existing reconstructions:

```
Input features:
    LineFit + MPEFit directions, MPEFit log L, n_doms, total charge,
    DM-Ice timing residual (t_obs - t_geo_MPEFit) - μ,
    d_perp of MPEFit track at DM-Ice position, TruncatedEnergy

Output: corrected (zenith, azimuth)
```

Faster to implement, publication-ready, and interpretable. **Recommended before
the GNN approach.**

---

## 4. Reading List

These 10 references specifically target reconstruction development. See also
`READING_LIST.md` for broader project context (not duplicated here).

---

### R1. Ahrens et al. — "Muon track reconstruction and data selection techniques in AMANDA"
**arXiv:** astro-ph/0407044 | NIM A 524 (2004)

Derives the analytic LineFit formula (Eq. 2.1–2.5), explains why it is biased
(constant-velocity assumption, no scattering), then introduces SPEFit/MPEFit.
**Read this to understand exactly what `pivot_linefit_ic()` is doing and its limits.**

---

### R2. IceCube Collaboration — "Energy Reconstruction Methods in the IceCube Neutrino Telescope"
**arXiv:** 1311.4767 | JINST 9 P03009 (2014)

Covers SplineMPE in detail: how photospline tables are constructed, how they replace
the analytic Pandel PDF, and why SplineMPE outperforms MPEFit when the ice model is
correct. **Essential reading before implementing Phase B.**

---

### R3. IceCube Collaboration — "Measurement of South Pole ice transparency with the IceCube LED calibration system"
**arXiv:** 1301.5361 | NIM A 711 (2013)

Describes SpiceMie, the layered ice model used by PPC. Explains why bulk-ice Pandel
parameters break down (absorption/scattering vary by ~4× between layers). **Read
this to understand why TABLED-01 failed.**

---

### R4. Lundberg et al. — "Light tracking through ice and water — Scattering and absorption in heterogeneous media with Photonics"
**arXiv:** astro-ph/0702108 | NIM A 581 (2007)

Describes the physics of photon propagation in layered ice — the same physics the
SplineMPE tables tabulate. Important background for Phase B.

---

### R5. Pandel — "Bestimmung von Wasser- und Detektoreigenschaften..."
**Source:** Diploma thesis, Humboldt University Berlin (1996). Available via IceCube internal library.

Original derivation of the Pandel PDF. Explains why the PDF has a Gamma-function
form and how absorption and scattering lengths enter. Read the key approximation
(Eq. 3.4) to understand why TABLED-01 fails: bulk parameters are ~50% off from
true layered-ice values at the depth of DM-Ice (−511 m, the clearest ice layer).

---

### R6. Abbasi et al. — "A convolutional neural network based cascade reconstruction for IceCube"
**arXiv:** 2101.11589 | JINST 16 P07041 (2021)

Describes IceCube's CNN-based reconstruction — input feature engineering (raw DOM
hits → fixed tensor), training procedure, and comparison to classical methods.
**Read Sections 3–5** for data format and training strategy applicable to Phase C.

---

### R7. IceCube Collaboration — "Graph Neural Networks for Low-Energy Event Classification and Reconstruction in IceCube"
**arXiv:** 2209.03561 | JINST 17 P11003 (2022)

The published DynEdge / GNN-Reco paper. Describes the architecture, k-NN graph
construction, EdgeConv layers, and results on DeepCore events.
**This is the architecture to adapt for Phase C.**

---

### R8. Eller et al. — "Bringing IceCube's GNN-based reconstruction to the rest of the collaboration with graphnet"
**arXiv:** 2210.12194 | PoS ICRC2023 (2023)

Documents the open-source graphnet library with concrete training examples.
**The most practical starting point for Phase C.** Read the README and worked
examples before writing any GNN code.

---

### R9. Cowan — *Statistical Data Analysis in Particle Physics* (Oxford, 1998)
**Chapters:** 7 (Maximum likelihood estimation), 9 (Hypothesis testing)

The standard particle-physics statistics reference. Chapter 7 covers maximum
likelihood estimation and profile likelihood — directly applicable to Phase B's
angular uncertainty estimation. Chapter 9 covers the Wilks' theorem approach for
confidence intervals from the log-likelihood surface.

---

### R10. Barlow — *Statistics: A Guide to the Use of Statistical Methods in the Physical Sciences* (Wiley, 1989)
**Chapters:** 5 (Maximum likelihood), 6 (Least squares)

Clear pedagogical treatment of the statistical foundations of all reconstructions
in this project. Better starting point than Cowan if you are less familiar with
frequentist statistics. Read before attempting to implement Phase B uncertainty
estimates.

---

## 5. Phased Implementation Plan

### Phase A — Fix the Existing Pivot (1–2 weeks, Cobalt)

**Priority: highest. Do these in order.**

#### A1. Diagnose μ correction (in progress)
Two ablation tests running now on Cobalt (`reco_step1only`, `reco_step2only`).
Results will tell us whether the μ correction in `pivot_linefit_ic()` helps or
hurts, and whether the vertex anchor change is the source of degradation.

#### A2. Fix MPEFit seed t₀

In `compute_pivot_lf()`, after computing the pivot direction, anchor the MPEFit
seed vertex at the DM-Ice position with t₀ derived from the corrected hit time:
```python
pp.pos  = dataclasses.I3Position(dm_pos[0], dm_pos[1], dm_pos[2])
s       = dot(dm_pos - lf_particle.pos, piv_dir)
pp.time = dm_t_corrected - s / C_M_NS
```
Expected: Pivot MPEFit drops from 0.55° to ~0.45°.

#### A3. Add IterMPE as pivot seed

Add `I3IterativePandelFitter` (n_iterations=3) seeded from LineFit, then use its
direction as the pivot anchor. The IterMPE direction (~0.4°) is much cleaner than
the LineFit direction (~4.5°), giving a much more precise projected DM-Ice time.

#### A4. Multi-seed timing uncertainty

Run Pivot MPEFit from seeds at dm_t_corrected ± σ and dm_t_corrected. Keep the
result with highest SPE log L. Adds 3× compute but handles single-photon timing
fluctuations.

#### A5. Re-run real data pipeline after fixes

The current −11.5° zenith shift may be a μ correction artifact. Re-run
`run_all_recos_real.py` after A1–A3 are applied (MU_NS = +280.0 for real data).

**Phase A deliverable:** Updated benchmark table. Expected Pivot LineFit: ~1.0°.
Expected Pivot MPEFit: ~0.45°.

---

### Phase B — Correct Combined Likelihood (2–6 weeks)

#### B1. Implement DM-Ice Gulliver service
Write `dmice_gulliver_llh.py` as described in Section 2.3.
Test in isolation: score MC truth tracks and verify log L peaks at truth direction.

#### B2. Combined SplineMPE + DM-Ice fit via Gulliver
Register DM-Ice service alongside SplineMPE. Use Minuit minimizer.
Seed from Pivot LineFit (not standard LineFit).

#### B3. Benchmark Phase B vs Phase A

| Fit | Expected median |
|-----|-----------------|
| LineFit | 4.5° |
| Pivot LineFit (A, fixed) | ~1.0° |
| IterMPE (std seed) | ~0.4° |
| Pivot MPEFit (A, fixed) | ~0.45° |
| Combined SplineMPE+DM-Ice (B) | ~0.3–0.35°? |

#### B4. Profile likelihood angular uncertainty
Scan (θ, φ) while minimizing over (x₀, y₀, z₀, t₀) for each direction hypothesis.
68% confidence region where ΔlogL < 1.15 is the standard IceCube angular error.

**Phase B deliverable:** Combined fit < 0.4° median on simulation, with per-event
uncertainty estimates. Publication-quality benchmark table.

---

### Phase C — Machine Learning (2–4 months)

#### C1. Generate 100k-event training set on WARD
```bash
python3 ~/dmice/simulate_muons_binned.py --n-per-bin 20000
```

#### C2. Install graphnet on Cobalt
```bash
pip install --user graphnet
```

#### C3. Convert BLO npz to graphnet SQLite format
Write `npz_to_graphnet.py`: IC hits as nodes (is_dmice=0), DM-Ice hit as node
(is_dmice=1), labels = MC truth (zenith, azimuth).

#### C4. Train DynEdge with DM-Ice conditioning
Use `nb_inputs=8` (x, y, z, t, charge, string, dom, is_dmice), DM-Ice time as
global conditioning variable after pooling.

#### C5. BDT shortcut (recommended first)
Before the full GNN: train a BDT on existing reco output features (Section 3.5).
2 weeks vs 2 months, immediately publishable.

**Phase C deliverable:** GNN achieving < 0.5° median on BLO test set, with DM-Ice
node ablation study quantifying the contribution of the timing anchor.

---

## 6. Infrastructure Notes

### Compute Resources

| Task | Machine | Estimate |
|------|---------|----------|
| Phase A reco runs | Cobalt-14 | 1–2 hours / 1000-event batch |
| Phase A real data re-run | Cobalt or NPX Condor | 2–4 hours (4851 events) |
| Phase B combined fit | Cobalt-14 | 3–5× slower than MPEFit alone |
| Phase C: 100k sim | WARD (GPU) | ~10 GPU-hours |
| Phase C: GNN training | Cobalt GPU nodes (Condor) | ~4–8 GPU-hours / run |

### SplineMPE Tables Location (Cobalt)
```
/cvmfs/icecube.opensciencegrid.org/data/photon-tables/splines/
    ems_mie_z20_a10.abs.fits
    ems_mie_z20_a10.prob.fits
```

### Validation Protocol for Every New Reco

1. Check valid event counts — not NaN-filtering away half the sample
2. Check per-bin medians — performance should scale with energy as expected
3. Verify μ correction sign: residual should centre near zero, not ±280 ns
4. Plot the angular error CDF, not just the median — tail behaviour matters
5. For real data: check zenith distribution is physically reasonable (downgoing muons)

### Known Issues to Monitor

| Issue | Status |
|-------|--------|
| μ correction effect on Pivot LF | Under investigation (ablation tests running) |
| Real data −11.5° zenith shift | Likely μ bug artifact; re-run after fix |
| BLO ↔ real-data timing mismatch (−42 ns vs +280 ns) | Always use correct μ per data source |
| Combined fit diverges (TABLED-01) | Phase B resolves completely |

---

## 7. Publication Strategy

**ICRC 2026 proceedings** (~July 2026 deadline): Phase A alone is sufficient.

**JINST/NIM methods paper**: requires Phase B + data/MC comparison.

**PRD/JCAP physics paper**: requires DM flux limit or astrophysical result.

**Hardware argument (PAPER_IDEAS.md):** simulate N = 1, 2, 4, 8 DM-Ice detectors,
compute fraction of events with ≥2 hits, show how two simultaneous hits fully
constrain track direction independently of IceCube. Clean self-contained result
that motivates the detector upgrade case.

---

## 8. Priority Queue

**This week:**
1. Interpret ablation test results (μ correction vs vertex anchor)
2. Apply whichever fix helps; re-run 1000-event benchmark
3. Add IterMPE as pivot seed
4. Re-run real data pipeline with corrected μ

**Next 2–3 weeks:**
5. Implement DM-Ice Gulliver service (Phase B1)
6. Test combined SplineMPE + DM-Ice on 200-event sim
7. Full 5000-event benchmark on combined fit

**1–3 months:**
8. BDT on reco features (quick ML win)
9. Generate 100k training set; set up graphnet
10. Draft ICRC proceedings

---

*Supersedes scattered notes in TOMORROW_PLAN.md, ALTERNATIVE_RECOS.md, and
DMICE_TIMING_IMPLEMENTATION.md for reconstruction-development planning. Those
files remain authoritative for their specific topics.*
