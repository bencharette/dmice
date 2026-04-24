# DM-Ice Paper Ideas and Discussion

Related: [[INDEX]] | [[RESULTS]] | [[RECONSTRUCTION_PLAN]] | [[READING_LIST]]

## Core Result

Improved muon angular reconstruction using DM-Ice NaI scintillator timing as a
track anchor. The "pivot" method seeds IceCube reconstructions (LineFit, MPEFit,
SPEFit) from the DM-Ice detection time, improving angular resolution.

| Method              | Median error (5000-ev BLO sim) |
|---------------------|-------------------------------|
| LineFit (std)       | ~5.75°                        |
| MPEFit (std)        | ~0.51°                        |
| MPEFit (LF-pivot)   | ~0.39°                        |
| MPEFit (MPE-pivot)  | TBD (benchmark running)       |

Real data (2012–2019, 4611 events, 585 with pivot): angular shift plots pending.

---

## Key Physical Insight

**DM-Ice NaI crystals only fire when a muon passes directly through them.**
This means d_perp ≈ 0 for every coincident event — the crystal detects direct
ionisation/scintillation, not Cherenkov photons propagated through ice.

Consequences:
- The DM-Ice hit time directly gives t_geo at the detector position (d_perp = 0)
- The only smearing is NaI scintillation + electronics: Gaussian(μ=+280ns, σ=81ns) on real data
- The timing model is already fully characterised from real data — no Pandel or photon
  propagation model needed
- The pivot gives a near-exact position+time anchor on the muon track

---

## Simulation vs Real Data Limitation

The BLO simulation treats DM-Ice as a Cherenkov DOM (via PPC photon propagation).
Real DM-Ice detects direct ionisation. This causes a mismatch:

| Quantity | BLO sim      | Real data    |
|----------|--------------|--------------|
| μ offset | −42 ns       | +280 ns      |
| σ spread | 47 ns        | 81 ns        |

**Implication:** Sim angular resolution numbers are slightly optimistic for the
pivot reconstructions (tighter timing constraint than reality). IceCube-only
reconstructions (MPEFit std, SPEFit std, LineFit) are unaffected — PPC correctly
models Cherenkov detection in IceCube DOMs.

**Paper strategy:**
- Use sim for IceCube reconstruction comparison (fully valid)
- Use real data angular shift plots as the primary pivot result (no MC truth needed)
- Sim motivates the method; real data demonstrates the effect
- The σ mismatch is a footnote, not a fundamental problem

---

## Argument for More DM-Ice Detectors

**Current limitation:** Two small crystals cover a tiny solid angle → low coincidence
rate (~10% of subruns have any coincident events).

**What more detectors would give:**
- Higher coincidence rate → more events with improved reconstruction
- Multiple simultaneous DM-Ice hits per muon → stronger track constraint
- Two hits at different positions constrains direction independently of IceCube
  (constrains all 6 track parameters, not just 4)
- Denser array → improved reconstruction applies to a larger fraction of events,
  including neutrino-induced muons

**Physics argument:**
One DM-Ice hit constrains 4 of 6 track parameters (x, y, z + arrival time).
Two hits at separated positions fully constrain direction — DM-Ice becomes an
independent directional detector. This is qualitatively different from the
current single-anchor approach.

**Quantifiable for a paper:** Simulate N detectors at various positions, compute
fraction of throughgoing muons with ≥2 hits, and angular resolution improvement
vs N. Clean self-contained result.

---

## Connection to Broader Astrophysics

Muons detected by IceCube are atmospheric (produced in Earth's atmosphere by
cosmic ray interactions — they cannot travel from intergalactic space).

However, **neutrino-induced muons** DO have extragalactic origins (AGN, GRBs,
galaxy clusters). These are also track events in IceCube. Better angular
resolution for muon tracks directly improves:
- Neutrino point source identification
- Neutrino-cosmic ray source correlation studies

IceCube's best reconstruction (SplineMPE) achieves ~0.3–0.4° for high-energy
throughgoing muons. The DM-Ice pivot achieves ~0.39° — competitive.

The pivot principle (additional timing anchors improve angular resolution) is
general. A denser NaI array inside IceCube would apply this to neutrino-induced
muons, potentially improving extragalactic source identification.

---

## DM Flux Limit (Future Work)

A DM flux limit is an upper bound on the dark matter interaction rate set when
no statistically significant excess is observed above background.

Expressed as an exclusion curve: DM-nucleon cross-section (cm²) vs DM mass (GeV).
Any point above the curve is ruled out by the experiment.

DM-Ice was designed to detect WIMPs scattering off NaI nuclei (sodium/iodine)
via nuclear recoils — the same signal DAMA/LIBRA claims. A flux limit requires:

1. Expected DM signal rate (requires assuming a DM model + detector efficiency)
2. Expected background rate (atmospheric muons, noise)
3. Observed event rate
4. Detection efficiency

**This is not yet part of the current analysis.** The current work improves
muon reconstruction; a flux limit requires a signal model and background
estimation on top of that. Discuss with advisor which DM model to target.

---

## Publication Pathways

| Venue | Difficulty | What's needed beyond current work |
|-------|-----------|-----------------------------------|
| ICRC proceedings | Low | Real data angular shift plots (done after job finishes) |
| IceCube tech note | Low | Write-up, internal review |
| JINST / NIM (methods) | Medium | Systematic uncertainties, data/MC comparison |
| PRD / JCAP (physics) | High | DM flux limit or astrophysical result |

**Nearest-term target: ICRC 2025** — current work as-is would be accepted.

**What's still missing for a full physics paper:**
- 2020–2021 data (step1 jobs running, ~6 more months of coincidences)
- Systematic uncertainties (timing calibration, ice model, energy)
- A physics claim (DM flux limit or neutrino source constraint)
- IceCube collaboration internal review (~months)

---

## Current Job Status (as of 2026-04-10)

| Job | Cluster | Status | Output |
|-----|---------|--------|--------|
| BLO benchmark (MPE-pivot) | 13046844 | Running | `splinempe_pivot_comparison.csv` |
| Real data all recos | 13072231 | **Completed** | `real_all_recos.csv` (4611 events, 585 with pivot, 2012–2019) |
| 2020–2021 step1 | multiple | Running (31,506 jobs) | step1_muons/2020, /2021 |

**Next steps:**
1. Plot real data results: `python plot_real_recos.py`
2. When benchmark finishes: `python replot_benchmark.py`
3. When 2020–2021 step1 done: run step2 → step3 → merge
4. Discuss paper target with advisor

---

## Key Files

| File | Purpose |
|------|---------|
| `run_splinempe_pivot.py` | BLO sim benchmark (7 reco methods) |
| `run_all_recos_real.py` | Real data pipeline (reads fixed 58MB merged file) |
| `plot_real_recos.py` | Real data angular shift plots |
| `replot_benchmark.py` | BLO benchmark plots from CSV |
| `ALTERNATIVE_RECOS.md` | Other reconstruction methods to consider |
| `DMICE_TIMING_IMPLEMENTATION.md` | Timing model and known issues |
