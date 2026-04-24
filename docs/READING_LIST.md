# Suggested Reading List — DM-Ice Pivot Reconstruction Project

Organised by topic. Chapters noted where a full textbook is cited.

Related: [[INDEX]] | [[RECONSTRUCTION_PLAN]] | [[DMICE_TIMING_IMPLEMENTATION]]

---

## 1. NaI Scintillator Physics

### 1.1 Knoll — *Radiation Detection and Measurement* (4th ed., Wiley, 2010)
**Chapters:** 8 (Scintillation detector principles), 10 (Photomultiplier tubes)

The standard detector physics textbook. Chapter 8 covers NaI(Tl) in detail: scintillation
mechanism, ~250 ns decay constant, light yield, energy resolution, and why the timing
response is an exponential decay convolved with PMT transit-time jitter (not a Gaussian).
Essential for understanding why our Gaussian(280 ns, 81 ns) is an approximation and where
it breaks down.

---

### 1.2 DM-Ice17 — *First search for dark matter with DM-Ice17*
**arXiv:** 1612.05093 | **Journal:** Phys. Rev. D 95, 032006 (2017)

The primary DM-Ice17 paper. Contains the detector description, NaI crystal timing
characterisation, and coincidence analysis methods that underpin the 280 ns / 81 ns
timing model used in this project.

---

### 1.3 DM-Ice17 — *Low-background NaI(Tl) detectors deployed in the deep ice at the South Pole*
**arXiv:** 1602.05765 | **Journal:** Astropart. Phys. 84 (2016)

Hardware paper. Details the PMT readout, crystal housing, deployment geometry, and
calibration procedure. Read alongside 1612.05093 to understand where the timing
parameters come from instrumentally.

---

## 2. IceCube Muon Reconstruction

### 2.1 Ahrens et al. — *Muon track reconstruction and data selection techniques in AMANDA*
**arXiv:** astro-ph/0407044 | **Journal:** NIM A 524 (2004)

Foundational paper for LineFit and the single/multi-photoelectron Pandel likelihood
(SPEFit, MPEFit) still used in IceCube today. Derives the analytical LineFit solution
and explains why it is fast but biased, and how the likelihood fits improve on it.
Read this before trying to modify or extend the reconstruction.

---

### 2.2 Pandel — *Bestimmung von Wasser- und Detektoreigenschaften und Rekonstruktion von Myonen bis 100 TeV mit dem Baikal-Neutrinoteleskop NT-72*
**Source:** Diploma thesis, Humboldt University Berlin (1996)

The original derivation of the Pandel PDF for photon arrival time distributions in
natural water/ice. The approximation used in SPEFit/MPEFit (`_PANDEL_LA=98m`,
`_PANDEL_LS=30m`) comes from this work. Hard to find but worth knowing it exists.

---

### 2.3 IceCube Collaboration — *Energy Reconstruction Methods in the IceCube Neutrino Telescope*
**arXiv:** 1311.4767 | **Journal:** JINST 9 P03009 (2014)

Covers SplineMPE, TruncatedEnergy, and MuEx — the reconstruction chain used in
IceCube analyses. Explains how spline tables replace the analytic Pandel PDF to
account for layered ice (SpiceMie model), which is the correct replacement for the
approximate uniform-ice Pandel used in this project's SPEFit/MPEFit.

---

### 2.4 IceCube Collaboration — *The IceCube Neutrino Observatory: Instrumentation and Online Systems*
**arXiv:** 1612.05093 wait — use: **arXiv:** 1612.05093 is DM-Ice. Use:
**arXiv:** 0602115 (IceCube detector JINST paper) | **Journal:** JINST 1 P11003 (2006)

Describes the IceCube detector geometry, DOM hardware, digitisation, and the
SMT8 Simple Majority Trigger used as the primary online trigger. Section on
triggering is directly relevant to the SMT8 efficiency study in this project.

---

### 2.5 IceCube Collaboration — *Measurement of South Pole ice transparency with the IceCube LED calibration system*
**arXiv:** 1301.5361 | **Journal:** NIM A 711 (2013)

Describes the SpiceMie layered ice model used by PPC for photon propagation
simulation. Important for understanding why the uniform-ice Pandel approximation
diverges from PPC-simulated data at large distances or steep angles.

---

## 3. Dark Matter Detection

### 3.1 Bertone, Hooper & Silk — *Particle Dark Matter: Evidence, Candidates and Constraints*
**arXiv:** hep-ph/0404175 | **Journal:** Phys. Rept. 405 (2005)

Comprehensive review of WIMP dark matter theory, detection methods (direct, indirect,
collider), and the landscape of candidates. Read Sections 4–5 for direct detection
and the connection to NaI-based experiments like DAMA/LIBRA and DM-Ice.

---

### 3.2 Lewin & Smith — *Review of mathematics, numerical factors, and corrections for dark matter experiments based on elastic nuclear recoil*
**Journal:** Astropart. Phys. 6, 87 (1996)

The standard reference for WIMP-nucleus scattering rate calculations — form factors,
velocity distributions, annual modulation signal. Required reading if the project
eventually moves toward a DM flux limit.

---

### 3.3 DAMA/LIBRA Collaboration — *New results from DAMA/LIBRA-phase2*
**arXiv:** 1805.10486 | **Journal:** Universe 4, 116 (2018)

The claimed annual modulation signal that DM-Ice was designed to test. Understanding
what DM-Ice is looking for (and why DAMA's result remains controversial) provides
the physics motivation for the entire project.

---

## 4. Statistical Methods

### 4.1 Barlow — *Statistics: A Guide to the Use of Statistical Methods in the Physical Sciences*
(Wiley, 1989) — **Chapters:** 5 (Maximum likelihood), 6 (Least squares)

Clear, concise introduction to the statistical foundations of the likelihood-based
reconstructions in this project. Good starting point before the more specialised
particle-physics treatments.

---

### 4.2 Cowan — *Statistical Data Analysis in Particle Physics*
(Oxford, 1998) — **Chapters:** 7 (Parameter estimation), 9 (Hypothesis testing)

The standard particle-physics statistics reference. Chapter 7 covers maximum likelihood
estimation and its properties; directly applicable to understanding SPEFit/MPEFit
convergence behaviour and why low-energy events (few hits) fail to converge.

---

## 5. Track Reconstruction / Timing Anchors

### 5.1 IceCube Collaboration — *Characteristics of the diffuse astrophysical muon neutrino flux*
**arXiv:** 1607.08006 | **Journal:** Phys. Rev. Lett. 117, 241101 (2016)

Shows SplineMPE angular resolution in the context of a physics result (~0.3–0.5° for
high-energy throughgoing muons). Provides the performance benchmark that the DM-Ice
pivot reconstruction is approaching (~0.55° median in simulation at 1 TeV+).

---

### 5.2 Abbasi et al. (IceCube) — *A muon track reconstruction algorithm for the IceCube neutrino telescope*
**arXiv:** 2308.05388 | **Journal:** submitted (2023)

Recent IceCube paper describing modern muon reconstruction improvements. Relevant for
understanding the state of the art that the DM-Ice pivot method should be compared to
in any publication.

---

## Reading Order Suggestion

For someone new to the project:

1. Knoll Ch. 8 → understand NaI timing physics
2. 1602.05765 + 1612.05093 → understand the DM-Ice detector
3. astro-ph/0407044 → understand LineFit, SPEFit, MPEFit
4. 1311.4767 → understand SplineMPE and why it beats Pandel
5. Bertone et al. + DAMA/LIBRA paper → understand the physics motivation
6. Cowan Ch. 7 → understand the statistics of the fits





    "Results from the DM-Ice17 dark matter
    experiment", Phys. Rev. D — this is your
    experiment's own paper. Describes the NaI
     detector response and timing.

    ---
    Machine learning for IceCube
    reconstruction (the fallback)

    - Abbasi et al. (2021), "A convolutional
    neural network based cascade
    reconstruction for the IceCube Neutrino
    Observatory", JINST — IceCube's own ML
    reco paper. Shows what features matter
    and how performance compares to
    likelihood methods.
    - IceCube's DNN_reco and retro GitHub
    repos — practical examples of feature
    engineering for track reconstruction.

    ---
    Gulliver framework (IceTray-specific)

    - The IceCube internal wiki page on
    Gulliver is the best reference, but
    requires IceCube credentials. The next
    best thing is reading the source of
    run_splinempe_pivot.py in your own
    project alongside the lilliput and
    gulliver_modules IceTray documentation —
    the existing code is itself a working
    tutorial.

    ---
    Suggested order if time is short: Ahrens
    2004 → Cowan chapters 6 & 9 → Knoll
    chapter 10 → Bradascio 2019 → Bergstra &
    Bengio 2012.


