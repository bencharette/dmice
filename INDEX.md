# DM-Ice Pivot Reconstruction — Knowledge Map

Central index for the project. Open in Obsidian and use **Graph View** (Ctrl+G) to see connections.

---

## Physics

- [[DMICE_TIMING_IMPLEMENTATION]] — NaI scintillation model, d_perp≡0, Gaussian(280,81 ns)
- [[READING_LIST]] — Key papers and textbooks
- [[dmice_photon_model_plan]] — BLO photon model notes

## Reconstruction Methods

- [[RECONSTRUCTION_PLAN]] — Phased development plan (Phase A/B/C)
- [[ALTERNATIVE_RECOS]] — Other reco approaches considered
- [[RESULTS]] — Run ledger: SIM-xx, ANALYSIS-xx, REAL-xx

## Analysis Results

| Method | Sim median | Notes |
|--------|-----------|-------|
| LineFit | 5.75° | Baseline |
| Pivot LineFit | 1.78° | DM spatial+time anchor |
| MPEFit std | 0.51° | Pandel likelihood |
| MPEFit pivot seed | 0.48° | Best full-sample reco |
| MPEFit spatial seed | 0.28° | Best reco (all-events fix pending) |

## Pipeline & Scripts

- [[SCRIPTS]] — All scripts, status flags, machine assignments
- [[COMMANDS]] — Command reference and workflow
- [[BLO_SIM_PIPELINE]] — WARD→Cobalt simulation pipeline

## Paper

- [[PAPER_IDEAS]] — Publication pathways, narrative, figures
- [[RESULTS]] — REAL-02: 586 pivot events, −11.5° zenith shift on real data

## Machines

| Machine | Role |
|---------|------|
| LOCAL | Dev, editing |
| NPX | HTCondor jobs |
| COBALT (`cobalt-14`) | IceTray runs |
| WARD | GPU BLO simulation |
