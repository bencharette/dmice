# DMice Papers Index

Quick reference for all papers and PDFs used in the DMice/IceCube coincidence project.
Claude Code: read this file first to find which paper is relevant before loading a full PDF.

---

## How to use this index

- **Full PDF**: load when you need exact quotes, figures, or detailed methodology
- **Notes file**: load for quick reference — pre-extracted key sections
- **Relevance tags**: scan these to find the right paper fast

---

## Papers

### Simulation & Software

| File | Notes | Relevance |
|------|-------|-----------|
| `prometheus_paper.pdf` | `prometheus_notes.md` | Prometheus framework, LI config, PPC photon propagation |
| `lepton_injector_paper.pdf` | `lepton_injector_notes.md` | LI API, injection modes (volume vs ranged), earth model |
| `proposal_paper.pdf` | `proposal_notes.md` | PROPOSAL muon propagation, energy loss |

### IceCube Detector & Reconstruction

| File | Notes | Relevance |
|------|-------|-----------|
| `icecube_detector.pdf` | `icecube_detector_notes.md` | IC86 geometry, DOM positions, coordinate system |
| `linefit_paper.pdf` | `linefit_notes.md` | LineFit algorithm, PoleMuonLinefit, direction conventions |
| `icecube_muon_reco.pdf` | `icecube_muon_reco_notes.md` | MPEFit, track reconstruction, angular resolution |

### DM-Ice Detector

| File | Notes | Relevance |
|------|-------|-----------|
| `dmice_detector.pdf` | `dmice_detector_notes.md` | DAQ timing (0.1 ns units), detector positions, scintillation response |
| `dmice_deployment.pdf` | `dmice_deployment_notes.md` | Physical installation, depth, IceCube coordinate offset (1948.07 m) |

### Dark Matter Physics

| File | Notes | Relevance |
|------|-------|-----------|
| `dm_search_icecube_2021.pdf` | `dm_search_notes.md` | DM search methodology, energy range, background estimation |
| `wimp_sensitivity.pdf` | `wimp_sensitivity_notes.md` | WIMP cross-section limits, slow particle signatures |

### Coincidence Analysis

| File | Notes | Relevance |
|------|-------|-----------|
| `coincidence_method.pdf` | `coincidence_notes.md` | Coincidence window (-10 to +60 µs), timing residuals, IQR method |

---

## Notes Template

When extracting notes from a new paper, use this structure.
Save as `~/dmice/papers/<shortname>_notes.md`.

```markdown
# [Paper Short Name] — Notes
**Full title:** ...
**Authors:** ...
**Year:** ...
**File:** ~/dmice/papers/<filename>.pdf

## Key facts (copy these into MEMORY.md if critical)
- ...

## Relevant sections
### Section X.X — [title]
[2-3 sentence summary of what this section says and why it matters for DMice]

### Section X.X — [title]
...

## Equations / values to remember
- ...

## Caveats / things that don't apply to our setup
- ...

## Quotes (exact, with page number)
- p.X: "..."
```

---

## Quick Lookup: Common Questions → Paper

| Question | Paper |
|----------|-------|
| What injection mode should I use? | `lepton_injector_notes.md` |
| What is the correct coincidence window? | `coincidence_notes.md` |
| Why is lf.time ≈ 12,000 ns? | `icecube_detector_notes.md` |
| What are the DM-Ice DAQ timing units? | `dmice_detector_notes.md` |
| How is the Z offset (1948.07 m) derived? | `dmice_deployment_notes.md` |
| What energy range is physically motivated? | `dm_search_notes.md` |
| How does PPC handle DM-Ice pseudo-DOMs? | `prometheus_notes.md` |

---

## Adding a new paper

1. Copy PDF to `~/dmice/papers/`
2. Tell Claude Code:
   ```
   Read ~/dmice/papers/<filename>.pdf and extract key sections relevant to 
   [your topic]. Save notes to ~/dmice/papers/<shortname>_notes.md using 
   the template in ~/dmice/papers/index.md
   ```
3. Add a row to the index table above
4. Add a row to the Quick Lookup table if it answers a common question
5. If it changes a critical convention, update `~/dmice/memory/local.md`
