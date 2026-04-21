# Tomorrow Plan 2
Created: 2026-04-20

## Context: Where We Are

The core DM-Ice reconstruction work (sim benchmarks, IterMPE charge-cap fix, pivot
reconstructions) is done and in good shape. Today's session pivoted to the
**coincidence selection problem**: how do we identify the ~4% of genuine muon-DM-Ice
hits among ~96% radioactivity accidentals in the 6000-event master dataset?

## What We Tried Today

### Event Display Plots
- Script: `~/dmice/plot_event_display_coinc.py`
- Ran on cobalt, synced back to `~/dmice_work/output/event_displays/`
- Key finding: FAIL events at Δt=527 ns and Δt=604 ns look **geometrically genuine**
  (LineFit track threads directly through DM-Ice crystal) but fail the 3σ timing cut
  by only 4 ns. These are likely real muon-DM-Ice hits clipped by MPEFit vertex-time bias.

### Three Discriminants Compared (compare_coinc_cuts_v2.py)
Run on 2012 (731 events with LineFit):

| Method | Pass | Notes |
|--------|------|-------|
| MPEFit Gaussian \|Δt−280\|<243 ns | 5 | Too few — MPEFit vertex-time bias |
| First-hit anchor \|Δt−280\|<243 ns | 26 | Better; 3.6% pass rate |
| LF vs PivotLF Δθ < 5° | 55 | Contaminated — see below |
| LF vs PivotLF Δθ < 10° | 140 | Too loose |
| First-hit AND Δθ < 10° | 6 | Conservative combined |

### Why Angular Cut Alone Fails
The PivotLineFit **degenerates** for large dm_t_corrected: when the DM-Ice timing is
far from IC hit times, the velocity formula collapses to pointing from the crystal
toward the IC DOM charge-centroid — which is approximately the same as LineFit.
So accidentals with large dm_t also get small Δθ(LF, PivotLF), making the angular
cut non-discriminating by itself.

### ML Classifier Prototype (ml_coinc_classifier.py)
- Features: 10 pairwise angular diffs between all reco pairs + scalar IC properties
- Models: IsolationForest + GradientBoosting
- **Problem:** n_hits_ic dominated (53% importance) because the label proxy
  (dm_t > 40 μs as background) correlates with IC event multiplicity —
  late accidentals are a biased sample of IC events
- Angular features (Δθ LF/MPE/IterMPE) appeared at 4th–7th place — promising
  but buried under correlated scalars

## Key Physics Insight (Established Today)
The MPEFit vertex-time bias is **energy-dependent**: larger for high-multiplicity events
(Δt=604 ns at N=337 hits vs Δt=518 ns at N=50 hits). A flat 3σ window around 280 ns
clips genuine high-energy muon-DM-Ice hits. The first-hit anchor approach removes
most of this bias.

## Tomorrow's Tasks (Priority Order)

### 1. Check DM-Ice energy/amplitude in the i3 file [HIGH PRIORITY]
The single most powerful discriminant would be the NaI pulse amplitude:
- Genuine muon: deposits ~30–40 MeV in 15 cm NaI → large amplitude
- Radioactivity (U/Th/K-40): 1–3 MeV → small amplitude
- Look in the i3 file for keys like `DMIce_amplitude`, `DMIce_energy`,
  `DMIce_charge`, or any waveform/ADC associated with the DM-Ice trigger
- If available, this single cut could reduce accidentals by 10–100×

```bash
# On cobalt — inspect what DM-Ice keys are in the frame:
/cvmfs/.../env-shell.sh python3 -c "
from icecube import dataio, icetray
f = dataio.I3File('/data/user/bcharett/dmice_coincidences_2011_2022/all_dmice_coincidences_2011_2022_fixed.i3')
while f.more():
    frame = f.pop_frame()
    if frame.Stop == icetray.I3Frame.Physics:
        keys = [k for k in frame.keys() if 'DM' in k or 'dm' in k or 'NaI' in k]
        print(keys)
        break
"
```

### 2. Fix ML classifier labeling [HIGH PRIORITY]
Current label proxy (dm_t > 40 μs as background) correlates with IC multiplicity.

**Better approaches:**
- **A. Use BLO simulation as genuine training examples** — we have
  `sim_all_recos_det1_spline_mpe.csv` with known-genuine events. Extract
  angular differences from it and use as positive class. Use the late-accidental
  real events as negative class. Cross-domain training — will it generalize?
- **B. Angular-features-only classifier** — remove ALL scalar IC properties
  (n_hits, lf_zen, energy, n_doms) from features. Use ONLY the pairwise angular
  differences. This is the purest test of "do all reconstructions agree?"
- **C. Anomaly score within the signal window** — fit a density model on the
  angular differences for dm_t > 40 μs events (background), then score events
  in the 7–16 μs window by how anomalous their angular differences are relative
  to the background model. No circular labeling.

Run from: `~/dmice/ml_coinc_classifier.py` (already written, just needs changes)

### 3. IterMPE cleaning → better PivotLineFit [COBALT]
After IterMPE converges, compute per-DOM Pandel time residuals and keep only
DOMs within ±200 ns. Rerun PivotLineFit on clean DOMs. This removes outlier
hits that currently corrupt LineFit and PivotLineFit equally (masking the
angular difference between them).

Add to `compare_coinc_cuts_v2.py`:
```python
# After getting iter_mpe track:
for omk, pulses in pmap.items():
    t_hit = min(p.time for p in pulses)
    d_hat_iter = [iter_mpe.dir.x, ...]
    s = dot(dom_pos - iter_mpe.pos, d_hat_iter)
    t_pandel = iter_mpe.time + s/C + pandel_peak(d_perp)
    if abs(t_hit - t_pandel) < 200:   # keep clean DOMs
        clean_hits.append(...)
# Then run PivotLineFit on clean_hits
```

### 4. Energy-corrected timing cut [MEDIUM]
The first-hit anchor Δt is much better than MPEFit Δt, but can still be improved.
Parameterize the residual bias as a function of n_hits:
- For the 26 first-hit passing events, compute Δt_1st vs n_hits
- Fit a linear correction: Δt_corrected = Δt_1st - f(n_hits)
- Apply corrected Gaussian cut

### 5. Run all-years compare_coinc_cuts_v2 [COBALT]
Currently only ran on 2012. Run without --year flag to get all years 2012–2021.
This takes ~10 min — use a screen session.

## Files Written Today

| File | Purpose |
|------|---------|
| `~/dmice/plot_event_display_coinc.py` | Event display: 10 pass + 10 fail timing cut |
| `~/dmice/compare_coinc_cuts_v2.py` | Three discriminants: MPEFit Gaussian, first-hit, Δθ |
| `~/dmice/ml_coinc_classifier.py` | ML classifier (IsoForest + GradientBoosting) |
| `~/dmice_work/output/event_displays/` | 16 event display PNGs + montage |
| `~/dmice_work/output/coinc_cuts_v2.csv` | 2012 results for all three discriminants |
| `~/dmice_work/output/ml_coinc_score.csv` | ML scores for all 6000 events |
| `~/dmice_work/output/ml_coinc_plots.png` | ML diagnostic plots |
| `~/dmice_work/output/lf_vs_pivot_ang_diff.png` | LF vs PivotLF Δθ distribution analysis |

## Key Numbers to Remember

- 2012 events: 1222 total, 731 with LineFit
- MPEFit available: ~35.6% of 2012 events
- First-hit anchor: 26/731 pass 3σ (3.6%) — best current timing discriminant
- Expected genuine rate: ~50/year from physics (~4%)
- Accidental rate: ~96% of all events
- MPEFit vertex-time bias: ~215 ns average, up to ~324 ns for high-energy events
- True signal window in dm_t_ns: ~8000–15000 ns (muon transit from top of IC to crystal)

## Longer-Term (Not for Tomorrow)

- Full paper: systematic uncertainties, DM flux limit, IceCube collab review
- ICRC 2026 proceedings deadline — check actual date
- Run IterMPE cleaning + pivot on full 6000-event dataset (needs cobalt Condor)
- Compare Gaussian timing cut yield vs geometric cut for all years
