# DM-Ice Timing in Pivot Reconstructions — Status & Next Steps

## What We Built

We added DM-Ice NaI hit timing as a constraint to two IceCube muon track reconstructions,
then benchmarked all methods on 200-event BLO simulation (tracks aimed through DM-Ice det1/det2).

### Timing Model

Built from BLO simulation hits (`muons_binned_200ev_repacked.npz`), stored at:
```
~/dmice_work/output/dmice_timing_model.npz          # BLO-derived (μ = -41.9 ns, σ = 47 ns)
~/dmice_work/output/dmice_timing_model_calibrated.npz  # real-data calibrated (μ = +280 ns, σ = 81 ns)
```

The model predicts:
```
t_DM = t_geo + Δt,   Δt ~ N(μ, σ)
t_geo = t_PCA + d_perp / (c · sin θ_C)
```

The +322 ns offset between BLO (μ = -42 ns) and real 2012 data (μ = +280 ns) is NaI
scintillation delay — the crystal fires on scintillation light, not the prompt Cherenkov photon.

---

## The Two New Reconstructions

### 1. Pivot LineFit (`PivotLineFit`)

**How it works** (`pivot_linefit_ic()` in `run_splinempe_pivot.py`):

1. Run standard `I3LineFit` on IC-only pulses (strings 1–86).
2. Project the DM-Ice detector position onto the LineFit track direction to get an
   estimated transit time `t_dm` at the DM-Ice location.
3. Re-anchor the weighted least-squares LineFit using `t_dm` as the time origin,
   so the reconstructed track is forced to pass through DM-Ice at the observed hit time.
4. Normalize the resulting velocity vector to get the pivot direction.

**Result:** median angular error **1.32°** vs 5.05° for standard LineFit.

### 2. Pivot MPEFit (`MPEFit_Pivot`)

**How it works:**

1. Compute Pivot LineFit (above).
2. Seed `I3SinglePandelFitter` (MPE likelihood) with the Pivot LineFit direction instead
   of the standard LineFit direction.
3. MPEFit then optimizes over IC pulses using the DM-Ice-constrained starting point.

**Result:** median angular error **0.36°** vs 0.50° for standard MPEFit seed, and vs 9.02°
for SplineMPE (which performs poorly on these deep bottom-of-detector tracks).

### Benchmark Summary

| Method              | Median Angular Error | N events |
|---------------------|---------------------|----------|
| LineFit             | 5.05°               | 179      |
| Pivot LineFit       | 1.32°               | 151      |
| MPEFit (std seed)   | 0.50°               | 172      |
| **Pivot MPEFit**    | **0.36°**           | 151      |
| SplineMPE (std)     | 9.02°               | 169      |
| SplineMPE (pivot)   | 4.57°               | 147      |

SplineMPE is excluded going forward — it performs much worse than MPEFit on these
deep tracks and the `default` configuration is too coarse to be useful here.

---

## Known Problems — What To Fix Tomorrow

### Problem 1: μ correction is missing (most important fix)

The pivot time calculation currently uses the raw DM-Ice hit time with no correction
for the timing model offset μ. The hit time includes ~280 ns of NaI scintillation delay
(real data) or ~42 ns of Cherenkov/geometric offset (BLO sim) that biases the pivot anchor.

**Fix:** subtract μ before anchoring the pivot:
```python
# In pivot_linefit_ic(), before computing t_dm:
dm_t_corrected = dm_t - MU_NS   # MU_NS = -41.9 (BLO) or +280.0 (real data)
# then use dm_t_corrected instead of dm_t
```

This is a one-line change in `run_splinempe_pivot.py` and should improve both
Pivot LF and Pivot MPEFit. Use the calibrated model for real data.

### Problem 2: MPEFit seed has wrong t₀

The Pivot LineFit particle written to the frame uses `pp.time = lf.time` (the original
LineFit t₀). MPEFit then re-optimizes t₀ freely, discarding the DM-Ice time constraint.

**Fix:** set the pivot particle's position and time to the DM-Ice-constrained values:
```python
# In PivotLFModule.Physics():
# compute t0_pivot = dm_t_corrected - (projection along track) / C_M_NS
pp.pos  = dataclasses.I3Position(*dm_pos)   # vertex at DM-Ice detector
pp.time = t0_pivot                           # DM-Ice-constrained t0
```

This gives MPEFit a fully DM-Ice-anchored starting point in both direction and time.

### Problem 3: Timing uncertainty σ is ignored

The pivot uses a single point estimate of dm_t. The model uncertainty σ (47 ns BLO,
81 ns real) is not propagated into the seed.

**Fix (optional, 3× compute cost):** multi-seed over timing uncertainty:
```python
for t_offset in [0.0, +SIGMA_NS, -SIGMA_NS]:
    # build pivot seed with dm_t_corrected + t_offset
    # run MPEFit from each seed
    # keep best-likelihood result
```

### Problem 4: Real data uses wrong μ (BLO model, not calibrated)

The scoring in `dmice_likelihood.py` and the pivot in `run_splinempe_pivot.py` both
load `dmice_timing_model.npz` (BLO, μ = -41.9 ns). For real coincidence data the
calibrated model should be used.

**Fix:** switch model path when processing real data:
```python
MODEL_FILE = os.path.expanduser("~/dmice_work/output/dmice_timing_model_calibrated.npz")
```

---

## Suggested Order of Work Tomorrow

1. **Apply μ correction** to `pivot_linefit_ic()` — one line, rerun `run_splinempe_pivot.py`
   on cobalt and check that Pivot LF and Pivot MPEFit medians improve.

2. **Fix MPEFit seed t₀** — set `pp.pos` and `pp.time` from the DM-Ice-constrained values
   rather than copying from LineFit. Rerun and compare.

3. **Re-benchmark** — the updated plot should show Pivot MPEFit improving below 0.36°.

4. **Apply to real data** — run the calibrated model on the real coincidence I3 file,
   deduplicating the 2013/2018 pipeline artifact first (each `dm_raw` value should appear
   at most once in the merged coincidence list).

5. **Phase D (longer term)** — implement DM-Ice timing as a gulliver `I3EventLogLikelihood`
   term so MPEFit minimizes over IC pulses + DM-Ice timing simultaneously. This is the
   ceiling of what DM-Ice can contribute to the reconstruction.

---

## Uncertainty Propagation — Questions for Math Physics Professor

Three questions to ask, ordered by relevance to the analysis:

### Q1: Nonlinear propagation through the geometric time formula (most important)

> "If the observed hit time is the sum of a geometric term plus a residual that's a
> convolution of an ice-scattering PDF (Pandel) and an exponential scintillation decay,
> how do I propagate uncertainty in the track parameters — direction angles θ, φ and
> vertex position — through the geometric time formula into the total timing residual?"

**Why this matters:** The pivot anchor uses
```
t_geo = t_PCA + d_perp / (c · sin θ_C)
```
where `d_perp` and `t_PCA` both depend nonlinearly on track direction and vertex.
First-order (Gaussian) error propagation — δt = |∂t/∂x| δx — is what we're implicitly
assuming when we subtract μ and use σ as a flat uncertainty. The professor can show
whether that approximation holds and what the correction terms look like.

### Q2: Why Gaussian propagation fails for non-Gaussian residuals

> "When the underlying PDF is non-Gaussian — say, a Pandel distribution with a long
> positive tail — why does standard first-order error propagation break down, and what
> is the correct approach for estimating parameter uncertainty in a maximum-likelihood fit?"

**Why this matters:** Our timing model fits a Gaussian to Δt residuals, but the true
distribution has a Pandel (ice scattering) tail on the positive side. This means σ
underestimates large-residual events and the multi-seed approach (Problem 3 above)
should really be sampling from the full PDF, not ±σ of a Gaussian.

### Q3: First-photon statistics and how τ_eff scales with photon count

> "For a scintillator with exponential decay constant τ ≈ 250 ns, if we trigger on the
> *first* photon out of n photoelectrons, what is the distribution of the first-arrival
> time, and how does the mean and variance scale with n?"

**Why this matters:** NaI fires on the first scintillation photon out of however many
are produced by the muon. By order statistics, the minimum of n i.i.d. Exp(1/τ) random
variables is Exp(n/τ), so the effective mean shifts to τ/n and variance to (τ/n)².
Our fixed-μ correction ignores this: high-energy muons produce more photons, so the
effective scintillation delay is shorter. This is part of why σ differs between BLO sim
(lower energy) and real data.

---

## Key Files

| File | Purpose |
|------|---------|
| `~/dmice/run_splinempe_pivot.py` | Main tray script: NPZInjector → Pivot LF → MPEFit × 2 → Scorer |
| `~/dmice/dmice_likelihood.py` | Timing model class + score_sim_npz() |
| `~/dmice/score_real_coincidences.py` | Real data sanity check (cobalt only) |
| `~/dmice_work/output/dmice_timing_model.npz` | BLO-derived model (use for sim) |
| `~/dmice_work/output/dmice_timing_model_calibrated.npz` | Real-data model (use for data) |
| `~/dmice_work/output/muons_binned_200ev_repacked.npz` | 200-event BLO sim input |
| `~/dmice_work/output/mpe_pivot_comparison.png` | Current benchmark plot |
