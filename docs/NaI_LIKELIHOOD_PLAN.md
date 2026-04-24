# DM-Ice NaI Likelihood Optimization — What We're Doing and Why

## The problem in one sentence

We can reconstruct muon tracks through IceCube using Cherenkov light (SplineMPE), but we're not yet using the fact that the muon **also physically hit the DM-Ice NaI crystal** at a known position and recorded time. This document explains how we're going to use that extra information to improve the reconstruction.

---

## Background: how SplineMPE works

SplineMPE is IceCube's standard muon track reconstruction. It fits a track hypothesis (position, direction, time) by maximizing a **likelihood** — a measure of how well the hypothesis explains the Cherenkov photon arrival times seen at each DOM.

Mathematically, it finds the track that maximizes:

```
log L_SplineMPE = Σ_DOMs  log p(t_hit | track)
```

where `p(t_hit | track)` is the probability of seeing a photon at that time, given the track, computed from spline-interpolated Monte Carlo tables.

---

## The new information: NaI timing

When a muon passes through the DM-Ice NaI scintillator, the crystal emits light (via ionization, not Cherenkov). The PMT records a hit time **T_meas**.

Key physics:
- The muon passes through the NaI crystal at some true transit time **T_transit**
- The scintillation adds a delay: NaI has a characteristic decay time around **~280 ns**
- There is also timing jitter of **~81 ns**
- So: `T_meas ≈ T_transit + Gaussian(280 ns, 81 ns)`

This means: if we know T_meas, and we know the NaI detector position exactly, we know **when the muon was at that specific (x, y, z) location** — with 81 ns uncertainty.

This is a powerful constraint because:
1. The DM-Ice position is precisely known
2. Every event in our dataset is a DM-Ice coincidence (the muon definitely crossed it)
3. Knowing the muon's position at one moment in time strongly constrains both direction and timing

---

## What we're adding to the likelihood

We add one extra term per event:

```
log L_total = log L_SplineMPE(track)  +  λ · log G(T_meas ; t_pred + 280ns, 81ns)
```

where:
- `t_pred(track)` = the time the track hypothesis predicts the muon reaches the DM-Ice crystal
- `log G(...)` = the Gaussian log-probability of the observed hit time given the prediction
- **λ** = a weight (hyperparameter) controlling how much we trust the NaI timing vs the Cherenkov timing

When λ = 0: standard SplineMPE, NaI ignored  
When λ → ∞: hard constraint — the track *must* be consistent with the NaI time  
When λ = λ_optimal: the best tradeoff between the two information sources

---

## Finding the optimal λ: hyperparameter optimization

We don't know the right λ in advance — it depends on the NaI timing resolution, the number of IceCube hits, and the muon energy. So we **scan over λ** and find the value that gives the best reconstruction across the full dataset.

### The loss function

For each λ, we run the combined reconstruction on all 5000 simulated events and compute:

| Loss | Formula | What it measures |
|------|---------|-----------------|
| `loss_ang_mean` | `mean(Δψ)` | Mean angular error vs MC truth |
| `loss_dp_mean` | `mean(d⊥)` | Mean distance from track to DM-Ice |
| `loss_combined` | `mean(Δψ) + mean(d⊥)/100` | Both together (normalised) |
| `loss_huber` | `mean(Huber(Δψ, δ=0.5°))` | Robust version — ignores events with already-tiny errors, strongly penalises large failures |

The **Huber loss** is particularly useful here. It is flat (zero penalty) for events where angular error is already small (below δ), and linear outside that region. This means it focuses the optimization on **fixing the badly reconstructed events** (the local-minimum failures we saw in bins 1 and 2) rather than trying to shave fractions of a degree off already-good events.

### The scan

λ values tested: `[0, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100]`

For each λ:
1. SplineMPE result (already computed) is used as seed
2. `scipy.optimize.minimize` (Nelder-Mead) refines the direction with the combined likelihood
3. Loss functions are evaluated across all events and per energy bin

### The upper bound run (`--true-time`)

We also run the scan using the **MC true transit time** (bypassing the 81 ns NaI timing uncertainty). This tells us: even with perfect NaI timing knowledge, how much can we improve? This is the theoretical upper bound.

---

## Why this helps with the failed events

Looking at the event displays from earlier, bins 1 and 2 (1 TeV and 2 TeV targeted at det2) had IterMPE reconstruction failures of 79° and 76°. These failed because:
- The hit cluster was compact near det2
- LineFit found a bad seed direction
- The reconstruction converged to a wrong local minimum

The NaI term fixes this by adding a penalty whenever the reconstructed track is inconsistent with the NaI hit time. A 79° wrong-direction track would predict a completely wrong transit time at DM-Ice — giving a huge penalty. This pulls the optimizer away from the wrong minimum.

---

## Implementation: what each script does

### `run_splinempe_pivot.py` (existing, modified)
Runs SplineMPE (and MPEFit, SPEFit) with standard and pivot seeds. Now also records SplineMPE angular errors in the output CSV (`smpe_std_ang_err`, `smpe_piv_ang_err`).

### `run_splinempe_nai_lambda.py` (new)
The λ scan script. Run on Cobalt with IceTray:
```bash
/cvmfs/.../env-shell.sh python3 ~/dmice/run_splinempe_nai_lambda.py
```
Output: `~/dmice_work/output/nai_lambda_scan.csv` + plot

For the upper bound:
```bash
/cvmfs/.../env-shell.sh python3 ~/dmice/run_splinempe_nai_lambda.py --true-time \
  --out ~/dmice_work/output/nai_lambda_scan_truetime.csv
```

### `plot_nai_lambda_scan.py` (to be created)
Generates ICRC2026-ready plots from the CSV.

---

## What we expect to see

If the NaI timing constraint helps:
- The loss function will have a clear minimum at some λ > 0
- The per-bin angular error plots will show improvement, especially for the low-energy bins targeting det2
- The `--true-time` run will show a deeper minimum (the 81 ns jitter is limiting)

If the NaI timing constraint doesn't help:
- The loss function will be flat or increasing for all λ > 0
- This would suggest the local-minimum problem is not fixed by adding a Gaussian penalty
- In that case: **constrained fit** (fix track to pass exactly through DM-Ice, reduce DOF from 6 to 2) or **ML approach** (MLP trained on SplineMPE output + NaI residual)

---

## Fallback plan: constrained fit

If λ optimization fails, the next step is a **hard constraint**: instead of a soft Gaussian penalty, we require the track to exactly pass through DM-Ice at the NaI-corrected time. This:
- Fixes the vertex at r_DM (the DM-Ice position)
- Fixes t0 = T_meas − 280 ns
- Leaves only (θ, φ) as free parameters — a 2D optimization instead of 6D

This is much harder for the optimizer to escape into wrong local minima because the solution space is much smaller.

---

## Fallback plan: ML correction

If the constrained fit also fails, we train a simple neural network (MLP) to correct the SplineMPE output. Key input features:
- SplineMPE direction (zenith, azimuth)
- **NaI time residual**: `T_meas − t_pred(SplineMPE) − 280 ns` ← the key new feature
- Number of hit DOMs, total charge, charge-weighted centroid position
- Time spread of IC hits
- Energy bin

Output: corrected angular direction. The NaI residual is the feature that tells the network "how wrong is SplineMPE relative to the DM-Ice constraint", which it can then learn to correct.
