# Tomorrow's Work Plan — DM-Ice Reconstruction

## Goal
Improve Pivot LineFit and Pivot MPEFit by fixing known timing bugs, add SPEFit for
low-energy events, and benchmark everything broken out by energy bin.

---

## Step 1: Fix the μ correction in pivot_linefit_ic()

**File:** `~/dmice/run_splinempe_pivot.py`

**What:** The pivot anchor currently uses raw DM-Ice hit time with no offset correction.
The NaI scintillation delay (μ) biases the pivot time systematically.

**Where:** In `PivotLFModule.Physics()`, before calling `pivot_linefit_ic()`:
```python
# Current (wrong):
dm_t = frame[DM_T_KEY].value

# Fixed:
dm_t = frame[DM_T_KEY].value - MU_NS   # MU_NS loaded from timing model
```

**Expected effect:** Pivot LineFit and Pivot MPEFit medians should both improve.

---

## Step 2: Fix the MPEFit seed t₀

**File:** `~/dmice/run_splinempe_pivot.py`

**What:** `PivotLFModule` currently copies `pp.time = lf.time` (original LineFit t₀).
MPEFit then re-optimizes t₀ freely, discarding the DM-Ice time constraint entirely.

**Where:** In `PivotLFModule.Physics()`, after computing the pivot direction, set the
vertex and time from the DM-Ice-constrained values:
```python
# Current (loses DM-Ice time info):
pp.pos  = lf.pos
pp.time = lf.time

# Fixed (anchor vertex at DM-Ice detector with constrained t0):
# t0_pivot = dm_t_corrected - (projection of dm_pos onto track) / C_M_NS
s = np.dot(dm_pos - np.array([lf.pos.x, lf.pos.y, lf.pos.z]), piv_dir)
t0_pivot = dm_t - s / C_M_NS   # dm_t already μ-corrected from Step 1
pp.pos  = dataclasses.I3Position(float(dm_pos[0]), float(dm_pos[1]), float(dm_pos[2]))
pp.time = float(t0_pivot)
```

**Expected effect:** MPEFit starts from a fully DM-Ice-anchored position + time,
not just a DM-Ice-anchored direction. Should improve Pivot MPEFit below 0.36°.

---

## Step 3: Add SPEFit (std + pivot seeds)

**File:** `~/dmice/run_splinempe_pivot.py`

**What:** MPEFit uses the multi-photon Pandel PDF, which is wrong at low energies
(< ~1 TeV) where most DOMs see only 1 photon. SPEFit uses the single-photon PDF
and should outperform MPEFit in the low-energy bins.

**Add two new frame keys:**
```python
SPE_STD_KEY = "SPEFit_Std"
SPE_PIV_KEY = "SPEFit_Pivot"
```

**Add to tray (after MPEFit blocks):**
```python
tray.Add(icecube.lilliput.segments.I3SinglePandelFitter,
    fitname = SPE_STD_KEY,
    domllh  = "SPE1st",
    pulses  = IC_PULSES,
    seeds   = [LF_KEY],
    If      = lambda f: LF_KEY in f and len(f[IC_PULSES]) >= 4,
)
tray.Add(icecube.lilliput.segments.I3SinglePandelFitter,
    fitname = SPE_PIV_KEY,
    domllh  = "SPE1st",
    pulses  = IC_PULSES,
    seeds   = [PIV_LF_KEY],
    If      = lambda f: PIV_LF_KEY in f and len(f[IC_PULSES]) >= 4,
)
```

**Add to Scorer rows dict:**
```python
spe_std_ang_err = ang(SPE_STD_KEY),
spe_piv_ang_err = ang(SPE_PIV_KEY),
```

---

## Step 4: Break benchmark out by energy bin

**What:** The 200-event sim has 5 energy bins (100 GeV–100 TeV, log-spaced).
The NPZ has `bin_id` per event (0–4). Report median angular error per bin
for each method to find the MPEFit/SPEFit crossover energy.

**Add to the summary block at the end of the script:**
```python
for bin_id in sorted(df.bin_id.unique()):
    sub = df[df.bin_id == bin_id]
    e_med = sub.mc_energy_GeV.median()
    print(f"\nBin {bin_id} (median E={e_med:.0f} GeV, n={len(sub)}):")
    for col, label in [
        ("lf_ang_err",      "LineFit"),
        ("piv_lf_ang_err",  "Pivot LF"),
        ("mpe_std_ang_err", "MPEFit (std)"),
        ("mpe_piv_ang_err", "MPEFit (piv)"),
        ("spe_std_ang_err", "SPEFit (std)"),
        ("spe_piv_ang_err", "SPEFit (piv)"),
    ]:
        vals = sub[col].dropna()
        if len(vals):
            print(f"  {label:20s}: {vals.median():.2f}° (n={len(vals)})")
```

**Also add a per-bin panel to the plot** — 5 subplots, one per energy bin, showing
MPEFit vs SPEFit angular error distributions.

---

## Step 5: Run on cobalt and collect results

```bash
scp ~/dmice/run_splinempe_pivot.py bcharett@cobalt-14:~/dmice/
ssh cobalt-14
screen -r splinempe_pivot   # or start new: screen -S run1
/cvmfs/icecube.opensciencegrid.org/py3-v4.3.0/RHEL_9_x86_64/metaprojects/icetray/v1.12.1/env-shell.sh \
  python3 ~/dmice/run_splinempe_pivot.py \
  > ~/dmice_work/output/splinempe_pivot_comparison.log 2>&1
```

Then `scp` the `.png` and `.csv` back locally.

---

## Step 6 (optional, if time): Multi-seed timing uncertainty

If Steps 1–2 improve results noticeably, try propagating σ into the pivot seed:
run MPEFit and SPEFit from three seeds — `dm_t ± σ` and `dm_t` — keep the
best-likelihood result. This handles events where the DM-Ice timing has a large
single-photon fluctuation.

---

## Expected Outcome After Steps 1–4

| Method           | Current | Expected after fixes |
|------------------|---------|----------------------|
| Pivot LineFit    | 1.32°   | ~1.0° (μ fix)        |
| MPEFit (std)     | 0.50°   | 0.50° (unchanged)    |
| Pivot MPEFit     | 0.36°   | < 0.30° (μ + t₀ fix) |
| SPEFit (std)     | —       | TBD, likely better at low E |
| Pivot SPEFit     | —       | TBD                  |

---

## Key Files

| File | Purpose |
|------|---------|
| `~/dmice/run_splinempe_pivot.py` | Main script to edit for Steps 1–4 |
| `~/dmice_work/output/muons_binned_200ev_repacked.npz` | 200-event BLO sim input (cobalt) |
| `~/dmice_work/output/dmice_timing_model.npz` | Timing model (μ=-41.9ns, σ=47ns) |
| `~/dmice_work/output/splinempe_pivot_comparison.{csv,png}` | Current benchmark output |
| `~/dmice_work/output/mpe_pivot_comparison.png` | Clean 4-method plot (local) |
| `~/dmice/DMICE_TIMING_IMPLEMENTATION.md` | Implementation reference |
