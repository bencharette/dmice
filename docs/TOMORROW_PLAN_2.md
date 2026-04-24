# Tomorrow's Plan 2 — DM-Ice NaI Amplitude Extraction

## Context

We discovered that DM-Ice NaI crystal amplitude data is accessible on Cobalt at:
```
/data/exp/DM-Ice/{year}/filtered/pole/data/tree/{month}/2021_processing/dmice_run{N}_LC23-*.root
```
These are pre-processed ROOT files (uproot-readable, no py2 env needed) created
by R. Clark's pipeline (`vetoRootMaster.py`).

---

## Background

Each ROOT file has Tree0 (det1: DM0+DM1 PMTs) and Tree1 (det2: DM2+DM3 PMTs).
Key amplitude branches (per PMT, two PMTs per crystal, take max):
- `DM{n}_max_atwd1`      — peak amplitude, calibrated (lowHV era: 2012–2014)
- `DM{n}_raw_max_atwd1`  — peak amplitude, raw
- `DM{n}_max_atwd2`      — peak amplitude, calibrated (highHV era: 2015+)
- `DM{n}_sum_128_atwd1`  — integrated charge over 128 bins (~energy)
- `DM{n}_thresh_128_atwd1` — threshold-based energy estimate
- `DM{n}_trigger_time`   — DAQ time (int64); matches `DMIce_detection_time`
                           in i3 file to within ±32 DAQ ticks (float64 precision)

HV era by year:
- 2012, 2013, 2014 → `*lowHV.root` → use `atwd1`
- 2015             → `*highHV2.root` → use `atwd2`
- 2016–2021        → `*highHV2-v9.root` → use `atwd2`

Matching verified: detection_time (i3 float64) vs trigger_time (ROOT int64)
differ by exactly 16 DAQ ticks for a test event. Tolerance ±64 is safe.

---

## Step 1: Check the overnight job result

The extraction script was launched on Cobalt in `screen dmamp`:
```bash
ssh cobalt-14
screen -r dmamp          # check progress
# or if done:
cat ~/dmice_work/output/dmamp.log
```

**Expected output:** `real_all_recos_with_dmamp.csv` at `~/dmice_work/output/`
with 6000 rows, new columns: `dm_amp`, `dm_raw_amp`, `dm_sum_128`, `dm_thresh_e`,
`dm_hv_era`, `dm_match`.

**Check:** `dm_match` should be True for most events. If < 80% matched, investigate:
- Might be events from years with different ROOT file naming (check years 2019–2021)
- Some step3 events might predate the 2021_processing reprocessing

Script: `~/dmice/extract_dm_amplitude.py`

---

## Step 2: Copy result locally

```bash
scp bcharett@cobalt-14.icecube.wisc.edu:~/dmice_work/output/real_all_recos_with_dmamp.csv \
    ~/dmice_work/output/
```

---

## Step 3: Explore amplitude distributions

Write `plot_dm_amplitude.py` to examine:

1. **Amplitude histogram** — `dm_amp` distribution. Muons through NaI should
   show a Landau-like peak. Background (noise, bipo) should peak at lower values.
   Compare lowHV vs highHV eras (units differ — ADC counts, but calibration differs).

2. **Amplitude vs IceCube energy** — scatter of `dm_amp` vs `energy_GeV` from the CSV.
   Should see correlation since both measure muon dE/dx.

3. **Amplitude vs d_perp** — check if pivot events (d⊥ < 15m) have systematically
   larger amplitude (closer track → more energy deposited in NaI crystal).

4. **Amplitude vs dm_t_ns** — genuine coincidences (dm_t_ns ≈ 8000–15000 ns)
   should have different amplitude distribution than accidentals (uniform dm_t_ns).

---

## Step 4: Use amplitude as discriminant

The DM-Ice NaI amplitude provides a fully independent discriminant:
- **Real muon through crystal:** large amplitude (Landau peak)
- **Accidental coincidence:** NaI triggered by noise/bipo/gamma → smaller amplitude

Proposed cut study:
```python
# Rough muon amplitude cut (calibrate from distribution)
muon_mask = df['dm_amp'] > MUON_THRESHOLD
```

Add `dm_amp > threshold` to the coincidence discriminant analysis
(alongside the existing Gaussian timing cut and geometric d_perp cut).
This is a third independent discriminant — potentially very powerful.

---

## Key Files

| File | Location |
|------|----------|
| Extraction script | `~/dmice/extract_dm_amplitude.py` |
| Extraction log | `~/dmice_work/output/dmamp.log` (on cobalt) |
| Amplitude CSV output | `~/dmice_work/output/real_all_recos_with_dmamp.csv` |
| ROOT data source | `/data/exp/DM-Ice/{year}/filtered/pole/data/tree/{month}/2021_processing/` |
| rclark pipeline ref | `/data/user/rclark/DMIce_sharable/vetoRootMaster.py` |

---

## Notes / Watch Out For

- Amplitude units are **not** in physical energy — they are calibrated ADC counts.
  The `thresh_128_atwd1` branch is R. Clark's energy estimator but units are unclear.
  Cross-check with published DM-Ice muon spectrum (arXiv:1509.02486).

- The ROOT files have DM0+DM1 paired by row index (same physical event = same row).
  We take `max(DM0_amp, DM1_amp)` as the crystal amplitude.

- A few step3 events may not have a corresponding ROOT file if the DM-Ice run
  was reprocessed differently. `dm_match=False` rows should be investigated.

- If the overnight job is still running, wait — it scans ~3600 ROOT files across
  10 years. Estimated runtime: 20–60 minutes.
