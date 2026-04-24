# Session: gauss-coincidence
Created: 2026-04-20

## Goal

Investigate using the NaI Gaussian timing model as a coincidence selection cut,
as an alternative to the current geometric d_perp < 15 m cut.

## Background

The existing coincidence finder (`find_dmice_coincidences.py`) uses a raw ±70 μs time
window to select IceCube events coincident with DM-Ice triggers. The NaI timing model
(Gaussian, μ=+280 ns, σ=81 ns) is not used in selection at all — only in reconstruction.

**Three cuts to compare:**
- **Geometric:** d⊥(reco → DM-Ice crystal) < 15 m
- **Gaussian:** |Δt − 280 ns| < 3σ = ±243 ns, where Δt = t_DM − t_geo(reco)
- **Both:** intersection

## Pilot Results (2012, 1222 events)

| | LineFit | MPEFit |
|--|---------|--------|
| Geometric | 156 | 262 |
| Gaussian | 93 | 50 |
| Both | 20 | 31 |

MPEFit Gaussian (50 events) is the highest-purity sample.
Only 35.6% of 2012 events have MPEFit available in the master i3 file.

## Key Script

**`~/dmice/compare_coinc_cuts.py`** — runs on cobalt with IceTray.

```bash
# On cobalt:
/cvmfs/icecube.opensciencegrid.org/py3-v4.3.0/RHEL_9_x86_64/metaprojects/icetray/v1.12.1/env-shell.sh \
  python3 -u ~/dmice/compare_coinc_cuts.py [--year YYYY]
```

Outputs: `~/dmice_work/output/coinc_cut_comparison.{csv,png}`

## VtxMPE2D Results — All Years (added 2026-04-19)

A complementary approach: **2D vertex-constrained MPE fit** that jointly optimizes IC direction
and DM-Ice timing. Script: `~/dmice/run_2d_vtx_mpe.py`.

**Run command on cobalt:**
```bash
/cvmfs/icecube.opensciencegrid.org/py3-v4.3.0/RHEL_9_x86_64/metaprojects/icetray/v1.12.1/env-shell.sh \
  python3 -u ~/dmice/run_2d_vtx_mpe.py \
    --out ~/dmice_work/output/vtx2d_mpe_all.csv
# CRITICAL: GCD file is needed (passed via --gcd, defaults to 2013 GCD)
```

**All-years results (2012-2021, 6518 events, 5042 valid fits):**

| Year | Events | Fits | 3σ candidates | 1σ candidates |
|------|--------|------|---------------|---------------|
| 2012 | 1222   | 670  | **6**         | 1             |
| 2013 | 1129   | 737  | 0             | 0             |
| 2014 | 851    | 688  | 0             | 0             |
| 2015 | 646    | 540  | 1             | 1             |
| 2016 | 418    | 358  | 1             | 1             |
| 2017 | 510    | 457  | 0             | 0             |
| 2018 | 493    | 445  | 0             | 0             |
| 2019 | 488    | 447  | 0             | 0             |
| 2020 | 398    | 370  | 1             | 1             |
| 2021 | 363    | 330  | 1             | 0             |

**Total:** 10/5042 (0.2%) pass 3σ, 4/5042 (0.1%) pass 1σ cut.

**Key diagnostic:** `vtx2d_dm_dt` = dm_t_corrected − t₀(fitted direction)
- Background peaks at ~−5000 ns (accidental DM-Ice fires early relative to IC track)
- Signal region: |dm_dt| < 243 ns
- 2012 excess (6 events) is anomalous — worth investigating

**Plot:** `~/dmice_work/output/vtx2d_dm_dt_dist.png`
**CSV:** `~/dmice_work/output/vtx2d_mpe_all.csv`

**Best signal candidate:** 2012 event run=121207 evt=17314026
- dm_dt = 6.2 ns (essentially perfect timing match)
- vtx2d_zen = 134.7° (upgoing — atmospheric neutrino?)
- 6 DOMs (low multiplicity)
- mpe_zen = 127.7° (consistent with standard reco)

**CRITICAL — 2012 anomaly explained:**
2012 has a dramatically wider dm_t_ns distribution: p75=13,853 ns vs ~7,700 ns for 2013-2021.
ALL 10 signal candidates have dm_t_ns > 6,000 ns. For a genuine downgoing muon-DM-Ice hit,
expected dm_t_ns ≈ 1,000-5,000 ns (transit time: 500-1500 m / 0.3 m/ns). The candidates
at 12,000-15,000 ns are almost certainly accidentals where DM-Ice radioactivity fires and
happens to match the IC track timing by chance.

**Implication:** Need an additional cut `dm_t_ns < 5000 ns` to isolate genuine muon hits.
The entire coincidence dataset (~1220 events/year) is dominated by accidentals. True
muon-DM-Ice hits are rare (estimated ~0.3/year from crystal cross-section + muon flux).

## Immediate Next Steps

1. Run `compare_coinc_cuts.py` on all years (no `--year` flag) — get full picture
2. Apply additional cut: `dm_t_ns < 5000 ns` to select genuine muon-DM-Ice hits
   - Rerun VtxMPE2D signal search with this pre-cut
   - Check how many events survive and what dm_dt they have
3. Understand 2012 anomaly: check if coincidence window was different in 2012
   - Look at the coincidence finder script for year-specific window settings
4. Test angular resolution: do the 50 MPE-Gaussian events reconstruct better than
   the 156-event geometric sample? Run pivot reco on both subsets and compare.
5. Investigate why 64.4% of 2012 events lack MPEFit — era/stream issue?
6. Try N_sigma = 1, 2, 3 to map yield vs purity

## Master Data File

`/data/user/bcharett/dmice_coincidences_2011_2022/all_dmice_coincidences_2011_2022_fixed.i3`
(on NPX/cobalt — use the `.i3` symlink, NOT `.i3.zst` despite the extension)

## Context

- This is a **sidebar** to the main DM-Ice reconstruction work — not for the upcoming meeting
- The meeting focus is using the Gaussian timing IN the reconstruction likelihood
- This session is about using it in the coincidence SELECTION (upstream step)
- See memory: `project_gauss_coincidence.md` and `project_dmice.md` for full context
