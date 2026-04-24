# Inbox for COBALT Machine

Messages from LOCAL and NPX for the Cobalt machine.

---

## 2026-04-20 — Event display plots for timing cut comparison

New script: `~/dmice/plot_event_display_coinc.py`

Run (after `git pull`):
```bash
/cvmfs/icecube.opensciencegrid.org/py3-v4.3.0/RHEL_9_x86_64/metaprojects/icetray/v1.12.1/env-shell.sh \
  python3 -u ~/dmice/plot_event_display_coinc.py --year 2012
```

**What it does:** Scans the master coincidence file for:
- 10 events that PASS the MPEFit Gaussian timing cut: |Δt_mpe − 280| < 243 ns
- 10 events that FAIL (accidentals, Δt_mpe ≫ 280 ns)

Plots x-z, y-z, x-y projections for each event:
- IC DOM hits (size ∝ charge, colour ∝ hit time)
- LineFit track (cyan line)
- DM-Ice crystal position (gold star at z ≈ −511 m)

Output: `~/dmice_work/output/event_displays/` (20 individual PNGs + 1 montage)
The montage `montage_pass_vs_fail.png` shows top 5 pass vs 5 fail side-by-side.

Copy back: `rsync -av cobalt-14:~/dmice_work/output/event_displays/ ~/dmice_work/output/event_displays/`

---

## 2026-04-21 — DM-Ice amplitude extraction job running

Script `~/dmice/extract_dm_amplitude.py` was launched in `screen -S dmamp`.
It scans all ROOT files at `/data/exp/DM-Ice/{year}/.../2021_processing/`
and matches each coincidence event (from merged i3) by trigger_time (±64 ticks).

**Monitor:** `screen -r dmamp` or `tail -f ~/dmice_work/output/dmamp.log`

**Expected output:** `~/dmice_work/output/real_all_recos_with_dmamp.csv`
New columns: dm_amp, dm_raw_amp, dm_sum_128, dm_thresh_e, dm_hv_era, dm_match

If the job fails partway: check the log for which year it reached.
Script can be re-run; output is overwritten.

See `docs/TOMORROW_PLAN_2.md` for next steps once the CSV is ready.

---

## 2026-04-23 — dmice/ directory reorganized

All loose files have been moved into subfolders. **Do a `git pull` before running anything new.**

Key path changes:

| Old path | New path |
|----------|----------|
| `~/dmice/run_splinempe_nai_lambda.py` | `~/dmice/reco/run_splinempe_nai_lambda.py` |
| `~/dmice/run_splinempe_pivot.py` | `~/dmice/reco/run_splinempe_pivot.py` |
| `~/dmice/run_all_recos_real.py` | `~/dmice/reco/run_all_recos_real.py` |
| `~/dmice/merge_2020_2021.py` | `~/dmice/analysis/merge_2020_2021.py` |
| `~/dmice/extract_dm_amplitude.py` | `~/dmice/analysis/extract_dm_amplitude.py` |
| `~/dmice/plot_event_display_coinc.py` | `~/dmice/plots/plot_event_display_coinc.py` |
| `~/dmice/overnight_run.sh` | `~/dmice/pipeline/overnight_run.sh` |
| `~/dmice/phase34_run.sh` | `~/dmice/pipeline/phase34_run.sh` |
| `~/dmice/steamshovel_artists.py` | `~/dmice/tools/steamshovel_artists.py` |
| `~/dmice/RESULTS.md` | `~/dmice/docs/RESULTS.md` |
| `~/dmice/TOMORROW_PLAN_2.md` | `~/dmice/docs/TOMORROW_PLAN_2.md` |

Full structure: `sim/` `reco/` `analysis/` `plots/` `pipeline/` `tools/` `docs/`

The **current lambda scan** (`reco/run_splinempe_nai_lambda.py`) is unaffected — it was already running before the move.

A newer local version exists as `reco/run_nai_lambda_pandel.py` (vertex-anchored Pandel approach, not yet run on Cobalt).

---
