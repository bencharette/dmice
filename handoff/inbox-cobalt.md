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
