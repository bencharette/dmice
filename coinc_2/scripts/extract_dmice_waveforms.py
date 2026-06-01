#!/usr/bin/env python3
"""
extract_dmice_waveforms.py

Run on Cobalt. Reads DM-Ice ROOT files for all years and extracts per-event
waveform variables needed to reproduce Hubbard Fig 7.11 (tau vs pulse height)
and Fig 7.15 (energy spectrum).

Output: ~/dmice_work/output/dmice_waveforms.csv

Columns:
  year, month, run_id, trigger_time, detector, pmt, hv_era,
  pulse_height   -- max_atwd1 (lowHV) or max_atwd2 (highHV) [ADC counts]
  tau            -- atwd1_jitter_tau (lowHV) or atwd2_jitter_tau (highHV) [ns]
  energy         -- sum_128_atwd1 (lowHV) or sum_128_atwd2 (highHV) [ADC·bin]
  is_muon        -- True if event passes Hubbard muon selection cuts

Muon cuts (Hubbard thesis §7.2.1, p.121):
  lowHV PMT-1a/1b (DM0/DM1):  height > 650  OR (height > 325 AND tau > 177)
  lowHV PMT-2a/2b (DM2/DM3):  height > 400  OR (height > 324 AND tau > 184)
  highHV: same logical structure, thresholds TBD.

ROOT files: /data/exp/DM-Ice/{year}/filtered/pole/data/tree/{month}/2021_processing/*.root
"""

import os
import glob
import numpy as np
import pandas as pd
import uproot

ROOT_BASE = "/data/exp/DM-Ice"
OUT_CSV   = os.path.expanduser("~/dmice_work/output/dmice_waveforms.csv")
os.makedirs(os.path.dirname(OUT_CSV), exist_ok=True)

YEARS  = range(2012, 2022)
MONTHS = [f"{m:02d}" for m in range(1, 13)]

# Muon cut thresholds (Hubbard §7.2.1)
# muon if height > height_high  OR  (height > height_psd AND tau > tau_min)
CUTS = {
    "DM0": (650, 325, 177),   # PMT-1a, det1
    "DM1": (650, 325, 177),   # PMT-1b, det1
    "DM2": (400, 324, 184),   # PMT-2a, det2
    "DM3": (400, 324, 184),   # PMT-2b, det2
}

def atwd_channel(filename):
    bn = os.path.basename(filename)
    return ("atwd1", "lowHV") if "lowHV" in bn else ("atwd2", "highHV")

def is_muon(height, tau, dm_name):
    h_hi, h_psd, tau_min = CUTS[dm_name]
    return bool((height > h_hi) or (height > h_psd and tau > tau_min))

records = []
n_files = 0
n_events = 0

for year in YEARS:
    for month in MONTHS:
        pattern = f"{ROOT_BASE}/{year}/filtered/pole/data/tree/{month}/2021_processing/*.root"
        files = glob.glob(pattern)
        if not files:
            continue

        for fpath in files:
            ch, hv_era = atwd_channel(fpath)
            run_id = os.path.basename(fpath).split("_")[1]

            try:
                rf = uproot.open(fpath)
            except Exception as e:
                print(f"  WARN: cannot open {fpath}: {e}")
                continue

            for tree_name, dm_pairs, det in [
                ("Tree0", [("DM0", "pmt1a"), ("DM1", "pmt1b")], "det1"),
                ("Tree1", [("DM2", "pmt2a"), ("DM3", "pmt2b")], "det2"),
            ]:
                if tree_name not in rf:
                    continue
                tree = rf[tree_name]

                for dm, pmt_label in dm_pairs:
                    try:
                        times  = tree[f"{dm}_trigger_time"].array(library="np")
                        height = tree[f"{dm}_max_{ch}"].array(library="np")
                        tau    = tree[f"{dm}_{ch}_jitter_tau"].array(library="np")
                        energy = tree[f"{dm}_sum_128_{ch}"].array(library="np")
                    except Exception as e:
                        print(f"  WARN: missing branch {dm} in {fpath}: {e}")
                        continue

                    for i in range(len(times)):
                        records.append({
                            "year":         year,
                            "month":        int(month),
                            "run_id":       run_id,
                            "trigger_time": int(times[i]),
                            "detector":     det,
                            "pmt":          pmt_label,
                            "hv_era":       hv_era,
                            "pulse_height": float(height[i]),
                            "tau":          float(tau[i]),
                            "energy":       float(energy[i]),
                            "is_muon":      is_muon(float(height[i]), float(tau[i]), dm),
                        })
                        n_events += 1

            n_files += 1
            if n_files % 500 == 0:
                print(f"  {n_files} files, {n_events} events so far...", flush=True)

    print(f"Year {year} done. {n_files} files total.", flush=True)

print(f"\nTotal: {n_files} files, {n_events} events")
df = pd.DataFrame(records)
df.to_csv(OUT_CSV, index=False)
print(f"Saved: {OUT_CSV}")
n_muons = df["is_muon"].sum()
print(f"Muon events: {n_muons} / {len(df)} ({100*n_muons/len(df):.2f}%)")
