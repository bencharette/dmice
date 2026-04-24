#!/usr/bin/env python3
"""
Extract DM-Ice NaI amplitude data for each coincidence event and merge
into the existing real_all_recos.csv.

For each event in the merged i3 file we have DMIce_detection_time (float64).
The ROOT files at /data/exp/DM-Ice/{year}/filtered/pole/data/tree/{month}/2021_processing/
store DM0/DM1/DM2/DM3 trigger_time (int64). Due to float64 precision limits
at 1.5e17, detection_time differs from trigger_time by at most 32 DAQ ticks.
We match with tolerance ±64.

DM channel → detector mapping:
  DM0, DM1  →  det1 (string 87)
  DM2, DM3  →  det2 (string 88)

Amplitude channel (from vetoRootMaster.py logic):
  *lowHV.root      →  atwd1
  *highHV*.root    →  atwd2

Output columns added:
  dm_amp        max(PMT0, PMT1) for the relevant atwd channel (calibrated)
  dm_raw_amp    same but raw (pre-droop-correction)
  dm_sum_128    max(PMT0, PMT1) sum_128_atwd{ch} (integrated charge ≈ energy)
  dm_thresh_e   max(PMT0, PMT1) thresh_128_atwd{ch}
  dm_hv_era     'lowHV' or 'highHV'
  dm_match      True/False whether a ROOT match was found
"""

import os
import sys
import glob
import numpy as np
import pandas as pd
from icecube import icetray, dataio, dataclasses

MERGED_I3 = "/data/user/bcharett/dmice_coincidences_2011_2022/all_dmice_coincidences_2011_2022_fixed.i3"
INPUT_CSV  = os.path.expanduser("~/dmice_work/output/real_all_recos.csv")
OUTPUT_CSV = os.path.expanduser("~/dmice_work/output/real_all_recos_with_dmamp.csv")
ROOT_BASE  = "/data/exp/DM-Ice"
TOLERANCE  = 64  # DAQ ticks; float64 precision gives at most 32 error

# ---------------------------------------------------------------------------
# Step 1: Read merged i3 file → dict (run_id, event_id) → detection_time_int
# ---------------------------------------------------------------------------
print("Reading merged i3 file...")
i3_data = {}  # (run_id, event_id) -> detection_time_int
det_lookup = {}  # (run_id, event_id) -> detector string ('det1' or 'det2')

f = dataio.I3File(MERGED_I3)
n_frames = 0
while f.more():
    frame = f.pop_frame()
    if frame.Stop != icetray.I3Frame.Physics:
        continue
    hdr = frame["I3EventHeader"]
    run_id   = hdr.run_id
    event_id = hdr.event_id
    det_time = int(frame["DMIce_detection_time"].value)
    detector = frame["DMIce_detector"].value
    key = (run_id, event_id)
    i3_data[key]   = det_time
    det_lookup[key] = detector
    n_frames += 1

print(f"  {n_frames} physics frames loaded.")

# Build reverse lookup: detection_time → (run_id, event_id)
# Bucket by detection_time // 10000 (10-µs bins) for fast lookup
from collections import defaultdict
time_bucket = defaultdict(list)
for key, t in i3_data.items():
    time_bucket[t // 10000].append((t, key))

# ---------------------------------------------------------------------------
# Step 2: Scan ROOT files and match trigger_times
# ---------------------------------------------------------------------------

import uproot

def atwd_channel(filename):
    """Return 'atwd1' for lowHV, 'atwd2' for highHV."""
    bn = os.path.basename(filename)
    if "lowHV" in bn:
        return "atwd1", "lowHV"
    else:
        return "atwd2", "highHV"

def find_match(trigger_time, bucket_dict, tol=TOLERANCE):
    """Find (run_id, event_id) whose detection_time matches trigger_time ±tol."""
    bin_key = trigger_time // 10000
    for bk in (bin_key - 1, bin_key, bin_key + 1):
        for (t, key) in bucket_dict.get(bk, []):
            if abs(t - trigger_time) <= tol:
                return key
    return None

print("Scanning ROOT files and matching amplitudes...")

amp_records = {}  # (run_id, event_id) -> dict of amplitude values

years = range(2012, 2022)
matched = 0
not_found_in_root = 0

for year in years:
    pattern = f"{ROOT_BASE}/{year}/filtered/pole/data/tree/*/2021_processing/*.root"
    root_files = glob.glob(pattern)
    print(f"  {year}: {len(root_files)} ROOT files", flush=True)
    for rfile in root_files:
        ch, hv_era = atwd_channel(rfile)
        try:
            rf = uproot.open(rfile)
        except Exception as e:
            print(f"    WARNING: could not open {rfile}: {e}")
            continue

        # det1: Tree0 has DM0 and DM1
        # det2: Tree1 has DM2 and DM3
        for tree_name, dm_a, dm_b, det_str in [
            ("Tree0", "DM0", "DM1", "det1"),
            ("Tree1", "DM2", "DM3", "det2"),
        ]:
            if tree_name not in rf:
                continue
            tree = rf[tree_name]

            try:
                times_a = tree[f"{dm_a}_trigger_time"].array(library="np")
                amp_a   = tree[f"{dm_a}_max_{ch}"].array(library="np")
                raw_a   = tree[f"{dm_a}_raw_max_{ch}"].array(library="np")
                sum_a   = tree[f"{dm_a}_sum_128_{ch}"].array(library="np")
                thr_a   = tree[f"{dm_a}_thresh_128_{ch}"].array(library="np")

                times_b = tree[f"{dm_b}_trigger_time"].array(library="np")
                amp_b   = tree[f"{dm_b}_max_{ch}"].array(library="np")
                raw_b   = tree[f"{dm_b}_raw_max_{ch}"].array(library="np")
                sum_b   = tree[f"{dm_b}_sum_128_{ch}"].array(library="np")
                thr_b   = tree[f"{dm_b}_thresh_128_{ch}"].array(library="np")
            except Exception as e:
                print(f"    WARNING: missing branch in {rfile}/{tree_name}: {e}")
                continue

            # Events in Tree0/Tree1 are paired (same index = same physical event).
            # DM0 and DM1 are the two PMTs of the same crystal; take max per event.
            for i in range(len(times_a)):
                t_a = int(times_a[i])
                key = find_match(t_a, time_bucket)
                if key is None:
                    continue
                # Confirm the detector matches
                if det_lookup.get(key) != det_str:
                    continue
                # Take max of two PMTs
                amp_records[key] = {
                    "dm_amp":      max(float(amp_a[i]), float(amp_b[i])),
                    "dm_raw_amp":  max(float(raw_a[i]), float(raw_b[i])),
                    "dm_sum_128":  max(float(sum_a[i]), float(sum_b[i])),
                    "dm_thresh_e": max(float(thr_a[i]), float(thr_b[i])),
                    "dm_hv_era":   hv_era,
                    "dm_match":    True,
                }
                matched += 1

print(f"Matched {matched} events with amplitude data.")

# ---------------------------------------------------------------------------
# Step 3: Merge with CSV
# ---------------------------------------------------------------------------
print("Merging with CSV...")
df = pd.read_csv(INPUT_CSV)

amp_cols = ["dm_amp", "dm_raw_amp", "dm_sum_128", "dm_thresh_e", "dm_hv_era", "dm_match"]
for col in amp_cols:
    df[col] = np.nan if col != "dm_hv_era" and col != "dm_match" else None

df["dm_match"] = False
df["dm_hv_era"] = ""

for i, row in df.iterrows():
    key = (int(row["run_id"]), int(row["event_id"]))
    if key in amp_records:
        rec = amp_records[key]
        for col in amp_cols:
            df.at[i, col] = rec[col]

n_matched = df["dm_match"].sum()
print(f"CSV rows with amplitude: {n_matched} / {len(df)}")

df.to_csv(OUTPUT_CSV, index=False)
print(f"Written: {OUTPUT_CSV}")
