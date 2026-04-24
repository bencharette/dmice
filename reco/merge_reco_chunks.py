#!/usr/bin/env python3
"""
merge_reco_chunks.py

Merge chunk CSVs from a Condor reco run into one file.

Usage:
    python3 ~/dmice/merge_reco_chunks.py --out-dir DIR [--output PATH]
"""

import os, sys, glob, csv, argparse

parser = argparse.ArgumentParser()
parser.add_argument("--out-dir",  required=True)
parser.add_argument("--output",   default=None, help="Merged CSV path (default: out-dir/merged.csv)")
args = parser.parse_args()

out_dir = os.path.expanduser(args.out_dir)
chunks  = sorted(glob.glob(os.path.join(out_dir, "chunk_*.csv")))

if not chunks:
    print(f"No chunk CSVs found in {out_dir}")
    sys.exit(1)

merged_path = args.output or os.path.join(out_dir, "merged.csv")

rows = []
fieldnames = None
missing = []

for path in chunks:
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        missing.append(os.path.basename(path))
        continue
    with open(path) as f:
        reader = csv.DictReader(f)
        if fieldnames is None:
            fieldnames = reader.fieldnames
        rows.extend(reader)

if missing:
    print(f"WARNING: {len(missing)} missing/empty chunks: {missing}")

# Sort by mc_energy_GeV to keep output deterministic
rows.sort(key=lambda r: float(r.get("mc_energy_GeV", 0)))

with open(merged_path, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fieldnames)
    w.writeheader()
    w.writerows(rows)

print(f"Merged {len(rows)} events from {len(chunks)-len(missing)} chunks")
print(f"Output: {merged_path}")
