#!/usr/bin/env python3
"""
repack_npz.py

Re-saves a BLO NPZ file using flat arrays + offsets instead of numpy object
arrays. The resulting file is compatible with numpy 1.x (no pickle needed).

Run on WARD (Python 3.14 / numpy 2.x) after simulation:
    python3 repack_npz.py input.npz output_repacked.npz

Usage:
    python3 repack_npz.py blo_dmice_targeted_det1det2_both_1000events.npz \
                          blo_dmice_targeted_det1det2_both_1000events_repacked.npz
"""

import sys
import numpy as np

RAGGED_KEYS = ["dom_x", "dom_y", "dom_z", "dom_t", "dom_nhits", "dom_string", "dom_sensor"]

def repack(input_path, output_path):
    print(f"[INFO] Loading {input_path}")
    d = np.load(input_path, allow_pickle=True)

    out = {}

    # Copy scalar/string arrays as-is
    for key in d.files:
        if key not in RAGGED_KEYS:
            out[key] = d[key]
            print(f"  copied: {key}  shape={d[key].shape}")

    # Flatten ragged arrays + save offsets
    for key in RAGGED_KEYS:
        if key not in d.files:
            continue
        arr = d[key]  # object array of variable-length arrays
        flat = np.concatenate([np.asarray(a, dtype=np.float64 if key not in ("dom_string", "dom_sensor", "dom_nhits") else np.int32) for a in arr])
        offsets = np.zeros(len(arr) + 1, dtype=np.int64)
        for i, a in enumerate(arr):
            offsets[i + 1] = offsets[i] + len(a)
        out[f"{key}_flat"]    = flat
        out[f"{key}_offsets"] = offsets
        print(f"  repacked: {key}  flat={flat.shape}  offsets={offsets.shape}")

    np.savez(output_path, **out)
    print(f"[INFO] Saved → {output_path}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python3 repack_npz.py input.npz output.npz")
        sys.exit(1)
    repack(sys.argv[1], sys.argv[2])
