#!/usr/bin/env python3
"""
parquet_to_npz.py  (Step 1 of 2)

Extracts Prometheus parquet photon data to an intermediate .npz file.
Run with system python3 (requires pyarrow, no IceTray needed).

Usage:
    python3 parquet_to_npz.py <input.parquet> <output.npz>
"""

import sys
import argparse
import json
import numpy as np
import pyarrow.parquet as pq


def parquet_to_npz(input_parquet: str, output_npz: str) -> None:
    print(f"[INFO] Reading: {input_parquet}")
    table = pq.read_table(input_parquet)
    n_events = len(table)
    print(f"[INFO] Events: {n_events}")

    # Extract mc_truth columns (one value per event)
    mc = table.column("mc_truth").combine_chunks()
    photons = table.column("photons").combine_chunks()

    # Pull per-event truth arrays
    def get_mc_field(field):
        try:
            arr = mc.field(field)
            # Some fields may themselves be ChunkedArrays after combine_chunks
            if hasattr(arr, "combine_chunks"):
                arr = arr.combine_chunks()
            return np.array(arr.to_pylist(), dtype=float)
        except Exception:
            return np.zeros(n_events)

    energy  = get_mc_field("initial_state_energy")
    zenith  = get_mc_field("initial_state_zenith")
    azimuth = get_mc_field("initial_state_azimuth")
    x       = get_mc_field("initial_state_x")
    y       = get_mc_field("initial_state_y")
    z       = get_mc_field("initial_state_z")

    # Try to get run number from parquet metadata
    run_number = 1001
    try:
        meta = table.schema.metadata
        if meta and b"config_prometheus" in meta:
            cfg = json.loads(meta[b"config_prometheus"])
            run_number = int(cfg["run"].get("run number", 1001))
    except Exception:
        pass

    arrays = dict(
        n_events=np.array(n_events),
        run_number=np.array(run_number),
        energy=energy,
        zenith=zenith,
        azimuth=azimuth,
        x=x, y=y, z=z,
    )

    # Per-event photon hit arrays
    string_ids_list = photons.field("string_id").to_pylist()
    sensor_ids_list = photons.field("sensor_id").to_pylist()
    t_list          = photons.field("t").to_pylist()
    sensor_pos_x_list = photons.field("sensor_pos_x").to_pylist()
    sensor_pos_y_list = photons.field("sensor_pos_y").to_pylist()
    sensor_pos_z_list = photons.field("sensor_pos_z").to_pylist()

    # Collect unique sensor positions: (string_id, sensor_id) -> (x, y, z)
    sensor_positions = {}

    total_hits = 0
    for i in range(n_events):
        sid = np.array(string_ids_list[i], dtype=np.int32)
        did = np.array(sensor_ids_list[i], dtype=np.int32)
        t   = np.array(t_list[i],          dtype=np.float64)
        arrays[f"ev_{i}_string_id"] = sid
        arrays[f"ev_{i}_sensor_id"] = did
        arrays[f"ev_{i}_t"]         = t
        total_hits += len(sid)

        for s, d, px, py, pz in zip(
            string_ids_list[i], sensor_ids_list[i],
            sensor_pos_x_list[i], sensor_pos_y_list[i], sensor_pos_z_list[i]
        ):
            key = (int(s), int(d))
            if key not in sensor_positions:
                sensor_positions[key] = (float(px), float(py), float(pz))

    print(f"[INFO] Total photon hits: {total_hits}")

    # Save unique sensor positions for I3Geometry construction
    if sensor_positions:
        keys = sorted(sensor_positions.keys())
        arrays["sensor_string_ids"] = np.array([k[0] for k in keys], dtype=np.int32)
        arrays["sensor_dom_ids"]    = np.array([k[1] for k in keys], dtype=np.int32)
        arrays["sensor_pos_x"]      = np.array([sensor_positions[k][0] for k in keys])
        arrays["sensor_pos_y"]      = np.array([sensor_positions[k][1] for k in keys])
        arrays["sensor_pos_z"]      = np.array([sensor_positions[k][2] for k in keys])
        print(f"[INFO] Unique sensors: {len(keys)}")

    np.savez(output_npz, **arrays)
    print(f"[INFO] Saved: {output_npz}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input",  help="Input .parquet file from Prometheus")
    parser.add_argument("output", help="Output .npz file")
    args = parser.parse_args()
    parquet_to_npz(args.input, args.output)


if __name__ == "__main__":
    main()
