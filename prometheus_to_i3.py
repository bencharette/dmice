#!/usr/bin/env python3
"""
prometheus_to_i3.py  (Step 2 of 2)

Reads intermediate .npz file (produced by parquet_to_npz.py) and writes
an .i3 file readable by Steamshovel.

Run inside IceTray env-shell:
    python prometheus_to_i3.py <input.npz> <output.i3>
"""

import sys
import argparse
import numpy as np
from collections import defaultdict

from icecube import icetray, dataclasses, dataio, simclasses
from icecube.icetray import I3Units

# Depth-coordinate → IceCube-coordinate offset (metres)
Z_OFFSET = 1948.07


def npz_to_i3(input_npz: str, output_i3: str) -> None:
    print(f"[INFO] Loading: {input_npz}")
    data = np.load(input_npz, allow_pickle=True)

    n_events = int(data["n_events"])
    print(f"[INFO] Events: {n_events}")

    outfile = dataio.I3File(output_i3, "w")

    # Build I3Geometry from sensor positions saved by parquet_to_npz.py
    geo_frame = icetray.I3Frame(icetray.I3Frame.Geometry)
    if "sensor_string_ids" in data:
        geo = dataclasses.I3Geometry()
        string_ids_geo = data["sensor_string_ids"]
        dom_ids_geo    = data["sensor_dom_ids"]
        pos_x          = data["sensor_pos_x"]
        pos_y          = data["sensor_pos_y"]
        pos_z          = data["sensor_pos_z"]
        for s, d, px, py, pz in zip(string_ids_geo, dom_ids_geo, pos_x, pos_y, pos_z):
            omkey = icetray.OMKey(int(s), int(d))
            omgeo = dataclasses.I3OMGeo()
            omgeo.position = dataclasses.I3Position(
                float(px), float(py), float(pz) + Z_OFFSET
            )
            omgeo.omtype = dataclasses.I3OMGeo.IceCube
            geo.omgeo[omkey] = omgeo
        geo_frame['I3Geometry'] = geo
        print(f"[INFO] Geometry: {len(string_ids_geo)} DOMs")
    else:
        print("[WARN] No sensor positions in npz — I3Geometry will be empty")
    outfile.push(geo_frame)

    for i in range(n_events):
        frame = icetray.I3Frame(icetray.I3Frame.Physics)

        # Event header
        header = dataclasses.I3EventHeader()
        header.event_id = i
        header.run_id = int(data["run_number"]) if "run_number" in data else 1001
        frame["I3EventHeader"] = header

        # Primary particle from mc_truth
        mc_tree = dataclasses.I3MCTree()
        primary = dataclasses.I3Particle()
        primary.type = dataclasses.I3Particle.MuMinus
        primary.location_type = dataclasses.I3Particle.InIce
        primary.shape = dataclasses.I3Particle.InfiniteTrack

        energy = float(data["energy"][i]) if "energy" in data else 1e4
        zenith = float(data["zenith"][i]) if "zenith" in data else 0.0
        azimuth = float(data["azimuth"][i]) if "azimuth" in data else 0.0
        x = float(data["x"][i]) if "x" in data else 0.0
        y = float(data["y"][i]) if "y" in data else 0.0
        z = float(data["z"][i]) if "z" in data else 0.0

        primary.energy = energy * I3Units.GeV
        primary.dir = dataclasses.I3Direction(np.pi - zenith, azimuth)
        primary.pos = dataclasses.I3Position(x * I3Units.m, y * I3Units.m,
                                              (z + Z_OFFSET) * I3Units.m)

        string_ids = data[f"ev_{i}_string_id"]
        sensor_ids = data[f"ev_{i}_sensor_id"]
        times      = data[f"ev_{i}_t"]

        primary.time = float(np.min(times)) if len(times) > 0 else 0.0

        mc_tree.add_primary(primary)
        frame["I3MCTree"] = mc_tree

        # Photon hits → I3MCPESeriesMap
        mcpe_map = simclasses.I3MCPESeriesMap()

        dom_hits = defaultdict(list)
        for s, d, t in zip(string_ids, sensor_ids, times):
            dom_hits[(int(s), int(d))].append(float(t))

        for (string_id, dom_id), hit_times in dom_hits.items():
            omkey = icetray.OMKey(string_id, dom_id)
            pe_series = simclasses.I3MCPESeries()
            for t in sorted(hit_times):
                pe = simclasses.I3MCPE()
                pe.time = t * I3Units.ns
                pe.npe = 1
                pe_series.append(pe)
            mcpe_map[omkey] = pe_series

        frame["I3MCPESeriesMap"] = mcpe_map

        # Also write as InIcePulses so linefit scripts can find them
        pulse_map = dataclasses.I3RecoPulseSeriesMap()
        for (string_id, dom_id), hit_times in dom_hits.items():
            omkey = icetray.OMKey(string_id, dom_id)
            pulses = dataclasses.I3RecoPulseSeries()
            for t in sorted(hit_times):
                p = dataclasses.I3RecoPulse()
                p.time = t  # already in ns
                p.charge = 1.0
                pulses.append(p)
            pulse_map[omkey] = pulses
        frame["InIcePulses"] = pulse_map

        outfile.push(frame)

        if (i + 1) % 10 == 0:
            print(f"[INFO] Written {i + 1}/{n_events} events...")

    outfile.close()
    print(f"[INFO] Done! Output: {output_i3}")
    print(f"[INFO] Open with: steamshovel {output_i3}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="Input .npz file (from parquet_to_npz.py)")
    parser.add_argument("output", help="Output .i3 file")
    args = parser.parse_args()
    npz_to_i3(args.input, args.output)


if __name__ == "__main__":
    main()
