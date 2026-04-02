#!/usr/bin/env python3
"""
blo_npz_to_i3.py

Convert BLO-format NPZ (from batch_dm_ice_targeted_sim.py or batch_dm_ice_sim.py)
to an .i3 file readable by Steamshovel and IceTray analysis scripts.

Coordinate conventions:
  - BLO NPZ dom_x/y/z: depth coordinates (m); z_icecube = z_depth + Z_OFFSET
  - BLO NPZ zenith_rad: Prometheus MOMENTUM convention (zenith > 90° = downgoing)
  - I3 direction: IceCube ANTI-MOMENTUM convention → zenith_i3 = π - zenith_blo

Run inside IceTray env-shell:
    /path/to/env-shell.sh python3 blo_npz_to_i3.py input.npz output.i3

Usage:
    python3 blo_npz_to_i3.py input.npz output.i3 [--geo icecube_with_dmice.geo]
"""

import sys
import os
import argparse
import numpy as np
from collections import defaultdict

from icecube import icetray, dataclasses, dataio, simclasses
from icecube.icetray import I3Units

# Depth-coordinate → IceCube-coordinate z offset (metres)
Z_OFFSET = 1948.07

# Default geo file (alongside this script)
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DEFAULT_GEO = os.path.join(_SCRIPT_DIR, "icecube_with_dmice.geo")


def load_geo(geo_path):
    """Parse a Prometheus .geo file → dict of (string, dom) → (x, y, z_icecube)."""
    doms = {}
    with open(geo_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 5:
                continue
            try:
                # Format: x y z_depth string_id dom_id
                x         = float(parts[0])
                y         = float(parts[1])
                z_depth   = float(parts[2])
                string_id = int(parts[3])
                dom_id    = int(parts[4])
                doms[(string_id, dom_id)] = (x, y, z_depth + Z_OFFSET)
            except ValueError:
                continue
    return doms


def blo_npz_to_i3(input_npz, output_i3, geo_path=DEFAULT_GEO, run_id=2000):
    print(f"[INFO] Loading: {input_npz}")
    d = np.load(input_npz, allow_pickle=True)

    n_events    = len(d["energy_GeV"])
    energy_GeV  = d["energy_GeV"]
    zenith_rad  = d["zenith_rad"]   # Prometheus momentum convention
    azimuth_rad = d["azimuth_rad"]
    det_id      = d["det_id"]   if "det_id"   in d else ["unknown"] * n_events
    dir_type    = d["dir_type"] if "dir_type" in d else ["unknown"] * n_events

    # Support both object-array format (numpy 2.x) and flat+offsets format (numpy 1.x compat)
    def load_ragged(key):
        if f"{key}_flat" in d:
            flat    = d[f"{key}_flat"]
            offsets = d[f"{key}_offsets"]
            return [flat[offsets[i]:offsets[i+1]] for i in range(n_events)]
        return d[key]

    dom_x      = load_ragged("dom_x")
    dom_y      = load_ragged("dom_y")
    dom_z      = load_ragged("dom_z")      # depth coordinates
    dom_t      = load_ragged("dom_t")      # ns
    dom_string = load_ragged("dom_string")
    dom_sensor = load_ragged("dom_sensor")

    print(f"[INFO] Events:   {n_events}")

    # ── Geometry ─────────────────────────────────────────────────────────────

    outfile = dataio.I3File(output_i3, "w")

    geo_frame = icetray.I3Frame(icetray.I3Frame.Geometry)
    geo = dataclasses.I3Geometry()

    if os.path.exists(geo_path):
        doms_geo = load_geo(geo_path)
        for (s, dom), (px, py, pz_ic) in doms_geo.items():
            omkey = icetray.OMKey(s, dom)
            omgeo = dataclasses.I3OMGeo()
            omgeo.position = dataclasses.I3Position(px, py, pz_ic)
            omgeo.omtype = dataclasses.I3OMGeo.IceCube
            geo.omgeo[omkey] = omgeo
        print(f"[INFO] Geometry: {len(doms_geo)} DOMs from {os.path.basename(geo_path)}")
    else:
        print(f"[WARN] Geo file not found: {geo_path} — I3Geometry will be empty")

    geo_frame["I3Geometry"] = geo
    outfile.push(geo_frame)

    # ── Events ───────────────────────────────────────────────────────────────

    for i in range(n_events):
        frame = icetray.I3Frame(icetray.I3Frame.Physics)

        # Event header
        header = dataclasses.I3EventHeader()
        header.run_id   = run_id
        header.event_id = i
        frame["I3EventHeader"] = header

        # MC truth — convert Prometheus momentum zenith → IceCube anti-momentum
        zen_ic  = float(np.pi - zenith_rad[i])
        azi_ic  = float(azimuth_rad[i])
        ene     = float(energy_GeV[i])

        mc_tree = dataclasses.I3MCTree()
        primary = dataclasses.I3Particle()
        primary.type          = dataclasses.I3Particle.MuMinus
        primary.location_type = dataclasses.I3Particle.InIce
        primary.shape         = dataclasses.I3Particle.InfiniteTrack
        primary.energy        = ene * I3Units.GeV
        primary.dir           = dataclasses.I3Direction(zen_ic, azi_ic)
        # Approximate vertex: DM-Ice depth z converted to IceCube z
        primary.pos           = dataclasses.I3Position(0.0, 0.0, 0.0)
        t0 = float(np.min(dom_t[i])) if len(dom_t[i]) > 0 else 0.0
        primary.time          = t0 * I3Units.ns
        mc_tree.add_primary(primary)
        frame["I3MCTree"] = mc_tree

        # Per-event hit arrays
        xs      = dom_x[i]
        ys      = dom_y[i]
        zs      = dom_z[i]          # depth coords
        ts      = dom_t[i]          # ns
        strings = dom_string[i].astype(int)
        sensors = dom_sensor[i].astype(int)

        # Build per-DOM hit lists (one entry per DOM, aggregated)
        dom_data = defaultdict(lambda: {"t": [], "x": 0.0, "y": 0.0, "z_ic": 0.0})
        for s, dom, x, y, z_dep, t in zip(strings, sensors, xs, ys, zs, ts):
            key = (int(s), int(dom))
            dom_data[key]["t"].append(float(t))
            dom_data[key]["x"]    = float(x)
            dom_data[key]["y"]    = float(y)
            dom_data[key]["z_ic"] = float(z_dep) + Z_OFFSET

        # I3MCPESeriesMap
        mcpe_map = simclasses.I3MCPESeriesMap()
        for (s, dom), info in dom_data.items():
            omkey     = icetray.OMKey(s, dom)
            pe_series = simclasses.I3MCPESeries()
            for t in sorted(info["t"]):
                pe       = simclasses.I3MCPE()
                pe.time  = t * I3Units.ns
                pe.npe   = 1
                pe_series.append(pe)
            mcpe_map[omkey] = pe_series
        frame["I3MCPESeriesMap"] = mcpe_map

        # InIcePulses (for linefit scripts)
        pulse_map = dataclasses.I3RecoPulseSeriesMap()
        for (s, dom), info in dom_data.items():
            omkey  = icetray.OMKey(s, dom)
            pulses = dataclasses.I3RecoPulseSeries()
            for t in sorted(info["t"]):
                p        = dataclasses.I3RecoPulse()
                p.time   = t  # ns
                p.charge = 1.0
                pulses.append(p)
            pulse_map[omkey] = pulses
        frame["InIcePulses"] = pulse_map

        # Metadata
        frame["BLO_DetId"]  = dataclasses.I3String(str(det_id[i]))
        frame["BLO_DirType"] = dataclasses.I3String(str(dir_type[i]))

        outfile.push(frame)

        if (i + 1) % 50 == 0 or (i + 1) == n_events:
            print(f"[INFO] Written {i + 1}/{n_events} events...")

    outfile.close()
    print(f"[INFO] Done → {output_i3}")
    print(f"[INFO] Open with: steamshovel {output_i3}")


def main():
    parser = argparse.ArgumentParser(description="Convert BLO NPZ to I3")
    parser.add_argument("input",          help="Input .npz file (from BLO batch sim)")
    parser.add_argument("output",         help="Output .i3 file")
    parser.add_argument("--geo",          default=DEFAULT_GEO, help="Path to .geo file")
    parser.add_argument("--run-id",       type=int, default=2000, help="I3 run ID (default: 2000)")
    args = parser.parse_args()
    blo_npz_to_i3(args.input, args.output, geo_path=args.geo, run_id=args.run_id)


if __name__ == "__main__":
    main()
