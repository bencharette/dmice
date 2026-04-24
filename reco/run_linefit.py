#!/usr/bin/env python3
"""
run_linefit.py

Post-processing script: reads a BLO i3 file and writes three I3Particle
reconstruction keys into every Physics frame:

    ICLineFit        — analytic IC-only LineFit
    DMIcePivotLineFit — DM-Ice Pivot LineFit (requires I3MCTree for hit time)
    MCTruthMuon      — primary muon copied from I3MCTree (for easy comparison)

All other frames (Geometry, DAQ, etc.) are passed through unchanged.

Usage (inside IceTray env-shell):
    python run_linefit.py -i input.i3.zst -o output_linefit.i3.zst

If -o is omitted, appends '_linefit' before the first dot in the input filename.
"""

import os
import sys
import math
import argparse

try:
    from icecube import icetray, dataclasses, dataio, simclasses
    from icecube.icetray import I3Units
except ImportError:
    sys.exit("ERROR: Load IceTray environment first.")

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

C_M_NS = 0.2998  # speed of light in m/ns

DMICE_POS = {
    "det1": (31.25,   -72.93,  -511.05),
    "det2": (-334.80, -424.50, -511.26),
}

PULSE_KEYS = [
    "SRTInIcePulses", "SplitInIcePulses", "InIcePulses",
    "TWCMuonPulseSeriesReco", "OfflinePulses",
    "OnlineL2_CleanedMuonPulses", "SplitUncleanedInIcePulses",
    "UncleanedInIcePulses",
]

# ---------------------------------------------------------------------------
# Pure-Python helpers (no numpy — matches steamshovel artist logic exactly)
# ---------------------------------------------------------------------------

def _wm(vals, ws):
    W = sum(ws)
    return sum(v * w for v, w in zip(vals, ws)) / W if W else 0.0

def _dot(a, b):
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]

def _sub(a, b):
    return (a[0]-b[0], a[1]-b[1], a[2]-b[2])

def _scale(s, a):
    return (s*a[0], s*a[1], s*a[2])

def _norm(a):
    return math.sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2])


def _get_dom_positions(geo_frame):
    dom_pos = {}
    if geo_frame is None or "I3Geometry" not in geo_frame:
        return dom_pos
    for omkey, omgeo in geo_frame["I3Geometry"].omgeo:
        p = omgeo.position
        dom_pos[omkey] = (p.x, p.y, p.z)
    return dom_pos


def _extract_hits(frame, dom_pos):
    pulse_map = None
    for key in PULSE_KEYS:
        if key in frame:
            raw = frame[key]
            pulse_map = raw.apply(frame) if hasattr(raw, "apply") else raw
            break
    if pulse_map is None:
        return None

    xs, ys, zs, ts, ws = [], [], [], [], []
    for omkey, pulses in pulse_map:
        if omkey not in dom_pos:
            continue
        px, py, pz = dom_pos[omkey]
        for p in pulses:
            c = p.charge if p.charge > 0 else 1.0
            xs.append(px); ys.append(py); zs.append(pz)
            ts.append(p.time); ws.append(c)

    return (xs, ys, zs, ts, ws) if xs else None


def _ic_linefit(xs, ys, zs, ts, ws):
    W = sum(ws)
    if not W:
        return None
    cx, cy, cz = _wm(xs, ws), _wm(ys, ws), _wm(zs, ws)
    tb = _wm(ts, ws)
    dts = [t - tb for t in ts]
    den = sum(w * dt * dt for w, dt in zip(ws, dts))
    if not den:
        return None
    vx = sum(w * dt * (x - cx) for w, dt, x in zip(ws, dts, xs)) / den
    vy = sum(w * dt * (y - cy) for w, dt, y in zip(ws, dts, ys)) / den
    vz = sum(w * dt * (z - cz) for w, dt, z in zip(ws, dts, zs)) / den
    spd = math.sqrt(vx*vx + vy*vy + vz*vz)
    if not spd:
        return None
    return dict(dx=vx/spd, dy=vy/spd, dz=vz/spd, speed=spd, cx=cx, cy=cy, cz=cz, t=tb)


def _pivot_linefit(xs, ys, zs, ts, ws, dm_pos, mc_dir):
    cx, cy, cz = _wm(xs, ws), _wm(ys, ws), _wm(zs, ws)
    tb = _wm(ts, ws)
    d_proj = _dot(_sub(dm_pos, (cx, cy, cz)), mc_dir)
    t_dm = tb + d_proj / C_M_NS

    dts  = [t - t_dm for t in ts]
    drxs = [x - dm_pos[0] for x in xs]
    drys = [y - dm_pos[1] for y in ys]
    drzs = [z - dm_pos[2] for z in zs]
    den = sum(w * dt * dt for w, dt in zip(ws, dts))
    if not den:
        return None
    vx = sum(w * dt * dr for w, dt, dr in zip(ws, dts, drxs)) / den
    vy = sum(w * dt * dr for w, dt, dr in zip(ws, dts, drys)) / den
    vz = sum(w * dt * dr for w, dt, dr in zip(ws, dts, drzs)) / den
    spd = math.sqrt(vx*vx + vy*vy + vz*vz)
    if not spd:
        return None
    return dict(dx=vx/spd, dy=vy/spd, dz=vz/spd, speed=spd, cx=cx, cy=cy, cz=cz, t=tb)


def _get_primary_muon(frame):
    for key in ("I3MCTree", "I3MCTree_preMuonProp"):
        if key not in frame:
            continue
        tree = frame[key]
        best, best_e = None, 0.0
        for p in tree:
            if p.type in (dataclasses.I3Particle.MuMinus, dataclasses.I3Particle.MuPlus):
                if p.energy > best_e:
                    best, best_e = p, p.energy
        if best is not None:
            return best
        prims = tree.primaries
        if prims:
            return prims[0]
    return None


def _make_particle(fit, dm_name=None):
    """Package a fit result dict as an I3Particle (InfiniteTrack)."""
    p = dataclasses.I3Particle()
    p.shape      = dataclasses.I3Particle.InfiniteTrack
    p.fit_status = dataclasses.I3Particle.OK
    p.pos        = dataclasses.I3Position(fit["cx"], fit["cy"], fit["cz"])
    p.dir        = dataclasses.I3Direction(fit["dx"], fit["dy"], fit["dz"])
    p.time       = fit["t"] * I3Units.ns
    p.speed      = fit["speed"] * I3Units.m / I3Units.ns
    return p


# ---------------------------------------------------------------------------
# Main processing loop
# ---------------------------------------------------------------------------

def process(input_path, output_path):
    print("Input:  {}".format(input_path))
    print("Output: {}".format(output_path))

    infile  = dataio.I3File(input_path, "r")
    outfile = dataio.I3File(output_path, "w")

    dom_pos   = {}
    n_total   = 0
    n_ic      = 0
    n_pivot   = 0

    while infile.more():
        frame = infile.pop_frame()

        # Cache geometry from the Geometry frame
        if frame.Stop == icetray.I3Frame.Geometry:
            dom_pos = _get_dom_positions(frame)
            print("  Geometry: {} DOMs".format(len(dom_pos)))
            outfile.push(frame)
            continue

        # Pass non-Physics frames through unchanged
        if frame.Stop != icetray.I3Frame.Physics:
            outfile.push(frame)
            continue

        n_total += 1

        # ── MC truth muon ────────────────────────────────────────────────
        muon = _get_primary_muon(frame)
        if muon is not None:
            # Copy as a standalone key for easy steamshovel access
            mc_copy = dataclasses.I3Particle(muon)
            frame["MCTruthMuon"] = mc_copy
            mc_dir = (muon.dir.x, muon.dir.y, muon.dir.z)
        else:
            mc_dir = None

        # ── IC-only LineFit ──────────────────────────────────────────────
        hits = _extract_hits(frame, dom_pos)
        if hits is not None and len(hits[0]) >= 4:
            xs, ys, zs, ts, ws = hits
            ic_fit = _ic_linefit(xs, ys, zs, ts, ws)
            if ic_fit is not None:
                frame["ICLineFit"] = _make_particle(ic_fit)
                n_ic += 1

                # ── DM-Ice Pivot LineFit ─────────────────────────────────
                if mc_dir is not None:
                    # Select DM-Ice detector from BLO tag or closest approach
                    if "BLO_DetId" in frame:
                        tag = str(frame["BLO_DetId"].value)
                        dm_name = tag if tag in DMICE_POS else "det1"
                    else:
                        def ca(dp):
                            cx = ic_fit["cx"]; cy = ic_fit["cy"]; cz = ic_fit["cz"]
                            d = _sub(dp, (cx, cy, cz))
                            proj = _dot(d, mc_dir)
                            return _norm(_sub(d, _scale(proj, mc_dir)))
                        dm_name = "det1" if ca(DMICE_POS["det1"]) <= ca(DMICE_POS["det2"]) else "det2"

                    piv_fit = _pivot_linefit(xs, ys, zs, ts, ws, DMICE_POS[dm_name], mc_dir)
                    if piv_fit is not None:
                        frame["DMIcePivotLineFit"] = _make_particle(piv_fit)
                        frame["DMIceClosestDet"]   = dataclasses.I3String(dm_name)
                        n_pivot += 1

        outfile.push(frame)

        if n_total % 100 == 0:
            print("  {} frames processed...".format(n_total))

    infile.close()
    outfile.close()

    print("\nDone.")
    print("  Physics frames:         {}".format(n_total))
    print("  ICLineFit written:      {}".format(n_ic))
    print("  DMIcePivotLineFit written: {}".format(n_pivot))
    print("Output: {}".format(output_path))


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _default_output(input_path):
    base = os.path.basename(input_path)
    dot  = base.find(".")
    stem = base[:dot] if dot != -1 else base
    ext  = base[dot:] if dot != -1 else ""
    return os.path.join(os.path.dirname(input_path), stem + "_linefit" + ext)


def main():
    parser = argparse.ArgumentParser(
        description="Add ICLineFit, DMIcePivotLineFit, MCTruthMuon to a BLO i3 file")
    parser.add_argument("-i", "--input",  required=True,  help="Input .i3 file")
    parser.add_argument("-o", "--output", default=None,   help="Output .i3 file (default: input_linefit.*)")
    args = parser.parse_args()

    output = args.output or _default_output(args.input)
    if os.path.abspath(output) == os.path.abspath(args.input):
        sys.exit("ERROR: output path is the same as input — aborting.")

    process(args.input, output)


if __name__ == "__main__":
    main()
