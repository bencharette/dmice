#!/usr/bin/env python3
"""
run_itermpe_events.py

Runs IterMPE (MPE + charge cap, 3 iterations) on selected events from a BLO NPZ,
and saves the reconstructed direction vectors alongside MC truth to a CSV.

Usage (on Cobalt with IceTray env):
    python3 run_itermpe_events.py [--npz PATH] [--out CSV] [--events 612 1471 ...]
"""

import os, sys, csv, math, argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("--npz", default=os.path.expanduser(
    "~/dmice_work/output/muons_binned_5bins_1000pbin_repacked.npz"))
parser.add_argument("--out", default=os.path.expanduser(
    "~/dmice_work/output/itermpe_events.csv"))
parser.add_argument("--events", type=int, nargs="+", default=None,
    help="Event indices to process (default: top nhits event per bin)")
args = parser.parse_args()

GEO_FILE = os.path.expanduser(
    "~/dmice/BlueLightOrchestra.jl/resources/geofiles/icecube_with_dmice.geo")
Z_OFFSET  = 1948.07
C_M_NS    = 0.2998
MAX_CHARGE = 5.0

from icecube import icetray, dataclasses, dataio, simclasses
from icecube import linefit, lilliput
import icecube.lilliput.segments
from icecube.icetray import I3Units, I3Tray

# ── Geometry ──────────────────────────────────────────────────────────────────

def load_geo(path):
    doms = {}
    with open(path) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 5:
                continue
            try:
                x, y, z_dep = float(parts[0]), float(parts[1]), float(parts[2])
                s, dom = int(parts[3]), int(parts[4])
                doms[(s, dom)] = (x, y, z_dep + Z_OFFSET)
            except ValueError:
                continue
    return doms

geo_doms = load_geo(GEO_FILE)
geo_obj  = dataclasses.I3Geometry()
for (s, dom), (px, py, pz) in geo_doms.items():
    omkey = icetray.OMKey(s, dom)
    omgeo = dataclasses.I3OMGeo()
    omgeo.position = dataclasses.I3Position(px, py, pz)
    omgeo.omtype   = dataclasses.I3OMGeo.IceCube
    geo_obj.omgeo[omkey] = omgeo

# ── Load NPZ ──────────────────────────────────────────────────────────────────

d = np.load(args.npz, allow_pickle=True)
N_TOTAL = len(d["energy_GeV"])

# Select events: specified indices, or top nhits per bin
if args.events:
    indices = args.events
else:
    bin_ids = d["bin_id"]
    n_hits  = d["n_hits"]
    indices = []
    for b in np.unique(bin_ids):
        mask = np.where(bin_ids == b)[0]
        top  = mask[np.argsort(n_hits[mask])[::-1][:1]]
        indices.extend(top.tolist())
    indices = sorted(indices)

print(f"Processing {len(indices)} events: {indices}")

def load_ragged(key, i):
    if f"{key}_flat" in d:
        flat    = d[f"{key}_flat"]
        offsets = d[f"{key}_offsets"]
        return flat[offsets[i]:offsets[i+1]]
    return d[key][i]

# ── Per-event IceTray run ─────────────────────────────────────────────────────

rows = []

def ang_err_deg(t1, t2):
    dot = max(-1.0, min(1.0, t1[0]*t2[0] + t1[1]*t2[1] + t1[2]*t2[2]))
    return math.degrees(math.acos(abs(dot)))

for ev_idx in indices:
    zen = float(d["zenith_rad"][ev_idx])
    azi = float(d["azimuth_rad"][ev_idx])
    ene = float(d["energy_GeV"][ev_idx])

    dx_mc = math.sin(zen) * math.cos(azi)
    dy_mc = math.sin(zen) * math.sin(azi)
    dz_mc = math.cos(zen)
    mc_dir = (dx_mc, dy_mc, dz_mc)

    dom_x   = np.array(load_ragged("dom_x",      ev_idx))
    dom_y   = np.array(load_ragged("dom_y",      ev_idx))
    dom_z   = np.array(load_ragged("dom_z",      ev_idx))
    dom_t   = np.array(load_ragged("dom_t",      ev_idx))
    dom_nh  = np.array(load_ragged("dom_nhits",  ev_idx))
    dom_str = np.array(load_ragged("dom_string", ev_idx), dtype=int)
    dom_sen = np.array(load_ragged("dom_sensor", ev_idx), dtype=int)

    if len(dom_x) < 4:
        print(f"  [ev {ev_idx}] skipped — too few DOMs")
        continue

    # Build single-event tray
    result = {}

    class Injector(icetray.I3Module):
        def __init__(self, ctx):
            super().__init__(ctx)
            self.done = False
        def Configure(self): pass
        def Process(self):
            if self.done:
                self.RequestSuspension(); return
            gframe = icetray.I3Frame(icetray.I3Frame.Geometry)
            gframe["I3Geometry"] = geo_obj
            self.PushFrame(gframe)

            frame = icetray.I3Frame(icetray.I3Frame.Physics)
            hdr = dataclasses.I3EventHeader()
            hdr.run_id = 1; hdr.event_id = ev_idx
            frame["I3EventHeader"] = hdr

            primary = dataclasses.I3Particle()
            primary.type       = dataclasses.I3Particle.MuMinus
            primary.shape      = dataclasses.I3Particle.InfiniteTrack
            primary.energy     = ene * I3Units.GeV
            primary.dir        = dataclasses.I3Direction(dx_mc, dy_mc, dz_mc)
            primary.fit_status = dataclasses.I3Particle.FitStatus.OK
            mc_tree = dataclasses.I3MCTree()
            mc_tree.add_primary(primary)
            frame["I3MCTree"] = mc_tree
            frame["MCTruth"]  = primary

            pulse_map = dataclasses.I3RecoPulseSeriesMap()
            for j in range(len(dom_x)):
                omkey = icetray.OMKey(int(dom_str[j]), int(dom_sen[j]))
                ps = dataclasses.I3RecoPulseSeries()
                p  = dataclasses.I3RecoPulse()
                p.time   = float(dom_t[j])
                p.charge = float(dom_nh[j])
                ps.append(p)
                pulse_map[omkey] = ps
            frame["InIcePulses"] = pulse_map

            # Charge-capped pulses for IterMPE
            capped = dataclasses.I3RecoPulseSeriesMap()
            for omk, plist in pulse_map:
                new_ps = dataclasses.I3RecoPulseSeries()
                for p in plist:
                    np_ = dataclasses.I3RecoPulse()
                    np_.time   = p.time
                    np_.charge = min(p.charge, MAX_CHARGE)
                    new_ps.append(np_)
                capped[omk] = new_ps
            frame["InIcePulsesCapped"] = capped

            self.PushFrame(frame)
            self.done = True

    def extract(frame):
        if "IterMPE" not in frame:
            return
        p = frame["IterMPE"]
        if p.fit_status != dataclasses.I3Particle.FitStatus.OK:
            result["ok"] = False
            return
        result["ok"]  = True
        result["dx"]  = p.dir.x
        result["dy"]  = p.dir.y
        result["dz"]  = p.dir.z
        result["zen"] = p.dir.zenith
        result["azi"] = p.dir.azimuth

    tray = I3Tray()
    tray.Add(Injector)
    tray.Add("I3LineFit",
        Name            = "LineFit",
        InputRecoPulses = "InIcePulses",
        AmpWeightPower  = 1.0,
    )
    tray.Add(icecube.lilliput.segments.I3IterativePandelFitter,
        fitname      = "IterMPE",
        domllh       = "MPE",
        pulses       = "InIcePulsesCapped",
        seeds        = ["LineFit"],
        n_iterations = 3,
        If           = lambda f: "LineFit" in f,
    )
    tray.Add(extract, Streams=[icetray.I3Frame.Physics])
    tray.Execute()
    tray.Finish()

    if result.get("ok"):
        reco_dir = (result["dx"], result["dy"], result["dz"])
        err = ang_err_deg(mc_dir, reco_dir)
        print(f"  [ev {ev_idx}]  E={ene/1e3:.2f} TeV  IterMPE ang_err={err:.2f}°")
        rows.append(dict(
            ev_idx      = ev_idx,
            mc_energy_GeV = ene,
            mc_zen_rad  = zen,
            mc_azi_rad  = azi,
            mc_dx       = dx_mc,
            mc_dy       = dy_mc,
            mc_dz       = dz_mc,
            itermpe_dx  = result["dx"],
            itermpe_dy  = result["dy"],
            itermpe_dz  = result["dz"],
            itermpe_zen = result["zen"],
            itermpe_azi = result["azi"],
            ang_err_deg = err,
        ))
    else:
        print(f"  [ev {ev_idx}]  IterMPE failed")

os.makedirs(os.path.dirname(args.out), exist_ok=True)
with open(args.out, "w", newline="") as f:
    if rows:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
print(f"Saved: {args.out}  ({len(rows)} events)")
