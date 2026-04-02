#!/usr/bin/env python3
"""
Targeted DM-Ice muon simulation using BLO (BlueLightOrchestra + PPC).

Simulates muons that physically pass through a DM-Ice NaI scintillator
detector (det1 or det2). Supports both downgoing and upgoing muons.
Injection geometry: start position back-projected 2 km from the DM-Ice
detector along the anti-momentum direction, so the muon is guaranteed
to traverse DM-Ice.

Physics:
  - Downgoing (zenith 130–170 deg): atmospheric muons, 1 TeV – 1 PeV, gamma=2
  - Upgoing   (zenith  10–50  deg): through-going,    100 GeV – 1 PeV, gamma=1
  - MuMinus (PDG 13)

DM-Ice positions (IceCube coords → converted to depth coords for BLO):
  det1: [31.25, -72.93,   -511.05] m  (string 87, DOM 1)
  det2: [-334.80, -424.50, -511.26] m  (string 88, DOM 1)

Output:
  NPZ: ~/dmice_work/output/blo_dmice_targeted_{det}_{direction}_{N}events.npz

Usage (on WARD):
    python ~/dmice/BLO/batch_dm_ice_targeted_sim.py
    python ~/dmice/BLO/batch_dm_ice_targeted_sim.py --nevents 50 --det 1 --direction down
    python ~/dmice/BLO/batch_dm_ice_targeted_sim.py --det both --direction both --no-gpu
"""

import os
import sys
import argparse
import numpy as np
import juliacall
from juliacall import Main as jl

# ── CLI ───────────────────────────────────────────────────────────────────────

parser = argparse.ArgumentParser(description="BLO targeted DM-Ice muon simulation")
parser.add_argument("--nevents",   type=int,   default=100,   help="Target accepted events (default: 100)")
parser.add_argument("--det",       type=str,   default="both", choices=["1", "2", "both"], help="DM-Ice detector (default: both)")
parser.add_argument("--direction", type=str,   default="both", choices=["up", "down", "both"], help="Muon direction (default: both)")
parser.add_argument("--outdir",    type=str,   default=os.path.expanduser("~/dmice_work/output/"), help="Output directory")
parser.add_argument("--no-gpu",    action="store_true", help="Use CPU PPC instead of GPU")
parser.add_argument("--seed",      type=int,   default=42,    help="RNG seed (default: 42)")
parser.add_argument("--logfile",   type=str,   default=None,  help="Tee progress to this file")
args = parser.parse_args()

USE_GPU    = not args.no_gpu
N_TARGET   = args.nevents
DET_ARG    = args.det
DIR_ARG    = args.direction
OUTPUT_DIR = args.outdir
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ── Logging ───────────────────────────────────────────────────────────────────

_logfile = open(args.logfile, "w", buffering=1) if args.logfile else None

def log(msg):
    print(msg, flush=True)
    if _logfile:
        print(msg, file=_logfile, flush=True)

# ── BLO environment ───────────────────────────────────────────────────────────

BLO_DIR = os.path.expanduser("~/.icevenv/BLO")

log("[INFO] Activating BLO Julia environment...")
jl.seval(f"""
using Pkg
Pkg.activate("{BLO_DIR}")
using BlueLightOrchestra
using BlueLightOrchestra.AstroParticleUnits
using BlueLightOrchestra.StaticArrays
using BlueLightOrchestra.Corpuscles
""")

BLO = jl.BlueLightOrchestra
jlx = getattr(jl, "*")

def make_particle(ene_GeV, pos_km, dir_xyz, pid=13, time_ns=0.0):
    """Create a BLO ParticleState from plain Python numbers."""
    return BLO.ParticleState(
        jlx(ene_GeV,  jl.GeV),
        jlx(pos_km,   jl.km),
        dir_xyz,
        jl.PDGID(pid),
        jlx(time_ns,  jl.ns),
    )

# ── Detector positions ────────────────────────────────────────────────────────

# IceCube coordinate → depth coordinate: depth_z = icecube_z - 1948.07 m
Z_OFFSET_M = 1948.07

DMICE_POS_M = {
    "det1": np.array([ 31.25,   -72.93,  -511.05]),
    "det2": np.array([-334.80, -424.50,  -511.26]),
}

# Convert to km depth coords (x, y unchanged; z shifted)
DMICE_POS_KM = {}
for det, pos in DMICE_POS_M.items():
    depth_z_km = (pos[2] - Z_OFFSET_M) / 1000.0
    DMICE_POS_KM[det] = np.array([pos[0] / 1000.0, pos[1] / 1000.0, depth_z_km])

# Which detectors to use
if DET_ARG == "both":
    det_list = ["det1", "det2"]
else:
    det_list = [f"det{DET_ARG}"]

# ── Simulation parameters ─────────────────────────────────────────────────────

BACKPROJECT_KM = 2.0    # back-project this far from DM-Ice along anti-momentum
PROP_KM        = 4.0    # total propagation distance (covers full detector volume)

# Direction-specific parameters
DIR_PARAMS = {
    "down": {
        "zen_min_deg": 130.0, "zen_max_deg": 170.0,   # downgoing in momentum convention
        "e_min_GeV":   1e3,   "e_max_GeV":   1e6,     # 1 TeV – 1 PeV
        "gamma":       2.0,                            # atmospheric spectrum
        "hit_cut":     50,                             # lower threshold (fewer side DOMs lit)
    },
    "up": {
        "zen_min_deg":  10.0, "zen_max_deg":  50.0,   # upgoing
        "e_min_GeV":   1e2,   "e_max_GeV":   1e6,     # 100 GeV – 1 PeV
        "gamma":       1.0,                            # log-flat
        "hit_cut":     200,
    },
}

if DIR_ARG == "both":
    dir_list = ["down", "up"]
else:
    dir_list = [DIR_ARG]

# ── Helper: sample energy from power-law ──────────────────────────────────────

def sample_energy(rng, e_min, e_max, gamma):
    """Sample energy from E^-gamma spectrum via inverse CDF."""
    if abs(gamma - 1.0) < 1e-6:
        # log-flat (gamma=1)
        return np.exp(rng.uniform(np.log(e_min), np.log(e_max)))
    g1 = 1.0 - gamma
    u = rng.uniform()
    return (e_min**g1 + u * (e_max**g1 - e_min**g1)) ** (1.0 / g1)

# ── Output filename ───────────────────────────────────────────────────────────

det_str = DET_ARG.replace("both", "det1det2")
dir_str = DIR_ARG
output_file = os.path.join(OUTPUT_DIR,
    f"blo_dmice_targeted_{det_str}_{dir_str}_{N_TARGET}events.npz")

# ── Storage ───────────────────────────────────────────────────────────────────

ev_energy_GeV  = []
ev_zenith_rad  = []
ev_azimuth_rad = []
ev_n_hits      = []
ev_n_doms      = []
ev_det_id      = []   # "det1" or "det2"
ev_dir_type    = []   # "up" or "down"
ev_dom_x       = []
ev_dom_y       = []
ev_dom_z       = []
ev_dom_t       = []
ev_dom_nhit    = []
ev_dom_str     = []
ev_dom_sen     = []

# ── Main loop ─────────────────────────────────────────────────────────────────

rng = np.random.default_rng(seed=args.seed)
n_simulated = 0
n_accepted  = 0

log(f"[INFO] Target:    {N_TARGET} accepted events")
log(f"[INFO] Detectors: {det_list}")
log(f"[INFO] Direction: {dir_list}")
log(f"[INFO] GPU:       {USE_GPU}")
log(f"[INFO] Output:    {output_file}")
log("")

while n_accepted < N_TARGET:

    # -- randomly pick detector and direction for this event
    det     = rng.choice(det_list)
    dir_typ = rng.choice(dir_list)
    p       = DIR_PARAMS[dir_typ]

    # -- sample energy
    ene_GeV = sample_energy(rng, p["e_min_GeV"], p["e_max_GeV"], p["gamma"])

    # -- sample direction (uniform in zenith range)
    zen_min = np.radians(p["zen_min_deg"])
    zen_max = np.radians(p["zen_max_deg"])
    zen = rng.uniform(zen_min, zen_max)
    azi = rng.uniform(0.0, 2.0 * np.pi)

    dx = np.sin(zen) * np.cos(azi)
    dy = np.sin(zen) * np.sin(azi)
    dz = np.cos(zen)   # <0 for downgoing (130–170 deg), >0 for upgoing

    # -- injection point: back-project from DM-Ice position
    dmice_km = DMICE_POS_KM[det]
    pos_km = dmice_km - BACKPROJECT_KM * np.array([dx, dy, dz])

    try:
        p_init     = make_particle(ene_GeV, pos_km.tolist(), [dx, dy, dz])
        ene_losses = BLO.propagate(p_init, jlx(PROP_KM, jl.km))
        hits       = BLO.run_ppc(p_init, ene_losses,
                                 suppress_error=True, use_gpu=USE_GPU)
        uhits      = BLO.process_hits(hits)
    except Exception as exc:
        log(f"  [WARN] sim {n_simulated+1} error: {exc}")
        n_simulated += 1
        continue

    n_simulated += 1

    nhits_arr  = np.array(uhits.nhits)
    total_hits = int(nhits_arr.sum())
    n_doms     = len(nhits_arr)

    if total_hits <= p["hit_cut"]:
        continue

    # -- accepted
    n_accepted += 1

    ev_energy_GeV.append(ene_GeV)
    ev_zenith_rad.append(zen)
    ev_azimuth_rad.append(azi)
    ev_n_hits.append(total_hits)
    ev_n_doms.append(n_doms)
    ev_det_id.append(det)
    ev_dir_type.append(dir_typ)

    _ustrip = jl.seval("(unit, arr) -> Float64.(ustrip.(unit, arr))")
    ev_dom_x.append(np.array(_ustrip(jl.m,  uhits.pos.x)))
    ev_dom_y.append(np.array(_ustrip(jl.m,  uhits.pos.y)))
    ev_dom_z.append(np.array(_ustrip(jl.m,  uhits.pos.z)))
    ev_dom_t.append(np.array(_ustrip(jl.ns, uhits.time)))
    ev_dom_nhit.append(nhits_arr)
    ev_dom_str.append(np.array(uhits.string_id))
    ev_dom_sen.append(np.array(uhits.sensor_id))

    log(f"  [{n_accepted:>3}/{N_TARGET}]  "
          f"tried={n_simulated:<5}  "
          f"{det}  {dir_typ:>4}  "
          f"E={ene_GeV/1e3:>7.2f} TeV  "
          f"zen={np.degrees(zen):>6.1f} deg  "
          f"hits={total_hits:>5}  doms={n_doms:>4}")

# ── Save ──────────────────────────────────────────────────────────────────────

log("")
efficiency = 100.0 * n_accepted / n_simulated
log(f"[INFO] {n_accepted}/{n_simulated} accepted  (efficiency {efficiency:.1f}%)")
log(f"[INFO] Saving to {output_file}")

np.savez(
    output_file,
    energy_GeV  = np.array(ev_energy_GeV),
    zenith_rad  = np.array(ev_zenith_rad),
    azimuth_rad = np.array(ev_azimuth_rad),
    n_hits      = np.array(ev_n_hits,   dtype=int),
    n_doms      = np.array(ev_n_doms,   dtype=int),
    det_id      = np.array(ev_det_id),
    dir_type    = np.array(ev_dir_type),
    dom_x       = np.array(ev_dom_x,    dtype=object),
    dom_y       = np.array(ev_dom_y,    dtype=object),
    dom_z       = np.array(ev_dom_z,    dtype=object),
    dom_t       = np.array(ev_dom_t,    dtype=object),
    dom_nhits   = np.array(ev_dom_nhit, dtype=object),
    dom_string  = np.array(ev_dom_str,  dtype=object),
    dom_sensor  = np.array(ev_dom_sen,  dtype=object),
)

log(f"[INFO] Done — {n_accepted} events saved to {output_file}")
log("")
log("To reload:")
log(f"  import numpy as np")
log(f"  d = np.load('{output_file}', allow_pickle=True)")
log(f"  d['energy_GeV']   # shape ({n_accepted},)")
log(f"  d['det_id']       # 'det1' or 'det2' per event")
log(f"  d['dir_type']     # 'up' or 'down' per event")
log(f"  d['dom_x'][0]     # x positions of hit DOMs in event 0  [m]")
