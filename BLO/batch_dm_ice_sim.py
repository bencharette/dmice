#!/usr/bin/env python3
"""
Batch DM-Ice muon simulation using BLO (BlueLightOrchestra + PPC GPU).

Simulates upgoing through-going muons and collects 100 events with
more than 200 total DOM hits.

Injection parameters match simulate_dm_ice_through.py:
  - MuMinus (PDG 13), ranged/through-going
  - Energy: 100 GeV – 1 PeV, log-flat (gamma=1)
  - Direction: upgoing only, zenith 0–90 deg (Earth-filtered)
  - Tracks aimed through IceCube centre (0, 0, -1950 m)

Output: ~/dmice_work/output/blo_muons_200hits.npz

Usage (on WARD):
    python ~/dmice/batch_dm_ice_sim.py
"""

import os
import sys
import numpy as np
import juliacall
from juliacall import Main as jl

# ── BLO environment ───────────────────────────────────────────────────────────

BLO_DIR = os.path.expanduser("~/.icevenv/BLO")

print("[INFO] Activating BLO Julia environment...")
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

# ── Simulation parameters ─────────────────────────────────────────────────────

E_MIN_GEV   = 1e2          # 100 GeV
E_MAX_GEV   = 1e6          # 1 PeV
PROP_KM     = 3.0          # propagation distance — traverses all of IceCube
INJECT_Z_KM = -2.5         # start below detector (IceCube bottom ~ -2.45 km)
CENTRE_Z_KM = -1.950       # IceCube centre depth — tracks aimed here
HIT_CUT     = 200          # minimum total DOM hits to accept an event
N_TARGET    = 100          # stop after this many accepted events
USE_GPU     = True

rng = np.random.default_rng(seed=42)

# ── Output paths ──────────────────────────────────────────────────────────────

output_dir  = os.path.expanduser("~/dmice_work/output/")
output_file = os.path.join(output_dir, "blo_muons_200hits.npz")
os.makedirs(output_dir, exist_ok=True)

# ── Storage for accepted events ───────────────────────────────────────────────

# per-event metadata
ev_energy_GeV  = []
ev_zenith_rad  = []
ev_azimuth_rad = []
ev_n_hits      = []
ev_n_doms      = []

# per-event hit arrays (ragged — stored as object arrays at save time)
ev_dom_x    = []   # DOM x position [m]
ev_dom_y    = []
ev_dom_z    = []
ev_dom_t    = []   # first-hit time [ns]
ev_dom_nhit = []   # number of photon hits on each DOM
ev_dom_str  = []   # string ID
ev_dom_sen  = []   # sensor ID

# ── Main loop ─────────────────────────────────────────────────────────────────

n_simulated = 0
n_accepted  = 0

print(f"[INFO] Target: {N_TARGET} events with >{HIT_CUT} DOM hits")
print(f"[INFO] Energy: {E_MIN_GEV:.0f} GeV – {E_MAX_GEV:.0e} GeV  (log-flat, gamma=1)")
print(f"[INFO] Upgoing (zenith 0–90 deg), aimed through IceCube centre")
print(f"[INFO] GPU: {USE_GPU}")
print()

while n_accepted < N_TARGET:

    # -- sample energy: log-flat (gamma=1 spectrum)
    ene_GeV = np.exp(rng.uniform(np.log(E_MIN_GEV), np.log(E_MAX_GEV)))

    # -- sample direction: isotropic upgoing hemisphere
    #    uniform in cos(zenith) → cos_zen in [0, 1] (0=horizontal, 1=straight up)
    cos_zen = rng.uniform(0.0, 1.0)
    zen     = np.arccos(cos_zen)
    azi     = rng.uniform(0.0, 2.0 * np.pi)

    dx = np.sin(zen) * np.cos(azi)
    dy = np.sin(zen) * np.sin(azi)
    dz = cos_zen                       # >0 = upgoing

    # -- aim track through IceCube centre (0, 0, CENTRE_Z_KM)
    #    step back along direction until z = INJECT_Z_KM
    if dz < 1e-6:
        # nearly horizontal — skip, unlikely to produce many hits
        continue
    t_to_centre = (CENTRE_Z_KM - INJECT_Z_KM) / dz   # km along track
    x0 = -dx * t_to_centre
    y0 = -dy * t_to_centre
    pos_km = [x0, y0, INJECT_Z_KM]

    try:
        p_init     = make_particle(ene_GeV, pos_km, [dx, dy, dz])
        ene_losses = BLO.propagate(p_init, jlx(PROP_KM, jl.km))
        hits       = BLO.run_ppc(p_init, ene_losses,
                                 suppress_error=True, use_gpu=USE_GPU)
        uhits      = BLO.process_hits(hits)
    except Exception as exc:
        print(f"  [WARN] sim {n_simulated+1} error: {exc}")
        n_simulated += 1
        continue

    n_simulated += 1

    nhits_arr = np.array(uhits.nhits)
    total_hits = int(nhits_arr.sum())
    n_doms     = len(nhits_arr)

    if total_hits <= HIT_CUT:
        continue

    # -- accepted: extract hit data to numpy
    n_accepted += 1

    ev_energy_GeV.append(ene_GeV)
    ev_zenith_rad.append(zen)
    ev_azimuth_rad.append(azi)
    ev_n_hits.append(total_hits)
    ev_n_doms.append(n_doms)

    # positions and times — divide out units
    ev_dom_x.append(np.array(uhits.pos.x) / float(jlx(1.0, jl.m)))
    ev_dom_y.append(np.array(uhits.pos.y) / float(jlx(1.0, jl.m)))
    ev_dom_z.append(np.array(uhits.pos.z) / float(jlx(1.0, jl.m)))
    ev_dom_t.append(np.array(uhits.time)  / float(jlx(1.0, jl.ns)))
    ev_dom_nhit.append(nhits_arr)
    ev_dom_str.append(np.array(uhits.string_id))
    ev_dom_sen.append(np.array(uhits.sensor_id))

    print(f"  [{n_accepted:>3}/{N_TARGET}]  "
          f"tried={n_simulated:<4}  "
          f"E={ene_GeV/1e3:>6.2f} TeV  "
          f"zen={np.degrees(zen):>5.1f} deg  "
          f"hits={total_hits:>5}  doms={n_doms:>4}")

# ── Save ──────────────────────────────────────────────────────────────────────

print()
efficiency = 100.0 * n_accepted / n_simulated
print(f"[INFO] {n_accepted}/{n_simulated} accepted  (efficiency {efficiency:.1f}%)")
print(f"[INFO] Saving to {output_file}")

np.savez(
    output_file,
    # per-event scalars
    energy_GeV  = np.array(ev_energy_GeV),
    zenith_rad  = np.array(ev_zenith_rad),
    azimuth_rad = np.array(ev_azimuth_rad),
    n_hits      = np.array(ev_n_hits, dtype=int),
    n_doms      = np.array(ev_n_doms, dtype=int),
    # per-event ragged DOM arrays (stored as object arrays)
    dom_x       = np.array(ev_dom_x,    dtype=object),
    dom_y       = np.array(ev_dom_y,    dtype=object),
    dom_z       = np.array(ev_dom_z,    dtype=object),
    dom_t       = np.array(ev_dom_t,    dtype=object),
    dom_nhits   = np.array(ev_dom_nhit, dtype=object),
    dom_string  = np.array(ev_dom_str,  dtype=object),
    dom_sensor  = np.array(ev_dom_sen,  dtype=object),
)

print(f"[INFO] Done — {n_accepted} events saved to {output_file}")
print()
print("To reload:")
print(f"  d = np.load('{output_file}', allow_pickle=True)")
print(f"  d['energy_GeV']   # shape ({n_accepted},)")
print(f"  d['dom_x'][0]     # x positions of hit DOMs in event 0  [m]")
