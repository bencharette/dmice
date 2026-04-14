#!/usr/bin/env python3
"""
Run blo_python.py pipeline on the 2 events in blo_muons_200hits.npz.
Reconstructs injection position from stored zenith/azimuth using the
same geometry as batch_dm_ice_sim.py, then runs propagate + run_ppc.

Output: blo_muons_200hits_rerun.npz  (same structure as input NPZ)
"""

import os, sys
import numpy as np

sys.path.insert(0, os.path.expanduser("~/dmice"))
import blo_python as blo

# ── Match batch_dm_ice_sim.py injection geometry ──────────────────────────────
INJECT_Z_KM = -2.5
CENTRE_Z_KM = -1.950
PROP_KM     = 3.0

# ── Load input ────────────────────────────────────────────────────────────────
NPZ_IN  = os.path.expanduser("~/dmice_work/output/blo_muons_200hits.npz")
NPZ_OUT = os.path.expanduser("~/dmice_work/output/blo_muons_200hits_rerun.npz")

d = np.load(NPZ_IN, allow_pickle=True)
n_events = len(d["energy_GeV"])
print(f"Loaded {n_events} events from {NPZ_IN}")

# ── Output storage (same structure as batch_dm_ice_sim.py) ────────────────────
ev_energy_GeV  = []
ev_zenith_rad  = []
ev_azimuth_rad = []
ev_n_hits      = []
ev_n_doms      = []
ev_dom_x, ev_dom_y, ev_dom_z = [], [], []
ev_dom_t, ev_dom_nhit = [], []
ev_dom_str, ev_dom_sen = [], []

# ── Process each event ────────────────────────────────────────────────────────
for i in range(n_events):
    ene_GeV = float(d["energy_GeV"][i])
    zen     = float(d["zenith_rad"][i])
    azi     = float(d["azimuth_rad"][i])

    dx = np.sin(zen) * np.cos(azi)
    dy = np.sin(zen) * np.sin(azi)
    dz = np.cos(zen)

    # reconstruct injection position (same as batch_dm_ice_sim.py)
    t_to_centre = (CENTRE_Z_KM - INJECT_Z_KM) / dz
    x0_m = -dx * t_to_centre * 1e3
    y0_m = -dy * t_to_centre * 1e3
    z0_m = INJECT_Z_KM * 1e3

    print(f"\n[Event {i}]  E={ene_GeV/1e3:.1f} TeV  "
          f"zen={np.degrees(zen):.1f}  azi={np.degrees(azi):.1f}")
    print(f"  injection: ({x0_m:.0f}, {y0_m:.0f}, {z0_m:.0f}) m")

    p = blo.ParticleState(
        energy_GeV = ene_GeV,
        pos_m      = [x0_m, y0_m, z0_m],
        dir        = [dx, dy, dz],
        pid        = 13,
        time_ns    = 0.0,
    )

    print("  propagating...")
    losses = blo.propagate(p, dist_km=PROP_KM)
    print(f"  {len(losses)} loss segments")

    print("  running PPC...")
    hits = blo.run_ppc(p, losses, suppress_error=True)
    print(f"  {len(hits)} photon hits")

    doms = blo.process_hits(hits)
    n_doms = len(doms["x"])
    n_hits = int(doms["nhits"].sum())
    print(f"  {n_doms} DOMs hit, {n_hits} total hits")

    ev_energy_GeV.append(ene_GeV)
    ev_zenith_rad.append(zen)
    ev_azimuth_rad.append(azi)
    ev_n_hits.append(n_hits)
    ev_n_doms.append(n_doms)
    ev_dom_x.append(doms["x"])
    ev_dom_y.append(doms["y"])
    ev_dom_z.append(doms["z"])
    ev_dom_t.append(doms["t"])
    ev_dom_nhit.append(doms["nhits"])
    ev_dom_str.append(doms["string_id"])
    ev_dom_sen.append(doms["sensor_id"])

# ── Save ──────────────────────────────────────────────────────────────────────
np.savez(
    NPZ_OUT,
    energy_GeV  = np.array(ev_energy_GeV),
    zenith_rad  = np.array(ev_zenith_rad),
    azimuth_rad = np.array(ev_azimuth_rad),
    n_hits      = np.array(ev_n_hits, dtype=int),
    n_doms      = np.array(ev_n_doms, dtype=int),
    dom_x       = np.array(ev_dom_x,    dtype=object),
    dom_y       = np.array(ev_dom_y,    dtype=object),
    dom_z       = np.array(ev_dom_z,    dtype=object),
    dom_t       = np.array(ev_dom_t,    dtype=object),
    dom_nhits   = np.array(ev_dom_nhit, dtype=object),
    dom_string  = np.array(ev_dom_str,  dtype=object),
    dom_sensor  = np.array(ev_dom_sen,  dtype=object),
)
print(f"\nSaved to {NPZ_OUT}")
