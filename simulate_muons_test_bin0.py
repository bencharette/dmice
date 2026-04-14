#!/usr/bin/env python3
"""
simulate_muons_test_bin0.py — 200-event test sim from lowest energy bin.

Simulates 200 downgoing muons in the 100–398 GeV range (Bin 0) through
IceCube + DM-Ice geometry.  Intended for quick angular resolution testing
with SPE/MPE fits before committing to a full 5000-event run.

Output: ~/dmice_work/output/muons_test_bin0_200ev.npz

Usage (on WARD with GPU PPC):
    BLO_PPC_EXE=~/.icevenv/BLO/resources/PPC_executables/PPC_CUDA/ppc \\
        python3 ~/dmice/simulate_muons_test_bin0.py

Usage (on Cobalt with CPU PPC):
    python3 ~/dmice/simulate_muons_test_bin0.py
"""

import os
import sys
import numpy as np

sys.path.insert(0, os.path.expanduser("~/dmice"))
import blo_python as blo

# ── Parameters ────────────────────────────────────────────────────────────────

N_EVENTS    = 200
E_MIN_GEV   = 1e2      # 100 GeV  (Bin 0 lower edge)
E_MAX_GEV   = 10 ** (np.log10(1e2) + (np.log10(1e5) - np.log10(1e2)) / 5)  # ~398 GeV
INJECT_Z_KM = -1.3
PROP_KM     = 3.0

DMICE = {
    "det1": np.array([ 0.03125,  -0.07293, -2.45912]),
    "det2": np.array([-0.33480,  -0.42450, -2.45933]),
}

C_M_NS_SIM  = 0.2998
MU_SCINT    = 280.0
SIGMA_SCINT =  81.0

output_dir  = os.path.expanduser("~/dmice_work/output/")
output_file = os.path.join(output_dir, "muons_test_bin0_200ev.npz")
os.makedirs(output_dir, exist_ok=True)

rng = np.random.default_rng(seed=99)

# ── Storage ───────────────────────────────────────────────────────────────────

ev_energy_GeV     = []
ev_zenith_rad     = []
ev_azimuth_rad    = []
ev_n_hits         = []
ev_n_doms         = []
ev_target_det     = []
ev_dom_x, ev_dom_y, ev_dom_z = [], [], []
ev_dom_t, ev_dom_nhit        = [], []
ev_dom_str, ev_dom_sen       = [], []
ev_dm_t_injected  = []
ev_dm_t_ppc       = []
ev_smt8_triggered = []

print(f"[INFO] Simulating {N_EVENTS} events, Bin 0: {E_MIN_GEV:.1f}–{E_MAX_GEV:.1f} GeV")
print(f"[INFO] PPC binary: {blo.PPC_EXE}")
print()

# ── Main loop ─────────────────────────────────────────────────────────────────

for ev in range(N_EVENTS):
    ene_GeV = 10 ** rng.uniform(np.log10(E_MIN_GEV), np.log10(E_MAX_GEV))

    cos_zen_std = rng.uniform(0.5, 1.0)
    dz          = -cos_zen_std
    sin_zen     = np.sqrt(1.0 - cos_zen_std**2)
    azi         = rng.uniform(0.0, 2.0 * np.pi)
    dx          = sin_zen * np.cos(azi)
    dy          = sin_zen * np.sin(azi)
    zen_blo     = np.arccos(dz)

    target    = "det1" if ev % 2 == 0 else "det2"
    target_id = 0 if target == "det1" else 1
    det_km    = DMICE[target]

    t_km  = (det_km[2] - INJECT_Z_KM) / dz
    x0_m  = (det_km[0] - dx * t_km) * 1e3
    y0_m  = (det_km[1] - dy * t_km) * 1e3
    z0_m  = INJECT_Z_KM * 1e3

    p = blo.ParticleState(
        energy_GeV = ene_GeV,
        pos_m      = [x0_m, y0_m, z0_m],
        dir        = [dx, dy, dz],
        pid        = 13,
        time_ns    = 0.0,
    )

    t_cross_ns  = (t_km * 1000.0) / C_M_NS_SIM
    dm_t_inject = t_cross_ns + rng.normal(MU_SCINT, SIGMA_SCINT)

    try:
        losses = blo.propagate(p, dist_km=PROP_KM)
        hits   = blo.run_ppc(p, losses, suppress_error=True)
        doms   = blo.process_hits(hits)
    except Exception as exc:
        print(f"  [WARN] ev {ev}: {exc}")
        doms = {"x": np.array([]), "y": np.array([]), "z": np.array([]),
                "t": np.array([]), "nhits": np.array([]),
                "string_id": np.array([], dtype=int),
                "sensor_id": np.array([], dtype=int)}

    triggered, _ = blo.smt8_trigger(doms)

    dm_str  = 87 if target == "det1" else 88
    dm_mask = np.asarray(doms.get("string_id", []), dtype=int) == dm_str
    dm_t_ppc = float(np.asarray(doms["t"])[dm_mask].min()) \
               if np.any(dm_mask) else float("nan")

    n_doms = len(doms["x"])
    n_hits = int(doms["nhits"].sum()) if n_doms else 0

    ev_energy_GeV.append(ene_GeV)
    ev_zenith_rad.append(zen_blo)
    ev_azimuth_rad.append(azi)
    ev_n_hits.append(n_hits)
    ev_n_doms.append(n_doms)
    ev_target_det.append(target_id)
    ev_dom_x.append(doms["x"])
    ev_dom_y.append(doms["y"])
    ev_dom_z.append(doms["z"])
    ev_dom_t.append(doms["t"])
    ev_dom_nhit.append(doms["nhits"])
    ev_dom_str.append(doms["string_id"])
    ev_dom_sen.append(doms["sensor_id"])
    ev_dm_t_injected.append(dm_t_inject)
    ev_dm_t_ppc.append(dm_t_ppc)
    ev_smt8_triggered.append(triggered)

    print(f"  [{ev+1:>3}/{N_EVENTS}]  "
          f"E={ene_GeV:>7.1f} GeV  "
          f"zen={np.degrees(zen_blo):>5.1f}°  "
          f"doms={n_doms:>4}  hits={n_hits:>5}  "
          f"smt8={'Y' if triggered else 'n'}")

# ── Save ──────────────────────────────────────────────────────────────────────

print()
print(f"[INFO] Saving to {output_file}")

np.savez(
    output_file,
    energy_GeV       = np.array(ev_energy_GeV),
    zenith_rad        = np.array(ev_zenith_rad),
    azimuth_rad       = np.array(ev_azimuth_rad),
    n_hits            = np.array(ev_n_hits,         dtype=int),
    n_doms            = np.array(ev_n_doms,         dtype=int),
    bin_id            = np.zeros(N_EVENTS,           dtype=int),  # all Bin 0
    target_det        = np.array(ev_target_det,      dtype=int),
    bin_edges         = np.array([E_MIN_GEV, E_MAX_GEV]),
    dom_x             = np.array(ev_dom_x,           dtype=object),
    dom_y             = np.array(ev_dom_y,           dtype=object),
    dom_z             = np.array(ev_dom_z,           dtype=object),
    dom_t             = np.array(ev_dom_t,           dtype=object),
    dom_nhits         = np.array(ev_dom_nhit,        dtype=object),
    dom_string        = np.array(ev_dom_str,         dtype=object),
    dom_sensor        = np.array(ev_dom_sen,         dtype=object),
    dm_t_injected_ns  = np.array(ev_dm_t_injected,   dtype=float),
    dm_t_ppc_ns       = np.array(ev_dm_t_ppc,        dtype=float),
    smt8_triggered    = np.array(ev_smt8_triggered,  dtype=bool),
)

smt8_arr = np.array(ev_smt8_triggered)
print(f"[INFO] Done — {N_EVENTS} events saved.")
print(f"[INFO] SMT8 trigger efficiency: {smt8_arr.sum()}/{N_EVENTS} = {smt8_arr.mean():.1%}")
