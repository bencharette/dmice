#!/usr/bin/env python3
"""
simulate_muons_offset.py — Muon sim with configurable d_perp offset from DM-Ice.

Extends simulate_muons_binned.py to inject tracks at a perpendicular offset from
the DM-Ice detector. Used to build the d_perp-dependent timing model for the
DM-Ice direct likelihood.

Offset is applied as a random-direction perpendicular displacement from the
track-through-detector injection point, so d_perp ≈ offset for all events.

Output: ~/dmice_work/output/muons_offset_{OFFSET}m_{N}ev.npz

Usage (on WARD with GPU PPC):
    BLO_PPC_EXE=~/.icevenv/BLO/resources/PPC_executables/PPC_CUDA/ppc \\
        python3 ~/dmice/simulate_muons_offset.py --offset 50 --n 500

    # on-axis (equivalent to simulate_muons_binned.py):
    python3 ~/dmice/simulate_muons_offset.py --offset 0 --n 1000
"""

import os
import sys
import argparse
import numpy as np

sys.path.insert(0, os.path.expanduser("~/dmice"))
import blo_python as blo

# ── CLI ───────────────────────────────────────────────────────────────────────

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--offset", type=float, default=0.0,
                    help="Perpendicular offset from DM-Ice [m] (default: 0)")
parser.add_argument("--n", type=int, default=500,
                    help="Total events to simulate (default: 500)")
parser.add_argument("--seed", type=int, default=None,
                    help="RNG seed (default: random)")
parser.add_argument("--out", type=str, default=None,
                    help="Output npz path (default: auto)")
args = parser.parse_args()

OFFSET_M  = args.offset
N_TOTAL   = args.n
N_BINS    = 5
N_PER_BIN = N_TOTAL // N_BINS   # round down; any remainder dropped

# ── Simulation parameters ─────────────────────────────────────────────────────

E_MIN_GEV   = 1e2
E_MAX_GEV   = 1e5
INJECT_Z_KM = -1.3
PROP_KM     = 3.0

DMICE = {
    "det1": np.array([ 0.03125,  -0.07293, -2.45912]),   # km, BLO coords
    "det2": np.array([-0.33480,  -0.42450, -2.45933]),
}

LOG_EDGES     = np.linspace(np.log10(E_MIN_GEV), np.log10(E_MAX_GEV), N_BINS + 1)
BIN_EDGES_GEV = 10 ** LOG_EDGES

rng = np.random.default_rng(seed=args.seed)

# ── Output path ───────────────────────────────────────────────────────────────

output_dir = os.path.expanduser("~/dmice_work/output/")
os.makedirs(output_dir, exist_ok=True)
if args.out:
    output_file = args.out
else:
    tag = f"muons_offset_{int(OFFSET_M)}m_{N_PER_BIN * N_BINS}ev"
    output_file = os.path.join(output_dir, f"{tag}.npz")

# ── Storage ───────────────────────────────────────────────────────────────────

ev_energy_GeV  = []
ev_zenith_rad  = []
ev_azimuth_rad = []
ev_n_hits      = []
ev_n_doms      = []
ev_bin_id      = []
ev_target_det  = []
ev_offset_m    = []
ev_dom_x, ev_dom_y, ev_dom_z = [], [], []
ev_dom_t, ev_dom_nhit        = [], []
ev_dom_str, ev_dom_sen       = [], []

# ── Helpers ───────────────────────────────────────────────────────────────────

def perp_offset_km(d_hat, offset_m, phi):
    """
    Return a perpendicular offset vector in km given track direction d_hat,
    offset magnitude in metres, and azimuthal angle phi around the track axis.
    """
    if offset_m == 0.0:
        return np.zeros(3)

    # Build two unit vectors perpendicular to d_hat
    ref = np.array([0.0, 0.0, 1.0])
    if abs(np.dot(d_hat, ref)) > 0.99:
        ref = np.array([1.0, 0.0, 0.0])
    v1 = np.cross(d_hat, ref)
    v1 /= np.linalg.norm(v1)
    v2 = np.cross(d_hat, v1)
    v2 /= np.linalg.norm(v2)

    offset_km = (offset_m / 1e3) * (np.cos(phi) * v1 + np.sin(phi) * v2)
    return offset_km

# ── Print header ──────────────────────────────────────────────────────────────

N_SIM = N_BINS * N_PER_BIN
print(f"[INFO] simulate_muons_offset.py")
print(f"[INFO]   Offset:  {OFFSET_M:.0f} m perpendicular to track")
print(f"[INFO]   Events:  {N_SIM}  ({N_PER_BIN}/bin × {N_BINS} bins)")
print(f"[INFO]   Output:  {output_file}")
print(f"[INFO]   PPC:     {blo.PPC_EXE}")
print()
for b in range(N_BINS):
    lo, hi = BIN_EDGES_GEV[b], BIN_EDGES_GEV[b + 1]
    print(f"  Bin {b}: {lo:.1f} – {hi:.1f} GeV")
print()

# ── Main loop ─────────────────────────────────────────────────────────────────

for bin_id in range(N_BINS):
    log_lo = LOG_EDGES[bin_id]
    log_hi = LOG_EDGES[bin_id + 1]
    print(f"── Bin {bin_id}: {10**log_lo:.1f} – {10**log_hi:.1f} GeV ──")

    for ev in range(N_PER_BIN):
        ene_GeV = 10 ** rng.uniform(log_lo, log_hi)

        # direction: downgoing, zenith 0–60° from vertical
        cos_zen = rng.uniform(0.5, 1.0)
        dz      = -cos_zen
        sin_zen = np.sqrt(1.0 - cos_zen**2)
        azi     = rng.uniform(0.0, 2.0 * np.pi)
        dx      = sin_zen * np.cos(azi)
        dy      = sin_zen * np.sin(azi)
        d_hat   = np.array([dx, dy, dz])
        zen_blo = np.arccos(dz)

        # target detector (alternating)
        target    = "det1" if (bin_id * N_PER_BIN + ev) % 2 == 0 else "det2"
        target_id = 0 if target == "det1" else 1
        det_km    = DMICE[target]

        # back-project from det to injection height (on-axis)
        t_km   = (det_km[2] - INJECT_Z_KM) / dz
        x0_km  = det_km[0] - dx * t_km
        y0_km  = det_km[1] - dy * t_km

        # add perpendicular offset (random azimuth around track axis)
        phi_off    = rng.uniform(0.0, 2.0 * np.pi)
        off_vec_km = perp_offset_km(d_hat, OFFSET_M, phi_off)

        x0_m = (x0_km + off_vec_km[0]) * 1e3
        y0_m = (y0_km + off_vec_km[1]) * 1e3
        z0_m = INJECT_Z_KM * 1e3

        p = blo.ParticleState(
            energy_GeV = ene_GeV,
            pos_m      = [x0_m, y0_m, z0_m],
            dir        = [dx, dy, dz],
            pid        = 13,
            time_ns    = 0.0,
        )

        try:
            losses = blo.propagate(p, dist_km=PROP_KM)
            hits   = blo.run_ppc(p, losses, suppress_error=True)
            doms   = blo.process_hits(hits)
        except Exception as exc:
            print(f"  [WARN] bin {bin_id} ev {ev}: {exc}")
            doms = {"x": np.array([]), "y": np.array([]), "z": np.array([]),
                    "t": np.array([]), "nhits": np.array([]),
                    "string_id": np.array([], dtype=int),
                    "sensor_id": np.array([], dtype=int)}

        n_doms = len(doms["x"])
        n_hits = int(doms["nhits"].sum()) if n_doms else 0

        ev_energy_GeV.append(ene_GeV)
        ev_zenith_rad.append(zen_blo)
        ev_azimuth_rad.append(azi)
        ev_n_hits.append(n_hits)
        ev_n_doms.append(n_doms)
        ev_bin_id.append(bin_id)
        ev_target_det.append(target_id)
        ev_offset_m.append(OFFSET_M)
        ev_dom_x.append(doms["x"])
        ev_dom_y.append(doms["y"])
        ev_dom_z.append(doms["z"])
        ev_dom_t.append(doms["t"])
        ev_dom_nhit.append(doms["nhits"])
        ev_dom_str.append(doms["string_id"])
        ev_dom_sen.append(doms["sensor_id"])

        total_done = bin_id * N_PER_BIN + ev + 1
        print(f"  [{total_done:>4}/{N_SIM}]  "
              f"E={ene_GeV/1e3:>7.3f} TeV  "
              f"zen={np.degrees(zen_blo):>5.1f}°  "
              f"doms={n_doms:>4}  hits={n_hits:>6}")

# ── Save ──────────────────────────────────────────────────────────────────────

print()
print(f"[INFO] Saving {N_SIM} events → {output_file}")

np.savez(
    output_file,
    energy_GeV  = np.array(ev_energy_GeV),
    zenith_rad  = np.array(ev_zenith_rad),
    azimuth_rad = np.array(ev_azimuth_rad),
    n_hits      = np.array(ev_n_hits,     dtype=int),
    n_doms      = np.array(ev_n_doms,     dtype=int),
    bin_id      = np.array(ev_bin_id,     dtype=int),
    target_det  = np.array(ev_target_det, dtype=int),
    offset_m    = np.array(ev_offset_m),
    bin_edges   = BIN_EDGES_GEV,
    dom_x       = np.array(ev_dom_x,    dtype=object),
    dom_y       = np.array(ev_dom_y,    dtype=object),
    dom_z       = np.array(ev_dom_z,    dtype=object),
    dom_t       = np.array(ev_dom_t,    dtype=object),
    dom_nhits   = np.array(ev_dom_nhit, dtype=object),
    dom_string  = np.array(ev_dom_str,  dtype=object),
    dom_sensor  = np.array(ev_dom_sen,  dtype=object),
)

print(f"[INFO] Done — {N_SIM} events saved.")
