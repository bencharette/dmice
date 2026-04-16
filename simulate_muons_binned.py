#!/usr/bin/env python3
"""
simulate_muons_binned.py — Binned energy muon simulation via BLO + GPU PPC.

Simulates 5000 downgoing muons through IceCube + DM-Ice geometry, split evenly
across 5 log-spaced energy bins from 100 GeV to 100 TeV (1000 events per bin).

Energy bins (log10 GeV):
    Bin 0:  100 –  398 GeV
    Bin 1:  398 GeV –  1.58 TeV
    Bin 2: 1.58 –  6.31 TeV
    Bin 3: 6.31 – 25.1 TeV
    Bin 4: 25.1 – 100 TeV

Direction:
    Downgoing: zenith 0–60° from vertical (cos_zen in [0.5, 1.0])
    Azimuth: uniform 0–360°

    NOTE: blo_python.py uses the same direction convention as BLO/PPC where
    zenith=0 is straight up and dz > 0 is upgoing. For downgoing (zenith 0–60°
    from vertical in the standard IceCube sense), we use dz < 0 and
    cos_zen in [-1.0, -0.5] in the BLO frame.

Tracks aimed through IceCube centre (0, 0, -1950 m), starting 2.5 km below.

Output: ~/dmice_work/output/muons_binned_5000ev.npz
    Arrays: energy_GeV, zenith_rad, azimuth_rad, n_hits, n_doms, bin_id,
            dom_x, dom_y, dom_z, dom_t, dom_nhits, dom_string, dom_sensor

Usage (on WARD with GPU PPC):
    BLO_PPC_EXE=~/.icevenv/BLO/resources/PPC_executables/PPC_CUDA/ppc \\
        python3 ~/dmice/simulate_muons_binned.py
"""

import os
import sys
import argparse
import numpy as np

sys.path.insert(0, os.path.expanduser("~/dmice"))
import blo_python as blo

# ── CLI args ──────────────────────────────────────────────────────────────────

parser = argparse.ArgumentParser()
parser.add_argument("--n-per-bin", type=int, default=1000,
                    help="Events per energy bin (default 1000)")
parser.add_argument("--n-bins", type=int, default=5,
                    help="Number of log-spaced energy bins (default 5)")
parser.add_argument("--e-min", type=float, default=1e2,
                    help="Minimum energy in GeV (default 100)")
parser.add_argument("--e-max", type=float, default=1e5,
                    help="Maximum energy in GeV (default 100000)")
parser.add_argument("--output", type=str, default=None,
                    help="Output npz path (default auto-named by n-per-bin and n-bins)")
parser.add_argument("--detector", type=str, default=None, choices=["det1", "det2", "det_center"],
                    help="Lock all events to one DM-Ice detector (default: alternate det1/det2)")
args = parser.parse_args()

# ── Simulation parameters ─────────────────────────────────────────────────────

N_BINS      = args.n_bins
N_PER_BIN   = args.n_per_bin
E_MIN_GEV   = args.e_min
E_MAX_GEV   = args.e_max
# Injection from ABOVE detector (downgoing muons enter from the top)
INJECT_Z_KM = -1.3    # km — above IC86 top (~-1.45 km)
PROP_KM     = 3.0     # km — propagates through full IC86 + DM-Ice depth

# DM-Ice detector positions in BLO coordinates [km]
# det_center: IceCube geometric center (0,0,0) IC coords → BLO z = (0 - 1948.07)/1000
DMICE = {
    "det1":       np.array([ 0.03125,  -0.07293, -2.45912]),
    "det2":       np.array([-0.33480,  -0.42450, -2.45933]),
    "det_center": np.array([ 0.0,       0.0,     -1.94807]),
}

# log10 bin edges
LOG_EDGES = np.linspace(np.log10(E_MIN_GEV), np.log10(E_MAX_GEV), N_BINS + 1)
BIN_EDGES_GEV = 10 ** LOG_EDGES

rng = np.random.default_rng(seed=42)

# ── Output ────────────────────────────────────────────────────────────────────

output_dir  = os.path.expanduser("~/dmice_work/output/")
_det_suffix  = f"_{args.detector}" if args.detector else ""
_default_out = os.path.join(output_dir, f"muons_binned_{N_BINS}bins_{N_PER_BIN}pbin{_det_suffix}.npz")
output_file = args.output if args.output else _default_out
os.makedirs(output_dir, exist_ok=True)

# ── NaI direct ionization timing model ───────────────────────────────────────
# Real DM-Ice detects muons via direct ionization/scintillation in NaI crystal,
# NOT via Cherenkov photons. This model injects a guaranteed hit regardless of
# PPC photon yield, using the calibrated timing parameters from real 2012 data.
C_M_NS_SIM   = 0.2998   # speed of light [m/ns]
MU_SCINT     = 280.0    # NaI mean scintillation delay [ns]
SIGMA_SCINT  =  81.0    # NaI scintillation jitter [ns]

# ── Storage ───────────────────────────────────────────────────────────────────

ev_energy_GeV  = []
ev_zenith_rad  = []
ev_azimuth_rad = []
ev_n_hits      = []
ev_n_doms      = []
ev_bin_id      = []
ev_target_det  = []   # 0=det1, 1=det2
ev_dom_x, ev_dom_y, ev_dom_z = [], [], []
ev_dom_t, ev_dom_nhit        = [], []
ev_dom_str, ev_dom_sen       = [], []
ev_dm_t_injected  = []   # analytically injected DM-Ice hit time [ns]
ev_dm_t_ppc       = []   # DM-Ice hit time from PPC (NaN if PPC missed)
ev_smt8_triggered = []   # bool: passed SMT8 HLC trigger (≥8 LC hits in 5 µs)

# ── Print bin table ───────────────────────────────────────────────────────────

print(f"[INFO] Simulating {N_BINS * N_PER_BIN} events: {N_PER_BIN}/bin × {N_BINS} bins")
print(f"[INFO] Energy bins (GeV):")
for b in range(N_BINS):
    lo, hi = BIN_EDGES_GEV[b], BIN_EDGES_GEV[b + 1]
    print(f"  Bin {b}: {lo:.1f} – {hi:.1f} GeV")
print(f"[INFO] Direction: downgoing, zenith 0–60° from vertical, azimuth 0–360°")
print(f"[INFO] PPC binary: {blo.PPC_EXE}")
print()

# ── Main loop ─────────────────────────────────────────────────────────────────

for bin_id in range(N_BINS):
    log_lo = LOG_EDGES[bin_id]
    log_hi = LOG_EDGES[bin_id + 1]

    print(f"── Bin {bin_id}: {10**log_lo:.1f} – {10**log_hi:.1f} GeV ──")

    for ev in range(N_PER_BIN):
        # energy: uniform in log within bin
        ene_GeV = 10 ** rng.uniform(log_lo, log_hi)

        # direction: downgoing, zenith 0–60° from vertical
        # dz < 0 → moving toward greater depth (downward)
        cos_zen_std = rng.uniform(0.5, 1.0)   # cos of angle from vertical (0–60°)
        dz  = -cos_zen_std                     # downgoing in BLO frame
        sin_zen = np.sqrt(1.0 - cos_zen_std**2)
        azi = rng.uniform(0.0, 2.0 * np.pi)
        dx  = sin_zen * np.cos(azi)
        dy  = sin_zen * np.sin(azi)

        # zenith in BLO frame for storage
        zen_blo = np.arccos(dz)   # in [pi/2, pi] for downgoing

        # aim track through DM-Ice detector
        if args.detector:
            target = args.detector
        else:
            target = "det1" if (bin_id * N_PER_BIN + ev) % 2 == 0 else "det2"
        target_id = 0 if target == "det1" else 1
        det_km = DMICE[target]

        # back-project from det position to injection height
        # t = (z_det - z_inject) / dz, dz < 0, z_det < z_inject → t > 0
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

        # ── Direct ionization DM-Ice hit (NaI scintillation model) ──────────
        # The muon always passes through the target DM-Ice crystal (by construction).
        # Compute the geometric crossing time and sample a scintillation delay.
        t_km_to_det = (det_km[2] - INJECT_Z_KM) / dz   # > 0: det is below injection
        t_cross_ns  = (t_km_to_det * 1000.0) / C_M_NS_SIM
        dm_t_inject = t_cross_ns + rng.normal(MU_SCINT, SIGMA_SCINT)

        try:
            losses = blo.propagate(p, dist_km=PROP_KM)
            hits   = blo.run_ppc(p, losses, suppress_error=True)
            doms   = blo.process_hits(hits)
        except Exception as exc:
            print(f"  [WARN] bin {bin_id} ev {ev}: {exc}")
            # still store the event with 0 hits so bin counts stay exact
            doms = {"x": np.array([]), "y": np.array([]), "z": np.array([]),
                    "t": np.array([]), "nhits": np.array([]),
                    "string_id": np.array([], dtype=int),
                    "sensor_id": np.array([], dtype=int)}

        # Extract PPC DM-Ice time (NaN if PPC didn't fire those strings)
        dm_str = 87 if target == "det1" else 88
        dm_mask = np.asarray(doms.get("string_id", []), dtype=int) == dm_str
        dm_t_ppc = float(np.asarray(doms["t"])[dm_mask].min()) \
                   if np.any(dm_mask) else float("nan")

        triggered, _ = blo.smt8_trigger(doms)

        n_doms = len(doms["x"])
        n_hits = int(doms["nhits"].sum()) if n_doms else 0

        ev_energy_GeV.append(ene_GeV)
        ev_zenith_rad.append(zen_blo)
        ev_azimuth_rad.append(azi)
        ev_n_hits.append(n_hits)
        ev_n_doms.append(n_doms)
        ev_bin_id.append(bin_id)
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

        total_done = bin_id * N_PER_BIN + ev + 1
        print(f"  [{total_done:>3}/{N_BINS*N_PER_BIN}]  "
              f"E={ene_GeV/1e3:>7.3f} TeV  "
              f"zen={np.degrees(zen_blo):>5.1f}°  "
              f"doms={n_doms:>4}  hits={n_hits:>5}  "
              f"smt8={'Y' if triggered else 'n'}")

# ── Save ──────────────────────────────────────────────────────────────────────

print()
print(f"[INFO] Saving to {output_file}")

np.savez(
    output_file,
    energy_GeV  = np.array(ev_energy_GeV),
    zenith_rad  = np.array(ev_zenith_rad),
    azimuth_rad = np.array(ev_azimuth_rad),
    n_hits      = np.array(ev_n_hits,  dtype=int),
    n_doms      = np.array(ev_n_doms,  dtype=int),
    bin_id      = np.array(ev_bin_id,    dtype=int),
    target_det  = np.array(ev_target_det, dtype=int),   # 0=det1, 1=det2
    bin_edges   = BIN_EDGES_GEV,
    dom_x       = np.array(ev_dom_x,   dtype=object),
    dom_y       = np.array(ev_dom_y,   dtype=object),
    dom_z       = np.array(ev_dom_z,   dtype=object),
    dom_t       = np.array(ev_dom_t,   dtype=object),
    dom_nhits   = np.array(ev_dom_nhit, dtype=object),
    dom_string       = np.array(ev_dom_str,       dtype=object),
    dom_sensor       = np.array(ev_dom_sen,       dtype=object),
    dm_t_injected_ns = np.array(ev_dm_t_injected,  dtype=float),
    dm_t_ppc_ns      = np.array(ev_dm_t_ppc,       dtype=float),
    smt8_triggered   = np.array(ev_smt8_triggered, dtype=bool),
)

print(f"[INFO] Done — {N_BINS * N_PER_BIN} events saved.")
print()
print("SMT8 trigger efficiency per energy bin:")
print(f"  {'Bin':>3}  {'Energy range':>22}  {'Triggered':>9}  {'Total':>5}  {'Efficiency':>10}")
smt8_arr = np.array(ev_smt8_triggered)
bin_arr  = np.array(ev_bin_id)
for b in range(N_BINS):
    lo, hi    = BIN_EDGES_GEV[b], BIN_EDGES_GEV[b + 1]
    bin_mask  = bin_arr == b
    n_total   = bin_mask.sum()
    n_trig    = smt8_arr[bin_mask].sum()
    eff       = n_trig / n_total if n_total else 0.0
    print(f"  {b:>3}  {lo:>8.1f} – {hi:>8.1f} GeV  {n_trig:>9}  {n_total:>5}  {eff:>9.1%}")
overall_eff = smt8_arr.mean()
print(f"  {'All':>3}  {'':>22}  {smt8_arr.sum():>9}  {len(smt8_arr):>5}  {overall_eff:>9.1%}")
print()
print("To reload:")
print(f"  d = np.load('{output_file}', allow_pickle=True)")
print(f"  d['energy_GeV']   # shape ({N_BINS * N_PER_BIN},)")
print(f"  d['bin_id']       # 0–{N_BINS-1}, 1000 events each")
print(f"  d['bin_edges']    # GeV bin edges, shape ({N_BINS+1},)")
