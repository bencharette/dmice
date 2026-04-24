#!/usr/bin/env python3
"""
inject_dmice_times.py

Patches an existing BLO simulation NPZ to add analytically computed DM-Ice
hit times based on direct ionization (not Cherenkov photon propagation).

For every event the muon track is geometrically guaranteed to pass through the
target DM-Ice detector (that's how simulate_muons_binned.py aims the tracks).
The crossing time is computed from the stored track parameters, then a NaI
scintillation delay is sampled from Gaussian(mu=+280ns, sigma=81ns).

The output NPZ is identical to the input but adds:
  dm_t_injected_ns  — per-event DM-Ice hit time [ns], always finite
  dm_t_ppc_ns       — per-event DM-Ice time from PPC (NaN if PPC missed it)

Usage:
  python3 ~/dmice/inject_dmice_times.py \
      --input  ~/dmice_work/output/muons_binned_5000ev_repacked.npz \
      --output ~/dmice_work/output/muons_binned_5000ev_repacked_injected.npz
"""

import os
import sys
import argparse
import numpy as np

# ── Constants ─────────────────────────────────────────────────────────────────
C_M_NS      =  0.2998          # speed of light [m/ns]
MU_SCINT    =  280.0           # NaI mean scintillation delay [ns]
SIGMA_SCINT =   81.0           # NaI scintillation jitter [ns]
INJECT_Z_KM =   -1.3          # injection z [km] — must match simulate_muons_binned.py

# DM-Ice positions in BLO coordinates [km]
DMICE_KM = {
    0: np.array([ 0.03125,  -0.07293, -2.45912]),   # det1
    1: np.array([-0.33480,  -0.42450, -2.45933]),   # det2
}

# DM-Ice OMKeys (string_id, sensor_id) — strings 87/88
DMICE_OMKEYS = {
    0: (87, 1),   # det1
    1: (88, 1),   # det2
}


def compute_crossing_time(zenith_rad, azimuth_rad, target_det):
    """
    Compute the time [ns] at which the muon crosses the target DM-Ice detector.
    Track starts at injection height INJECT_Z_KM with time_ns = 0.
    """
    # Reconstruct direction (BLO convention: dz < 0 for downgoing)
    cos_zen = np.cos(zenith_rad)   # zenith_rad is angle from vertical in BLO frame
    sin_zen = np.sin(zenith_rad)
    dz = -abs(cos_zen)             # always downgoing
    dx = sin_zen * np.cos(azimuth_rad)
    dy = sin_zen * np.sin(azimuth_rad)

    det_km = DMICE_KM[int(target_det)]

    # t_km: how far along the track direction until we reach detector z
    # z_injection + t_km * dz = det_km[2]  =>  t_km = (det_km[2] - INJECT_Z_KM) / dz
    t_km = (det_km[2] - INJECT_Z_KM) / dz   # dz < 0, det below injection → t_km > 0

    # Distance in metres → time in ns
    t_ns = (t_km * 1000.0) / C_M_NS
    return t_ns


def extract_ppc_dm_time(dom_string, dom_t, target_det):
    """
    Extract the earliest DM-Ice hit time from PPC DOM arrays.
    Returns NaN if PPC produced no hit on the DM-Ice strings.
    """
    dm_s, _ = DMICE_OMKEYS[int(target_det)]
    mask = np.asarray(dom_string, dtype=int) == dm_s
    if not np.any(mask):
        return float("nan")
    return float(np.asarray(dom_t)[mask].min())


def unpack_event(flat, offsets, i):
    """Unpack one event from a flat+offsets array pair."""
    lo = int(offsets[i])
    hi = int(offsets[i + 1]) if i + 1 < len(offsets) else len(flat)
    return flat[lo:hi]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input",  required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--seed",   type=int, default=42)
    args = parser.parse_args()

    print(f"Loading {args.input}")
    d = np.load(args.input, allow_pickle=True)
    n = len(d["energy_GeV"])
    print(f"  {n} events")

    rng = np.random.default_rng(args.seed)

    zenith    = d["zenith_rad"]
    azimuth   = d["azimuth_rad"]
    target    = d["target_det"]

    # Flat+offsets format
    str_flat    = d["dom_string_flat"].astype(int)
    str_offsets = d["dom_string_offsets"].astype(int)
    t_flat      = d["dom_t_flat"].astype(float)
    t_offsets   = d["dom_t_offsets"].astype(int)

    dm_t_injected = np.zeros(n, dtype=float)
    dm_t_ppc      = np.full(n, float("nan"))

    n_ppc_hit = 0
    for i in range(n):
        # Geometric crossing time
        t_cross = compute_crossing_time(zenith[i], azimuth[i], target[i])

        # Sample NaI scintillation delay
        scint_delay = rng.normal(MU_SCINT, SIGMA_SCINT)
        dm_t_injected[i] = t_cross + scint_delay

        # PPC time (may be NaN)
        dom_str_i = unpack_event(str_flat, str_offsets, i)
        dom_t_i   = unpack_event(t_flat,   t_offsets,   i)
        ppc_t = extract_ppc_dm_time(dom_str_i, dom_t_i, target[i])
        dm_t_ppc[i] = ppc_t
        if not np.isnan(ppc_t):
            n_ppc_hit += 1

        if (i + 1) % 500 == 0:
            print(f"  [{i+1}/{n}]  t_cross={t_cross:.1f}ns  "
                  f"injected={dm_t_injected[i]:.1f}ns  "
                  f"ppc={'NaN' if np.isnan(ppc_t) else f'{ppc_t:.1f}ns'}")

    print(f"\nPPC hit coverage: {n_ppc_hit}/{n} ({100*n_ppc_hit/n:.1f}%)")
    print(f"Injected coverage: {n}/{n} (100%)")

    # Save patched NPZ
    arrays = {k: d[k] for k in d.files}
    arrays["dm_t_injected_ns"] = dm_t_injected
    arrays["dm_t_ppc_ns"]      = dm_t_ppc

    print(f"\nSaving to {args.output}")
    np.savez(args.output, **arrays)
    print("Done.")


if __name__ == "__main__":
    main()
