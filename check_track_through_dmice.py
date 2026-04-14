#!/usr/bin/env python3
"""
Check whether simulated BLO muon tracks actually pass through DM-Ice det1/det2.

For each event, reconstructs the track from stored zenith/azimuth and the
known injection geometry (back-projected 2 km from DM-Ice), then computes
the closest-approach distance to each detector's position.
"""

import numpy as np

# ── DM-Ice positions (IceCube coords, meters) ─────────────────────────────────
DMICE_POS_M = {
    "det1": np.array([ 31.25,   -72.93,  -511.05]),
    "det2": np.array([-334.80, -424.50,  -511.26]),
}

# ── Load data ─────────────────────────────────────────────────────────────────
import sys, os
npz_path = os.path.expanduser(
    "~/dmice/WARD_dm_sim/blo_dmice_targeted_det1det2_both_1000events_repacked.npz"
)
if len(sys.argv) > 1:
    npz_path = sys.argv[1]

print(f"Loading {npz_path}")
d = np.load(npz_path, allow_pickle=True)

energy_GeV  = d["energy_GeV"]
zenith_rad  = d["zenith_rad"]
azimuth_rad = d["azimuth_rad"]
det_id      = d["det_id"]     # "det1" or "det2"
dir_type    = d["dir_type"]   # "up" or "down"
N = len(energy_GeV)
print(f"  {N} events")
print(f"  det_id counts: det1={np.sum(det_id=='det1')}, det2={np.sum(det_id=='det2')}")
print(f"  dir_type counts: down={np.sum(dir_type=='down')}, up={np.sum(dir_type=='up')}")
print()

# ── Reconstruct track and compute closest approach to each detector ───────────
BACKPROJECT_KM = 2.0   # same as simulation

# Results arrays
ca_targeted = np.zeros(N)   # closest approach to the targeted detector [m]
ca_other    = np.zeros(N)   # closest approach to the other detector [m]

for i in range(N):
    zen = zenith_rad[i]
    azi = azimuth_rad[i]
    det = det_id[i]

    # Unit direction vector
    dx = np.sin(zen) * np.cos(azi)
    dy = np.sin(zen) * np.sin(azi)
    dz = np.cos(zen)
    d_hat = np.array([dx, dy, dz])

    # Injection point: DM-Ice pos − 2 km * d_hat  (in meters)
    dmice = DMICE_POS_M[det]
    pos_m = dmice - BACKPROJECT_KM * 1000.0 * d_hat  # metres

    # Closest approach of line (pos_m + t*d_hat) to a point P:
    #   t* = (P - pos_m) · d_hat
    #   closest_approach = ||(P - pos_m) - t* * d_hat||
    def closest_approach_m(point):
        v = point - pos_m
        t_star = np.dot(v, d_hat)
        ca_vec = v - t_star * d_hat
        return np.linalg.norm(ca_vec)

    other_det = "det2" if det == "det1" else "det1"
    ca_targeted[i] = closest_approach_m(DMICE_POS_M[det])
    ca_other[i]    = closest_approach_m(DMICE_POS_M[other_det])

# ── Report ────────────────────────────────────────────────────────────────────
print("=== Closest-approach distance to TARGETED detector ===")
print(f"  mean:   {ca_targeted.mean():.4f} m")
print(f"  median: {np.median(ca_targeted):.4f} m")
print(f"  max:    {ca_targeted.max():.4f} m")
print(f"  frac < 1 m:  {np.mean(ca_targeted < 1.0):.3f}")
print(f"  frac < 10 m: {np.mean(ca_targeted < 10.0):.3f}")
print()
print("=== Closest-approach distance to OTHER detector ===")
print(f"  mean:   {ca_other.mean():.1f} m")
print(f"  median: {np.median(ca_other):.1f} m")
print(f"  min:    {ca_other.min():.1f} m")
print(f"  frac < 100 m: {np.mean(ca_other < 100.0):.3f}")
print()

# Per-detector breakdown
for det in ["det1", "det2"]:
    mask = det_id == det
    if mask.sum() == 0:
        continue
    print(f"--- {det} events (n={mask.sum()}) ---")
    print(f"  CA to {det}:  mean={ca_targeted[mask].mean():.4f} m  max={ca_targeted[mask].max():.4f} m")
    other = "det2" if det == "det1" else "det1"
    print(f"  CA to {other}: mean={ca_other[mask].mean():.1f} m  min={ca_other[mask].min():.1f} m")
print()

# Check: are any events accidentally also passing through the other detector?
DMICE_RADIUS_M = 0.5   # rough NaI crystal radius
accidental = ca_other < DMICE_RADIUS_M
print(f"Events accidentally passing through the OTHER detector (CA < {DMICE_RADIUS_M} m): {accidental.sum()}")
