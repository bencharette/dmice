#!/usr/bin/env python3
"""
DMice Targeted Muon Simulation
================================
Simulates downgoing muons guaranteed to pass through one of the DM-Ice NaI
detectors.  Physics:
  - DM-Ice is a NaI scintillator that detects DIRECT muon traversal only.
  - The muon physically crosses the crystal, giving a precise timing anchor
    (no Cherenkov propagation uncertainty) — this is the "pivot point".
  - Each event is aimed at det1 or det2 with zero miss distance.
  - The muon starts 1000 m back along the track from the DM-Ice position so
    it traverses a substantial fraction of IC86 before reaching DM-Ice.

Energy range: 1 TeV – 1 PeV (gamma=2)
  Based on Prometheus paper (Lazar et al. 2023, Fig. 3 example: 100 GeV–1 PeV
  gamma=1).  We raise the minimum to 1 TeV so muons produce enough Cherenkov
  light throughout IC86 for a reliable LineFit, and use gamma=2 to reflect the
  realistic atmospheric muon spectrum shape at depth.

Direction convention (Prometheus momentum convention):
  zenith=0   → straight up   (upgoing)
  zenith=90  → horizontal
  zenith=180 → straight down (downgoing)
  Downgoing atmospheric muons have zenith 100–170°.

Run on NPX (or Cobalt) inside the Prometheus environment:
    python simulate_muons.py [--run RUN_NUMBER] [--nevents N] [--det {1,2,both}]

Usage for Condor:
    python simulate_muons.py --run $(Process) --nevents 500
"""

import sys
import os
import argparse
import numpy as np
import h5py

# ── Paths ──────────────────────────────────────────────────────────────────────
prometheus_path = os.path.expanduser("~/prometheus")
if prometheus_path not in sys.path:
    sys.path.insert(0, prometheus_path)

import prometheus
from prometheus import config

prometheus_base = '/'.join(prometheus.__path__[0].split('/')[:-1])
resource_dir    = f"{prometheus_base}/resources/"

# ── DM-Ice positions in Prometheus depth coordinates (metres) ─────────────────
# IceCube coord z → depth z = z_icecube − 1948.07
DMICE_DEPTH = {
    "det1": np.array([ 31.25,  -72.93, -2459.12]),   # string 87, dom 1
    "det2": np.array([-334.80, -424.50, -2459.33]),   # string 88, dom 1
}

# PDG / Prometheus particle codes
PDG_NUMU      =  14
PDG_MUON      =  13
PDG_HADRONS   = -2000001006   # Prometheus hadronic shower code

# ── Physics constants ──────────────────────────────────────────────────────────
C_M_NS = 0.2998   # speed of light m/ns


def power_law_energies(rng, n, e_min, e_max, gamma):
    """Sample n energies from E^{-gamma} between e_min and e_max (GeV)."""
    alpha = 1.0 - gamma
    u = rng.uniform(0.0, 1.0, n)
    return (u * (e_max**alpha - e_min**alpha) + e_min**alpha) ** (1.0 / alpha)


def make_targeted_injection(out_h5, n_events, det_target, rng):
    """
    Write a LeptonInjector-format HDF5 injection file.

    Each muon:
      - Is aimed at det_target (one of 'det1', 'det2')
      - Has start position = dm_pos - 1000 m * direction  (back-projected)
      - Has random downgoing direction: zenith uniformly in [100°, 170°],
        azimuth uniformly in [0°, 360°]
      - Has energy drawn from E^{-2} between 1 TeV and 1 PeV

    The injection file follows the VolumeInjector0 schema that Prometheus
    injection_from_LI_output() expects.
    """
    dm_pos = DMICE_DEPTH[det_target]

    # Energy spectrum: 1 TeV – 3 TeV, gamma=2
    # Upper limit kept at 3 TeV so PPC (CPU-only) finishes in reasonable time
    # (~5-15 min per event on Condor CPU nodes).  The power-law with gamma=2
    # means ~75% of events fall below 1.7 TeV, keeping the batch tractable.
    # 1 TeV minimum ensures the muon survives the 1500 m track to DM-Ice.
    E_min, E_max, gamma = 1e3, 3e3, 2.0   # GeV
    e_muon = power_law_energies(rng, n_events, E_min, E_max, gamma)

    # Random downgoing directions (Prometheus momentum convention)
    # zenith in [130°, 170°] → muon traverses most of IC86 before DM-Ice.
    # Restricted to |dz| > 0.5 so the back-projected start point is above
    # IC86 without requiring an unreasonably long (> ~3 km) track.
    zen_min, zen_max = np.radians(130.0), np.radians(170.0)
    cos_min, cos_max = np.cos(zen_max), np.cos(zen_min)   # note reversal
    cos_zen = rng.uniform(cos_min, cos_max, n_events)
    zen     = np.arccos(cos_zen)
    azi     = rng.uniform(0.0, 2.0 * np.pi, n_events)

    # 3-D direction unit vectors (Prometheus momentum convention)
    dx = np.sin(zen) * np.cos(azi)
    dy = np.sin(zen) * np.sin(azi)
    dz = np.cos(zen)   # < 0 for downgoing (zenith > 90°)
    directions = np.column_stack([dx, dy, dz])   # (N, 3)

    # Start positions: back-project 1500 m along the track from DM-Ice.
    # With zenith 130-170° and 1500 m back-projection, start positions
    # land at IceCube z ≈ +450 to +1000 m — at or above the top of IC86
    # — so muons traverse the full detector volume before reaching DM-Ice.
    TRACK_BACK_M = 1500.0
    starts = dm_pos[None, :] - TRACK_BACK_M * directions   # (N, 3)

    # Neutrino energies (bjorken_y ≈ 0.5 → E_nu = 2 * E_muon)
    bjorken_x = rng.uniform(0.1, 0.5, n_events)
    bjorken_y = rng.uniform(0.3, 0.7, n_events)
    e_nu      = e_muon / (1.0 - bjorken_y)
    e_had     = bjorken_y * e_nu

    # Structured dtype matching existing LI output
    particle_dtype = np.dtype([
        ('initial',      'u1'),
        ('ParticleType', '<i4'),
        ('Position',     '<f8', (3,)),
        ('Direction',    '<f8', (2,)),
        ('Energy',       '<f8'),
    ], align=False)

    properties_dtype = np.dtype([
        ('totalEnergy',      '<f8'),
        ('zenith',           '<f8'),
        ('azimuth',          '<f8'),
        ('finalStateX',      '<f8'),
        ('finalStateY',      '<f8'),
        ('finalType1',       '<i4'),
        ('finalType2',       '<i4'),
        ('initialType',      '<i4'),
        ('x',                '<f8'),
        ('y',                '<f8'),
        ('z',                '<f8'),
        ('totalColumnDepth', '<f8'),
    ], align=False)

    initial_arr = np.zeros(n_events, dtype=particle_dtype)
    final1_arr  = np.zeros(n_events, dtype=particle_dtype)
    final2_arr  = np.zeros(n_events, dtype=particle_dtype)
    props_arr   = np.zeros(n_events, dtype=properties_dtype)

    for i in range(n_events):
        d2 = np.array([zen[i], azi[i]])   # [zenith, azimuth]
        pos = starts[i]

        # Initial neutrino
        initial_arr[i]['initial']      = 1
        initial_arr[i]['ParticleType'] = PDG_NUMU
        initial_arr[i]['Position']     = pos
        initial_arr[i]['Direction']    = d2
        initial_arr[i]['Energy']       = e_nu[i]

        # Final state 1: muon
        final1_arr[i]['initial']      = 0
        final1_arr[i]['ParticleType'] = PDG_MUON
        final1_arr[i]['Position']     = pos
        final1_arr[i]['Direction']    = d2
        final1_arr[i]['Energy']       = e_muon[i]

        # Final state 2: hadronic shower
        final2_arr[i]['initial']      = 0
        final2_arr[i]['ParticleType'] = PDG_HADRONS
        final2_arr[i]['Position']     = pos
        final2_arr[i]['Direction']    = d2
        final2_arr[i]['Energy']       = e_had[i]

        # Properties summary
        props_arr[i]['totalEnergy']      = e_nu[i]
        props_arr[i]['zenith']           = zen[i]
        props_arr[i]['azimuth']          = azi[i]
        props_arr[i]['finalStateX']      = bjorken_x[i]
        props_arr[i]['finalStateY']      = bjorken_y[i]
        props_arr[i]['finalType1']       = PDG_MUON
        props_arr[i]['finalType2']       = PDG_HADRONS
        props_arr[i]['initialType']      = PDG_NUMU
        props_arr[i]['x']                = pos[0]
        props_arr[i]['y']                = pos[1]
        props_arr[i]['z']                = pos[2]
        props_arr[i]['totalColumnDepth'] = 100000.0   # placeholder MWE

    with h5py.File(out_h5, 'w') as f:
        grp = f.create_group('VolumeInjector0')
        grp.create_dataset('initial',    data=initial_arr)
        grp.create_dataset('final_1',    data=final1_arr)
        grp.create_dataset('final_2',    data=final2_arr)
        grp.create_dataset('properties', data=props_arr)

    print(f"[INFO] Injection file written: {out_h5} ({n_events} events)")
    print(f"[INFO]   Target: {det_target} at depth {dm_pos}")
    print(f"[INFO]   Energy: {E_min:.0e}–{E_max:.0e} GeV (gamma={gamma})")
    print(f"[INFO]   Zenith: {np.degrees(zen_min):.0f}°–{np.degrees(zen_max):.0f}° (downgoing)")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run',     type=int, default=2000,
                        help='Run number (also used as RNG seed, default 2000)')
    parser.add_argument('--nevents', type=int, default=500,
                        help='Number of events to simulate (default 500)')
    parser.add_argument('--det',     choices=['1', '2', 'both'], default='both',
                        help='Target DM-Ice detector (default: both, alternating)')
    args = parser.parse_args()

    rng        = np.random.default_rng(args.run)
    run_number = args.run
    n_events   = args.nevents

    output_dir = f"/data/user/bcharett/dmice_sim_output/run_{run_number:05d}/"
    os.makedirs(output_dir, exist_ok=True)

    # ── Choose target detector ──────────────────────────────────────────────
    if args.det == '1':
        target = 'det1'
    elif args.det == '2':
        target = 'det2'
    else:
        # alternate events between det1 and det2
        # generate two halves, interleave
        target = None

    # ── Generate injection file ─────────────────────────────────────────────
    inj_file = os.path.join(output_dir, f"dmice_injection_{run_number:05d}.h5")

    if target is not None:
        make_targeted_injection(inj_file, n_events, target, rng)
    else:
        # Split evenly: half aimed at det1, half at det2
        n1 = n_events // 2
        n2 = n_events - n1
        rng1 = np.random.default_rng(run_number * 2)
        rng2 = np.random.default_rng(run_number * 2 + 1)

        inj1 = inj_file.replace('.h5', '_det1.h5')
        inj2 = inj_file.replace('.h5', '_det2.h5')
        make_targeted_injection(inj1, n1, 'det1', rng1)
        make_targeted_injection(inj2, n2, 'det2', rng2)

        # Merge the two injection files into one
        import shutil
        shutil.copy(inj1, inj_file)   # start with det1
        # Append det2 events into the same file
        with h5py.File(inj_file, 'a') as fdst, h5py.File(inj2, 'r') as fsrc:
            grp_dst = fdst['VolumeInjector0']
            grp_src = fsrc['VolumeInjector0']
            for key in ['initial', 'final_1', 'final_2', 'properties']:
                combined = np.concatenate([grp_dst[key][:], grp_src[key][:]])
                del grp_dst[key]
                grp_dst.create_dataset(key, data=combined)
        os.remove(inj1)
        os.remove(inj2)
        print(f"[INFO] Merged det1+det2 injection → {inj_file}")

    # ── Configure Prometheus ────────────────────────────────────────────────
    print(f"\n[INFO] Configuring Prometheus (run {run_number}, {n_events} events)...")

    config["run"]["run number"]     = run_number
    config["run"]["nevents"]        = n_events
    config["run"]["storage prefix"] = output_dir

    # Full IC86 + DM-Ice geometry
    geofile = f"{resource_dir}geofiles/icecube_with_dmice.geo"
    config["detector"]["geo file"] = geofile
    print(f"[INFO] Geometry: {geofile}")

    # Use pre-generated injection file — skip LeptonInjector generation step.
    # Setting inject=False makes Prometheus read the file directly via
    # injection_from_LI_output(), bypassing LI so it cannot overwrite our file.
    config["injection"]["name"] = "LeptonInjector"
    config["injection"]["LeptonInjector"]["inject"] = False
    config["injection"]["LeptonInjector"]["paths"]["injection file"] = inj_file

    # PPC photon propagator (CPU — no CUDA on Condor nodes for this binary)
    config['photon propagator']['name'] = 'PPC_CUDA'
    config['photon propagator']['PPC_CUDA']['paths']['ppc_exe'] = \
        os.path.expanduser('~/prometheus/resources/PPC_executables/PPC_CUDA/ppc')
    ppc_tmpdir = os.path.join(output_dir, 'ppc_tmp')
    # Pre-clean so Prometheus's os.rmdir doesn't fail on a non-empty dir
    import shutil
    if os.path.exists(ppc_tmpdir):
        shutil.rmtree(ppc_tmpdir)
    os.makedirs(ppc_tmpdir, exist_ok=True)
    config['photon propagator']['PPC_CUDA']['paths']['ppc_tmpdir'] = ppc_tmpdir
    config['photon propagator']['PPC_CUDA']['paths']['ppctables'] = \
        os.path.expanduser('~/prometheus/resources/PPC_tables/south_pole/')
    config['photon propagator']['PPC_CUDA']['paths']['force'] = True
    config['photon propagator']['PPC_CUDA']['simulation']['device'] = -1   # CPU
    config['photon propagator']['PPC_CUDA']['simulation']['supress_output'] = True

    print(f"[INFO] Output dir: {output_dir}")

    # ── Run simulation ──────────────────────────────────────────────────────
    from prometheus import Prometheus
    print("[INFO] Initializing Prometheus...")
    p = Prometheus(config)
    print("[INFO] Starting simulation...")
    p.sim()
    print(f"\n[INFO] Done! Parquet: {output_dir}{run_number}_photons.parquet")


if __name__ == "__main__":
    main()
