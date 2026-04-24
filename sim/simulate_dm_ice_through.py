#!/usr/bin/env python3
"""
DMice Dark Matter Ice-Through Simulation
Simulates dark matter through-going events in IceCube using Prometheus.

"Ice-through" = ranged injection: DM enters the detector volume from outside,
travels through the ice, and exits. Useful for sensitivity studies on
DM-nucleon scattering or exotic stable heavy particles.
"""

import os
import prometheus
from prometheus import Prometheus, config

# Set up resource paths
prometheus_base = '/'.join(prometheus.__path__[0].split('/')[:-1])
resource_dir = f"{prometheus_base}/resources/"
output_dir = os.path.join(os.path.expanduser("~/dmice_work/output/"))
os.makedirs(output_dir, exist_ok=True)

print("[INFO] Setting up Prometheus DM ice-through simulation...")
print(f"[INFO] Resource directory: {resource_dir}")
print(f"[INFO] Output directory: {output_dir}")

# Configure run parameters
config["run"]["run_number"] = 2001
config["run"]["nevents"] = 100
config["run"]["storage_prefix"] = output_dir

print(f"[INFO] Run number: {config['run']['run_number']}")
print(f"[INFO] Events to simulate: {config['run']['nevents']}")

# Set detector
geofile = f"{resource_dir}geofiles/icecube.geo"
config["detector"]["geo file"] = geofile

print(f"[INFO] Detector: IceCube (ice-based)")
print(f"[INFO] Geometry file: {geofile}")

# Configure injection via LeptonInjector
injector = "LeptonInjector"
config["injection"]["name"] = injector
injection_config = config["injection"][injector]

# Ranged injection: DM enters from outside the detector volume and
# travels through the ice — the defining characteristic of ice-through events.
injection_config["simulation"]["injection_mode"] = "Ranged"

# Upgoing only (zenith 0-90 deg), matching the existing dataset which
# uses Earth-filtering to suppress atmospheric muon background.
degrees = 3.1415926536 / 180
injection_config["simulation"]["min_zenith"] = 0  * degrees
injection_config["simulation"]["max_zenith"] = 90 * degrees

print("[INFO] Injection configuration:")
print(f"  - Mode: Ranged (DM enters from outside, travels through ice)")
print(f"  - Directions: Upgoing only (zenith 0-90 deg, Earth-filtered)")

# Energy range matching the existing data (~400 TeV observed).
# Keeping the same range as the muon sim: 100 GeV - 1 PeV.
injection_config["simulation"]["minimal_energy"] = 1e2   # 100 GeV
injection_config["simulation"]["maximal_energy"] = 1e6   # 1 PeV
injection_config["simulation"]["gamma"] = 1

print(f"  - Energy range: 1e2 - 1e6 GeV (matches existing dataset range)")
print(f"  - Spectral index (gamma): {injection_config['simulation']['gamma']}")

# Final states matching the existing data exactly:
#   final_state_type: [13.0, -2000001010.0] = MuMinus + Hadrons (NuMu CC)
# This produces the same track topology as what you already have.
injection_config["simulation"]["final_state_1"] = "MuMinus"
injection_config["simulation"]["final_state_2"] = "Hadrons"

print(f"  - Final state 1: MuMinus  (PDG 13, matches existing data)")
print(f"  - Final state 2: Hadrons  (PDG -2000001010, matches existing data)")

# Run simulation
print("\n[INFO] Initializing Prometheus...")
p = Prometheus(config)

print("[INFO] Starting DM ice-through simulation...")
print("[INFO] This may take a few minutes...\n")

p.sim()

print("\n[INFO] Simulation complete!")
print(f"[INFO] Output files saved to: {output_dir}")
print("[INFO] Done!")
