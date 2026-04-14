#!/usr/bin/env python3
"""
plot_blo_truth_distributions.py

Read MC truth energy and zenith from a BLO-converted I3 file and
produce distribution plots + a CSV.

Usage (on Cobalt):
    /cvmfs/icecube.opensciencegrid.org/py3-v4.3.0/RHEL_9_x86_64/metaprojects/icetray/v1.12.1/env-shell.sh \
        python plot_blo_truth_distributions.py \
        -i ~/dmice_work/output/blo_dmice_targeted_det1det2_both_1000events.i3.zst \
        -o ~/dmice_work/output/blo_truth_distributions
"""

import argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from icecube import icetray, dataclasses, dataio


def get_primary_muon(frame):
    """Return the highest-energy muon from I3MCTree, or the first primary."""
    tree = None
    for key in ('I3MCTree', 'I3MCTree_preMuonProp'):
        if key in frame:
            tree = frame[key]
            break
    if tree is None:
        return None
    primaries = tree.primaries
    if not primaries:
        return None

    best = None
    best_energy = 0
    for p in tree:
        if abs(p.type) in (13, dataclasses.I3Particle.MuMinus,
                           dataclasses.I3Particle.MuPlus):
            if p.energy > best_energy:
                best = p
                best_energy = p.energy

    if best is None:
        best = primaries[0]
    return best


def main():
    parser = argparse.ArgumentParser(
        description='Plot MC truth energy and zenith distributions from BLO I3 file')
    parser.add_argument('-i', '--input', required=True, help='Input I3 file')
    parser.add_argument('-o', '--output-prefix', default='blo_truth_distributions',
                        help='Output prefix for .png and .csv (default: blo_truth_distributions)')
    args = parser.parse_args()

    zeniths = []
    energies = []
    n_frames = 0
    n_skipped = 0

    f = dataio.I3File(args.input)
    while f.more():
        frame = f.pop_frame()
        if frame.Stop != icetray.I3Frame.Physics:
            continue
        n_frames += 1

        muon = get_primary_muon(frame)
        if muon is None:
            n_skipped += 1
            continue

        zeniths.append(np.degrees(muon.dir.zenith))
        energies.append(muon.energy)

    f.close()

    zeniths = np.array(zeniths)
    energies = np.array(energies)

    print(f"Physics frames: {n_frames}, skipped (no muon): {n_skipped}, used: {len(zeniths)}")
    print(f"Zenith range:  {zeniths.min():.1f} – {zeniths.max():.1f} deg")
    print(f"Energy range:  {energies.min():.1f} – {energies.max():.1f} GeV")

    # Save CSV
    csv_path = args.output_prefix + '.csv'
    np.savetxt(csv_path,
               np.column_stack([zeniths, energies]),
               delimiter=',',
               header='mc_zenith_deg,mc_energy_GeV',
               comments='')
    print(f"Saved CSV: {csv_path}")

    # Plot
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle(f'BLO Simulated Muons — MC Truth  (n={len(zeniths)})', fontsize=13)

    # Zenith
    ax = axes[0]
    ax.hist(zeniths, bins=30, color='steelblue', edgecolor='white', linewidth=0.5)
    ax.set_xlabel('MC truth zenith (deg)', fontsize=12)
    ax.set_ylabel('Events', fontsize=12)
    ax.set_title('Zenith distribution')

    # Energy (log scale)
    ax = axes[1]
    log_bins = np.logspace(np.log10(max(energies.min(), 1)), np.log10(energies.max()), 30)
    ax.hist(energies, bins=log_bins, color='darkorange', edgecolor='white', linewidth=0.5)
    ax.set_xscale('log')
    ax.set_xlabel('MC truth energy (GeV)', fontsize=12)
    ax.set_ylabel('Events', fontsize=12)
    ax.set_title('Energy distribution')

    plt.tight_layout()
    png_path = args.output_prefix + '.png'
    plt.savefig(png_path, dpi=150)
    print(f"Saved plot: {png_path}")


if __name__ == '__main__':
    main()
