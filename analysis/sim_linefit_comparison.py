#!/usr/bin/env python3
"""
Compare IC-only LineFit vs IC+DM-Ice pivot LineFit vs MC truth on simulated muons.

Supports BLO-format I3 files (from blo_npz_to_i3.py).  BLO files store the
momentum direction directly in primary.dir (no anti-momentum flip).

For each event:
  1. Get MC truth direction from I3MCTree primary
  2. Compute IC-only LineFit analytically from InIcePulses
  3. Compute IC+DM-Ice Pivot LineFit using MC truth DM-Ice transit time
  4. Compare angular errors of both reconstructions vs MC truth

Run inside IceTray environment:
    /cvmfs/icecube.opensciencegrid.org/py3-v4.3.0/RHEL_9_x86_64/metaprojects/icetray/v1.12.1/env-shell.sh \
        python sim_linefit_comparison.py -i blo_events.i3 -g gcdfile.i3.gz [--max-events N]
"""

import os
import sys
import argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

try:
    from icecube import icetray, dataio, dataclasses, simclasses
except ImportError:
    sys.exit("ERROR: Load IceTray environment first.")

# DM-Ice detector positions in IceCube coordinates (meters)
DMICE_POS = {
    "det1": np.array([ 31.25,  -72.93, -511.05]),
    "det2": np.array([-334.80, -424.50, -511.26]),
}

# Pulse series keys to try
PULSE_KEYS = ["SRTInIcePulses", "SplitInIcePulses", "InIcePulses",
              "TWCMuonPulseSeriesReco", "OfflinePulses",
              "OnlineL2_CleanedMuonPulses", "SplitUncleanedInIcePulses",
              "UncleanedInIcePulses"]

C_M_NS = 0.2998  # speed of light in m/ns
C_ICE_M_NS = 0.22  # c / n_ice


def run_ic_linefit(x, y, z, t, w):
    """Standard analytic IC-only LineFit (centre-of-gravity pivot)."""
    W = np.sum(w)
    if W == 0:
        return None
    x_bar = np.dot(w, x) / W
    y_bar = np.dot(w, y) / W
    z_bar = np.dot(w, z) / W
    t_bar = np.dot(w, t) / W
    dt = t - t_bar
    denom = np.dot(w, dt**2)
    if denom == 0:
        return None
    vx = np.dot(w * dt, x - x_bar) / denom
    vy = np.dot(w * dt, y - y_bar) / denom
    vz = np.dot(w * dt, z - z_bar) / denom
    speed = np.sqrt(vx**2 + vy**2 + vz**2)
    if speed == 0:
        return None
    dx, dy, dz = vx / speed, vy / speed, vz / speed
    zenith  = np.degrees(np.arccos(np.clip(-dz, -1.0, 1.0)))
    azimuth = np.degrees(np.arctan2(-dy, -dx) % (2 * np.pi))
    return dict(dx=dx, dy=dy, dz=dz, speed_m_ns=speed,
                zenith_deg=zenith, azimuth_deg=azimuth)


def angular_diff_deg(d1, d2):
    """Great-circle angle between two direction vectors, in degrees."""
    dot = float(np.dot(d1, d2))
    return np.degrees(np.arccos(np.clip(dot, -1.0, 1.0)))


def compute_dmice_hit_time(x_dom, y_dom, z_dom, t_dom, w_dom, dm_pos, mc_dir):
    """
    Compute the expected DM-Ice hit time in the same time frame as the IC pulses.

    Uses the charge-weighted IC centroid as a reference point, then projects
    (dm_pos - centroid) onto the MC truth direction to get the travel time offset.
    This is independent of the absolute DAQ time offset.
    """
    W = np.sum(w_dom)
    r_bar = np.array([np.dot(w_dom, x_dom), np.dot(w_dom, y_dom), np.dot(w_dom, z_dom)]) / W
    t_bar = np.dot(w_dom, t_dom) / W
    d = np.dot(dm_pos - r_bar, mc_dir)  # signed distance along track (m)
    return t_bar + d / C_M_NS           # same time frame as t_dom


def run_dmice_pivot_linefit_iterative(x_dom, y_dom, z_dom, t_dom, w_dom,
                                      x_dm, y_dm, z_dm, t_dm_ns,
                                      huber_sigma=200.0, max_iter=10):
    """
    Iterative DM-Ice pivot LineFit with Huber-style outlier down-weighting.

    Matches PoleMuonLinefit characteristics:
      1. Fit pivot LineFit to get direction + speed
      2. Compute time residual for each hit: t_residual = dt_i - dot(dr_i, v_hat)/speed
      3. Re-weight: w_new = w * min(1, huber_sigma / |residual|)
      4. Repeat until direction converges

    huber_sigma: residual scale in ns above which hits are down-weighted (~200 ns)
    """
    w = w_dom.copy()
    prev_dir = None

    for _ in range(max_iter):
        dt = t_dom - t_dm_ns
        dr_x = x_dom - x_dm
        dr_y = y_dom - y_dm
        dr_z = z_dom - z_dm

        denom = np.dot(w, dt**2)
        if denom == 0:
            return None
        vx = np.dot(w * dt, dr_x) / denom
        vy = np.dot(w * dt, dr_y) / denom
        vz = np.dot(w * dt, dr_z) / denom
        speed = np.sqrt(vx**2 + vy**2 + vz**2)
        if speed == 0:
            return None

        # Time residual: actual dt vs expected dt from track
        t_expected = (dr_x * vx + dr_y * vy + dr_z * vz) / speed**2
        residual = dt - t_expected

        # Huber re-weighting
        w = w_dom * np.minimum(1.0, huber_sigma / (np.abs(residual) + 1e-6))

        # Check convergence
        v_hat = np.array([vx, vy, vz]) / speed
        if prev_dir is not None:
            shift = np.degrees(np.arccos(np.clip(np.dot(v_hat, prev_dir), -1.0, 1.0)))
            if shift < 0.01:
                break
        prev_dir = v_hat

    dx, dy, dz = v_hat
    zenith  = np.degrees(np.arccos(np.clip(-dz, -1.0, 1.0)))
    azimuth = np.degrees(np.arctan2(-dy, -dx) % (2 * np.pi))
    return dict(dx=dx, dy=dy, dz=dz, speed_m_ns=speed,
                zenith_deg=zenith, azimuth_deg=azimuth)


def run_dmice_pivot_linefit(x_dom, y_dom, z_dom, t_dom, w_dom,
                            x_dm, y_dm, z_dm, t_dm_ns):
    """
    Pivot LineFit: the DM-Ice detector is the fixed reference in space AND time.

    Positions are measured from r_dm and times from t_dm, so the fit is
    constrained to pass through DM-Ice at t_dm.  This is the correct pivot
    formulation:

        dt_i = t_i - t_dm
        dr_i = r_i - r_dm
        v = Σ(w_i * dt_i * dr_i) / Σ(w_i * dt_i²)
    """
    dt = t_dom - t_dm_ns
    dr_x = x_dom - x_dm
    dr_y = y_dom - y_dm
    dr_z = z_dom - z_dm

    denom = np.dot(w_dom, dt**2)
    if denom == 0:
        return None

    vx = np.dot(w_dom * dt, dr_x) / denom
    vy = np.dot(w_dom * dt, dr_y) / denom
    vz = np.dot(w_dom * dt, dr_z) / denom

    speed = np.sqrt(vx**2 + vy**2 + vz**2)
    if speed == 0:
        return None

    dx, dy, dz = vx / speed, vy / speed, vz / speed
    zenith  = np.degrees(np.arccos(np.clip(-dz, -1.0, 1.0)))
    azimuth = np.degrees(np.arctan2(-dy, -dx) % (2 * np.pi))
    return dict(dx=dx, dy=dy, dz=dz, speed_m_ns=speed,
                zenith_deg=zenith, azimuth_deg=azimuth)


def closest_approach_distance(pos, direction, point):
    """Distance from a point to an infinite line defined by pos+t*direction."""
    dp = point - pos
    proj = np.dot(dp, direction)
    closest = pos + proj * direction
    return np.linalg.norm(point - closest)


def extract_hits(frame, dom_pos):
    """Extract charge-weighted pulses. Returns (x, y, z, t_abs_ns, w, n_doms)."""
    pulse_map = None
    for key in PULSE_KEYS:
        if key in frame:
            raw = frame[key]
            if hasattr(raw, 'apply'):
                pulse_map = raw.apply(frame)
            else:
                pulse_map = raw
            break

    if pulse_map is None:
        return None, None, None, None, None, 0

    event_start_ns = frame['I3EventHeader'].start_time.utc_daq_time * 0.1

    x_list, y_list, z_list, t_list, w_list = [], [], [], [], []
    n_doms = 0
    for omkey, pulses in pulse_map:
        if omkey not in dom_pos:
            continue
        x, y, z = dom_pos[omkey]
        dom_fired = False
        for pulse in pulses:
            charge = pulse.charge if pulse.charge > 0 else 1.0
            t_abs = event_start_ns + pulse.time
            x_list.append(x); y_list.append(y); z_list.append(z)
            t_list.append(t_abs); w_list.append(charge)
            dom_fired = True
        if dom_fired:
            n_doms += 1

    if not x_list:
        return None, None, None, None, None, 0

    return (np.array(x_list), np.array(y_list),
            np.array(z_list), np.array(t_list),
            np.array(w_list), n_doms)


def get_primary_muon(frame):
    """Get the primary muon from I3MCTree. Returns (direction, position, time) or None."""
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

    # Find the highest-energy muon
    best = None
    best_energy = 0
    for p in tree:
        if abs(p.type) in (13, dataclasses.I3Particle.MuMinus,
                           dataclasses.I3Particle.MuPlus):
            if p.energy > best_energy:
                best = p
                best_energy = p.energy

    if best is None:
        # Fall back to primary
        best = primaries[0]

    return best


def main():
    parser = argparse.ArgumentParser(
        description='Compare IC-only LineFit vs IC+DM-Ice LineFit vs MC truth')
    parser.add_argument('-i', '--input', required=True, help='Input i3 file')
    parser.add_argument('-g', '--gcd', default=None,
                        help='GCD file (optional; if omitted, geometry is read from the input file)')
    parser.add_argument('--max-events', type=int, default=0)
    parser.add_argument('--output', default='sim_linefit_results.csv')
    parser.add_argument('--plot', default=None, metavar='OUTPUT.png',
                        help='Save angular resolution histogram to this file')
    args = parser.parse_args()

    # Load geometry — from separate GCD or from first Geometry frame in input
    dom_pos = {}
    geo_source = args.gcd if args.gcd else args.input
    print("Loading geometry from {}...".format(geo_source))
    f = dataio.I3File(geo_source)
    while f.more():
        frame = f.pop_frame()
        if frame.Stop == icetray.I3Frame.Geometry:
            geo = frame['I3Geometry']
            for omkey, omgeo in geo.omgeo:
                pos = omgeo.position
                dom_pos[omkey] = np.array([pos.x, pos.y, pos.z])
            break
    f.close()
    print("  Loaded {} DOM positions".format(len(dom_pos)))

    # Process events
    print("Processing {}...".format(args.input))
    f = dataio.I3File(args.input)
    results = []
    n_total = 0
    n_with_linefit = 0
    n_near_dmice = 0

    current_daq = None
    while f.more():
        frame = f.pop_frame()
        if frame.Stop == icetray.I3Frame.DAQ:
            current_daq = frame
            continue
        if frame.Stop != icetray.I3Frame.Physics:
            continue

        n_total += 1

        # MC truth — I3MCTree lives in DAQ frame; try both Physics and cached DAQ
        muon = get_primary_muon(frame)
        if muon is None and current_daq is not None:
            muon = get_primary_muon(current_daq)
        if muon is None:
            continue

        # BLO files store momentum direction directly in primary.dir (no anti-momentum flip).
        mc_dir = np.array([muon.dir.x, muon.dir.y, muon.dir.z])   # travel direction
        mc_pos = np.array([muon.pos.x, muon.pos.y, muon.pos.z])
        mc_zenith = np.degrees(muon.dir.zenith)
        mc_azimuth = np.degrees(muon.dir.azimuth)
        mc_energy = muon.energy

        # DM-Ice detector selection: prefer BLO_DetId tag, else use closest approach
        if 'BLO_DetId' in frame:
            blo_det = str(frame['BLO_DetId'].value)
            closest_det = blo_det if blo_det in DMICE_POS else 'det1'
            ca_det1 = closest_approach_distance(mc_pos, mc_dir, DMICE_POS['det1'])
            ca_det2 = closest_approach_distance(mc_pos, mc_dir, DMICE_POS['det2'])
            ca_min = min(ca_det1, ca_det2)
        else:
            ca_det1 = closest_approach_distance(mc_pos, mc_dir, DMICE_POS['det1'])
            ca_det2 = closest_approach_distance(mc_pos, mc_dir, DMICE_POS['det2'])
            ca_min = min(ca_det1, ca_det2)
            closest_det = 'det1' if ca_det1 < ca_det2 else 'det2'

        if ca_min < 500:
            n_near_dmice += 1

        dm_pos = DMICE_POS[closest_det]

        # Compute IC-only LineFit analytically from pulses
        x, y, z, t_abs, w, n_doms = extract_hits(frame, dom_pos)
        if x is None or n_doms < 4:
            continue

        lf_ic = run_ic_linefit(x, y, z, t_abs, w)
        if lf_ic is None:
            continue
        n_with_linefit += 1

        ic_dir = np.array([lf_ic['dx'], lf_ic['dy'], lf_ic['dz']])
        ic_ang_err = angular_diff_deg(mc_dir, ic_dir)
        ic_speed = lf_ic['speed_m_ns']

        cfit_ang_err = np.nan
        cfit_ang_diff = np.nan
        cfit_speed = np.nan
        cfit_iter_ang_err = np.nan

        # mc_dir is travel direction (BLO convention); use directly for DM-Ice timing
        t_dm_ns = compute_dmice_hit_time(x, y, z, t_abs, w, dm_pos, mc_dir)

        cfit = run_dmice_pivot_linefit(
            x, y, z, t_abs, w,
            dm_pos[0], dm_pos[1], dm_pos[2], t_dm_ns,
        )
        if cfit is not None:
            cdir = np.array([cfit['dx'], cfit['dy'], cfit['dz']])
            cfit_ang_err = angular_diff_deg(mc_dir, cdir)
            cfit_ang_diff = angular_diff_deg(ic_dir, cdir)
            cfit_speed = cfit['speed_m_ns']

        cfit_iter = run_dmice_pivot_linefit_iterative(
            x, y, z, t_abs, w,
            dm_pos[0], dm_pos[1], dm_pos[2], t_dm_ns,
        )
        if cfit_iter is not None:
            cdir_iter = np.array([cfit_iter['dx'], cfit_iter['dy'], cfit_iter['dz']])
            cfit_iter_ang_err = angular_diff_deg(mc_dir, cdir_iter)

        results.append(dict(
            mc_zenith_deg=mc_zenith,
            mc_azimuth_deg=mc_azimuth,
            mc_energy_GeV=mc_energy,
            mc_ca_det1_m=ca_det1,
            mc_ca_det2_m=ca_det2,
            mc_ca_min_m=ca_min,
            closest_det=closest_det,
            ic_speed_m_ns=ic_speed,
            ic_ang_err_deg=ic_ang_err,
            cfit_ang_err_deg=cfit_ang_err,
            cfit_iter_ang_err_deg=cfit_iter_ang_err,
            cfit_ang_diff_deg=cfit_ang_diff,
            cfit_speed_m_ns=cfit_speed,
            n_doms=n_doms,
        ))

        if n_total % 100 == 0:
            print("  {} events processed...".format(n_total))

        if args.max_events and n_total >= args.max_events:
            break

    f.close()

    print("\n════════════════════════════════════════════════════════")
    print("Total Physics frames:    {}".format(n_total))
    print("With LineFit reco:       {}".format(n_with_linefit))
    print("MC track near DM-Ice:    {} (CA < 500m)".format(n_near_dmice))
    print("════════════════════════════════════════════════════════")

    if not results:
        print("No results — check that L2 processing produced LineFit.")
        return

    import pandas as pd
    df = pd.DataFrame(results)

    print("\n── IC-only LineFit angular error vs MC truth (deg) ─────")
    print(df['ic_ang_err_deg'].describe().to_string())

    has_cfit = df['cfit_ang_err_deg'].notna()
    if has_cfit.any():
        print("\n── DM-Ice Pivot LineFit angular error vs MC truth (deg) ──")
        print(df.loc[has_cfit, 'cfit_ang_err_deg'].describe().to_string())

        print("\n── Angular shift from adding DM-Ice pivot (deg) ────────")
        print(df.loc[has_cfit, 'cfit_ang_diff_deg'].describe().to_string())

        has_both = has_cfit & df['ic_ang_err_deg'].notna()
        improved = df.loc[has_both, 'cfit_ang_err_deg'] < df.loc[has_both, 'ic_ang_err_deg']
        print("\n── DM-Ice improves direction? ──────────────────────────")
        print("  {} / {} ({:.1f}%)".format(
            improved.sum(), has_both.sum(), 100 * improved.mean()))
        print("  IC-only LineFit median:          {:.2f} deg".format(
            df.loc[has_both, 'ic_ang_err_deg'].median()))
        print("  DM-Ice Pivot LineFit median:     {:.2f} deg".format(
            df.loc[has_cfit, 'cfit_ang_err_deg'].median()))
        has_iter = df['cfit_iter_ang_err_deg'].notna()
        if has_iter.any():
            print("  DM-Ice Pivot Iterative median:   {:.2f} deg".format(
                df.loc[has_iter, 'cfit_iter_ang_err_deg'].median()))

    # Break down by distance to DM-Ice
    print("\n── Angular error by distance to DM-Ice ─────────────────")
    bins = [0, 100, 200, 500, 1000, 5000]
    df['ca_bin'] = pd.cut(df['mc_ca_min_m'], bins)
    for ca_bin, grp in df.groupby('ca_bin', observed=True):
        has = grp['cfit_ang_err_deg'].notna()
        if has.any():
            print("  CA {:>15s}:  n={:3d}  IC={:.1f}°  IC+DM={:.1f}°  shift={:.1f}°".format(
                str(ca_bin), has.sum(),
                grp.loc[has, 'ic_ang_err_deg'].median(),
                grp.loc[has, 'cfit_ang_err_deg'].median(),
                grp.loc[has, 'cfit_ang_diff_deg'].median()))

    # Save
    df.to_csv(args.output, index=False)
    print("\nSaved {} rows to {}".format(len(df), args.output))

    # ── Overlaid angular resolution plot ──────────────────────────────────────
    if args.plot:
        ic_err        = df['ic_ang_err_deg'].dropna()
        cfit_err      = df['cfit_ang_err_deg'].dropna()
        cfit_iter_err = df['cfit_iter_ang_err_deg'].dropna()

        max_err = max(ic_err.max() if len(ic_err) else 0,
                      cfit_err.max() if len(cfit_err) else 0,
                      cfit_iter_err.max() if len(cfit_iter_err) else 0)
        bins = np.linspace(0, min(max_err * 1.05, 90), 46)

        fig, ax = plt.subplots(figsize=(9, 6))
        fig.suptitle(
            'Angular error vs MC truth  (n={} events)'.format(len(df)),
            fontsize=13)

        if len(ic_err) > 0:
            ax.hist(ic_err, bins=bins, histtype='stepfilled', alpha=0.5,
                    color='steelblue', edgecolor='steelblue',
                    label='IC-only LineFit  median={:.1f}°'.format(ic_err.median()))
            ax.axvline(ic_err.median(), color='navy', linewidth=2, linestyle='--')

        if len(cfit_err) > 0:
            ax.hist(cfit_err, bins=bins, histtype='stepfilled', alpha=0.4,
                    color='tomato', edgecolor='tomato',
                    label='DM-Ice Pivot LineFit  median={:.1f}°'.format(cfit_err.median()))
            ax.axvline(cfit_err.median(), color='darkred', linewidth=2, linestyle='--')

        if len(cfit_iter_err) > 0:
            ax.hist(cfit_iter_err, bins=bins, histtype='step', linewidth=2,
                    color='darkorange',
                    label='DM-Ice Pivot Iterative  median={:.1f}°'.format(cfit_iter_err.median()))
            ax.axvline(cfit_iter_err.median(), color='darkorange', linewidth=2, linestyle='--')

        ax.set_xlabel('Angular error from MC truth (deg)', fontsize=12)
        ax.set_ylabel('Events', fontsize=12)
        ax.legend(fontsize=11)
        plt.tight_layout()
        plt.savefig(args.plot, dpi=150)
        print("Saved histogram plot to {}".format(args.plot))


if __name__ == '__main__':
    main()
