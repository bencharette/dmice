#!/usr/bin/env python3
"""
Re-run PoleMuonLinefit on IC pulses + DM-Ice pivot hit, compare to standard.

For each event:
  1. Get MC truth muon transit time at the closest DM-Ice detector
  2. Create a new pulse series with a synthetic DM-Ice hit appended
  3. Run linefit.simple on the combined pulses -> "DM-Ice pivot PoleMuonLinefit"
  4. Compare angular error vs MC truth for both standard and DM-Ice versions

Usage (inside IceTray env):
    python dmice_pivot_refit.py -i dmice_muons_filtered.i3 -g gcdfile.i3.zst \
        --plot pivot_refit_comparison.png
"""

import sys
import argparse
import numpy as np
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    HAS_MPL = True
except Exception:
    HAS_MPL = False

try:
    from icecube import icetray, dataio, dataclasses, simclasses, linefit
    from icecube.icetray import I3Tray
except ImportError:
    sys.exit("ERROR: Load IceTray environment first.")

# DM-Ice detector positions in IceCube coordinates (meters)
DMICE_POS = {
    "det1": np.array([ 31.25,  -72.93, -511.05]),
    "det2": np.array([-334.80, -424.50, -511.26]),
}

# DM-Ice OMKeys matching icecube_with_dmice.geo (strings 87-88, dom 1)
DMICE_OMKEY = {
    "det1": icetray.OMKey(87, 1),
    "det2": icetray.OMKey(88, 1),
}

C_M_NS = 0.2998  # speed of light m/ns

PULSE_KEYS = ["OnlineL2_CleanedMuonPulses", "SRTInIcePulses",
              "SplitInIcePulses", "InIcePulses"]

COMBINED_PULSES_KEY = "DMIcePivotPulses"
COMBINED_LINEFIT_KEY = "DMIcePivotLineFit"


def angular_diff_deg(d1, d2):
    dot = float(np.dot(d1, d2))
    return np.degrees(np.arccos(np.clip(dot, -1.0, 1.0)))


def closest_approach_distance(pos, direction, point):
    dp = point - pos
    proj = np.dot(dp, direction)
    closest = pos + proj * direction
    return np.linalg.norm(point - closest)


def get_primary_muon(frame):
    for key in ('I3MCTree', 'I3MCTree_preMuonProp'):
        if key in frame:
            tree = frame[key]
            primaries = tree.primaries
            if not primaries:
                continue
            best, best_e = None, 0
            for p in tree:
                if abs(p.type) in (13,) or p.type in (
                        dataclasses.I3Particle.MuMinus,
                        dataclasses.I3Particle.MuPlus):
                    if p.energy > best_e:
                        best, best_e = p, p.energy
            if best is not None:
                return best
            return primaries[0]
    return None


class AddDMIceHit(icetray.I3Module):
    """
    Appends a synthetic DM-Ice pulse at the MC truth transit time to the
    existing IC pulse series and writes it as COMBINED_PULSES_KEY.
    """
    def __init__(self, context):
        super().__init__(context)
        self.current_daq = None

    def DAQ(self, frame):
        self.current_daq = frame
        self.PushFrame(frame)

    def Physics(self, frame):
        muon = get_primary_muon(frame)
        if muon is None and self.current_daq is not None:
            muon = get_primary_muon(self.current_daq)

        # Find existing pulse series
        pulse_map = None
        used_key = None
        for key in PULSE_KEYS:
            if key in frame:
                raw = frame[key]
                pulse_map = raw.apply(frame) if hasattr(raw, 'apply') else raw
                used_key = key
                break

        if pulse_map is None or muon is None:
            self.PushFrame(frame)
            return

        mc_pos = np.array([muon.pos.x, muon.pos.y, muon.pos.z])
        mc_dir = np.array([muon.dir.x, muon.dir.y, muon.dir.z])

        # Find closest DM-Ice detector
        ca = {k: closest_approach_distance(mc_pos, mc_dir, v)
              for k, v in DMICE_POS.items()}
        closest = min(ca, key=ca.get)
        dm_pos = DMICE_POS[closest]
        dm_omkey = DMICE_OMKEY[closest]

        # Compute MC truth transit time at DM-Ice
        # Anchor to IC pulse time centroid to stay in same time frame
        x_list, t_list, w_list = [], [], []
        for omkey, pulses in pulse_map.items():
            for p in pulses:
                x_list.append(np.array([0., 0., 0.]))  # position not needed
                t_list.append(p.time)
                w_list.append(p.charge if p.charge > 0 else 1.0)

        if not t_list:
            self.PushFrame(frame)
            return

        w = np.array(w_list)
        t = np.array(t_list)
        # Charge-weighted centroid position (need spatial centroid)
        x_s, y_s, z_s = [], [], []
        for omkey, pulses in pulse_map.items():
            for p in pulses:
                charge = p.charge if p.charge > 0 else 1.0
                x_s.append(charge)  # placeholder - redo below
                break

        # Redo with geometry - but we don't have dom_pos here.
        # Use only time centroid + MC direction projection approach:
        # t_dm = t_bar + dot(dm_pos - r_bar, mc_dir) / c
        # We don't have r_bar without dom_pos, so use muon.time directly
        # but corrected to pulse time frame via: t_ic_ref = first pulse time
        # and t_muon_at_ic_ref via projection.
        # Simplest: use t_bar of pulses as IC reference, then project.
        # Position centroid requires dom_pos which isn't loaded here.
        # We'll pass t_dm from outside via frame key instead.
        t_dm_ns = frame.get("DMIceTransitTime_ns", None)
        if t_dm_ns is None:
            self.PushFrame(frame)
            return

        # Build new pulse series = existing + one DM-Ice pulse
        new_map = dataclasses.I3RecoPulseSeriesMap()
        for omkey, pulses in pulse_map.items():
            new_map[omkey] = pulses

        dm_pulse = dataclasses.I3RecoPulse()
        dm_pulse.time = t_dm_ns
        dm_pulse.charge = 1.0
        dm_pulse.flags = 0
        new_map[dm_omkey] = dataclasses.I3RecoPulseSeries([dm_pulse])

        frame[COMBINED_PULSES_KEY] = new_map
        self.PushFrame(frame)


class ComputeDMIceTransitTime(icetray.I3Module):
    """
    Computes MC truth DM-Ice transit time in the IC pulse time frame
    and stores it in the frame as DMIceTransitTime_ns.
    """
    def __init__(self, context):
        super().__init__(context)
        self.dom_pos = {}
        self.current_daq = None

    def Configure(self):
        pass

    def set_dom_pos(self, dom_pos):
        self.dom_pos = dom_pos

    def DAQ(self, frame):
        self.current_daq = frame
        self.PushFrame(frame)

    def Physics(self, frame):
        muon = get_primary_muon(frame)
        if muon is None and self.current_daq is not None:
            muon = get_primary_muon(self.current_daq)
        if muon is None:
            self.PushFrame(frame)
            return

        mc_pos = np.array([muon.pos.x, muon.pos.y, muon.pos.z])
        mc_dir = np.array([muon.dir.x, muon.dir.y, muon.dir.z])

        ca = {k: closest_approach_distance(mc_pos, mc_dir, v)
              for k, v in DMICE_POS.items()}
        closest = min(ca, key=ca.get)
        dm_pos = DMICE_POS[closest]
        frame['DMIceClosestDet'] = dataclasses.I3String(closest)

        # Get IC pulse time centroid (charge-weighted)
        pulse_map = None
        for key in PULSE_KEYS:
            if key in frame:
                raw = frame[key]
                pulse_map = raw.apply(frame) if hasattr(raw, 'apply') else raw
                break

        if pulse_map is None or not self.dom_pos:
            self.PushFrame(frame)
            return

        x_list, y_list, z_list, t_list, w_list = [], [], [], [], []
        for omkey, pulses in pulse_map.items():
            if omkey not in self.dom_pos:
                continue
            pos = self.dom_pos[omkey]
            for p in pulses:
                charge = p.charge if p.charge > 0 else 1.0
                x_list.append(pos[0]); y_list.append(pos[1]); z_list.append(pos[2])
                t_list.append(p.time); w_list.append(charge)

        if not t_list:
            self.PushFrame(frame)
            return

        w = np.array(w_list)
        W = np.sum(w)
        r_bar = np.array([np.dot(w, x_list), np.dot(w, y_list), np.dot(w, z_list)]) / W
        t_bar = np.dot(w, t_list) / W

        # Distance from IC centroid to DM-Ice along MC track
        d = np.dot(dm_pos - r_bar, mc_dir)
        t_dm = t_bar + d / C_M_NS

        frame['DMIceTransitTime_ns'] = dataclasses.I3Double(t_dm)
        self.PushFrame(frame)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', required=True)
    parser.add_argument('-g', '--gcd',   required=True)
    parser.add_argument('--output', default='pivot_refit_results.csv')
    parser.add_argument('--plot',   default=None)
    args = parser.parse_args()

    # Load GCD for DOM positions
    print("Loading GCD...")
    dom_pos = {}
    f = dataio.I3File(args.gcd)
    while f.more():
        frame = f.pop_frame()
        if frame.Stop == icetray.I3Frame.Geometry:
            geo = frame['I3Geometry']
            for omkey, omgeo in geo.omgeo.items():
                pos = omgeo.position
                dom_pos[omkey] = np.array([pos.x, pos.y, pos.z])
            break
    f.close()
    print("  Loaded {} DOM positions".format(len(dom_pos)))

    # Add fake DM-Ice OMKeys to dom_pos so linefit.simple can find them
    for det, omkey in DMICE_OMKEY.items():
        dom_pos[omkey] = DMICE_POS[det]

    # Run IceTray pipeline
    outfile = '/tmp/dmice_pivot_refit_out.i3'
    tray = I3Tray()
    # Only prepend the GCD if it is a separate file; passing the same file
    # twice causes duplicate G frames mid-stream which corrupts the output.
    filelist = [args.input] if args.gcd == args.input else [args.gcd, args.input]
    tray.Add('I3Reader', filenamelist=filelist)

    # Step 1: compute DM-Ice transit time and store in frame
    def add_transit_time(frame):
        if frame.Stop != icetray.I3Frame.Physics:
            return True
        # inline version (module approach has state issues with dom_pos)
        muon = get_primary_muon(frame)
        if muon is None:
            return True
        mc_pos = np.array([muon.pos.x, muon.pos.y, muon.pos.z])
        mc_dir = np.array([muon.dir.x, muon.dir.y, muon.dir.z])
        ca = {k: closest_approach_distance(mc_pos, mc_dir, v)
              for k, v in DMICE_POS.items()}
        closest = min(ca, key=ca.get)
        dm_pos = DMICE_POS[closest]
        frame['DMIceClosestDet'] = dataclasses.I3String(closest)

        pulse_map = None
        for key in PULSE_KEYS:
            if key in frame:
                raw = frame[key]
                pulse_map = raw.apply(frame) if hasattr(raw, 'apply') else raw
                break
        if pulse_map is None:
            return True

        x_list, y_list, z_list, t_list, w_list = [], [], [], [], []
        for omkey, pulses in pulse_map.items():
            if omkey not in dom_pos:
                continue
            pos = dom_pos[omkey]
            for p in pulses:
                charge = p.charge if p.charge > 0 else 1.0
                x_list.append(pos[0]); y_list.append(pos[1]); z_list.append(pos[2])
                t_list.append(p.time); w_list.append(charge)
        if not t_list:
            return True

        w = np.array(w_list)
        W = np.sum(w)
        r_bar = np.array([np.dot(w, x_list), np.dot(w, y_list), np.dot(w, z_list)]) / W
        t_bar = np.dot(w, t_list) / W
        d = np.dot(dm_pos - r_bar, mc_dir)
        t_dm = t_bar + d / C_M_NS
        frame['DMIceTransitTime_ns'] = dataclasses.I3Double(t_dm)
        return True

    # Run standard LineFit first so PoleMuonLinefit exists for comparison.
    # (Real L2 data already has this; sim output does not.)
    tray.Add(linefit.simple, 'StandardLineFit',
             inputResponse='InIcePulses',
             fitName='PoleMuonLinefit')

    tray.Add(add_transit_time, 'ComputeTransitTime',
             Streams=[icetray.I3Frame.Physics])

    # Step 2: build combined pulse series with DM-Ice hit
    def build_combined_pulses(frame):
        if frame.Stop != icetray.I3Frame.Physics:
            return True
        if 'DMIceTransitTime_ns' not in frame:
            return True

        t_dm = frame['DMIceTransitTime_ns'].value
        closest = frame['DMIceClosestDet'].value
        dm_omkey = DMICE_OMKEY[closest]

        pulse_map = None
        for key in PULSE_KEYS:
            if key in frame:
                raw = frame[key]
                pulse_map = raw.apply(frame) if hasattr(raw, 'apply') else raw
                break
        if pulse_map is None:
            return True

        new_map = dataclasses.I3RecoPulseSeriesMap()
        for omkey, pulses in pulse_map.items():
            new_map[omkey] = pulses

        dm_pulse = dataclasses.I3RecoPulse()
        dm_pulse.time = float(t_dm)
        dm_pulse.charge = 1.0
        dm_pulse.flags = 0
        new_map[dm_omkey] = dataclasses.I3RecoPulseSeries([dm_pulse])
        frame[COMBINED_PULSES_KEY] = new_map
        return True

    tray.Add(build_combined_pulses, 'BuildCombinedPulses',
             Streams=[icetray.I3Frame.Physics])

    # Step 3: run linefit.simple on combined pulses
    tray.Add(linefit.simple, 'DMIcePivotLineFit',
             inputResponse=COMBINED_PULSES_KEY,
             fitName=COMBINED_LINEFIT_KEY)

    tray.Add('I3Writer', Filename=outfile,
             Streams=[icetray.I3Frame.Geometry,
                      icetray.I3Frame.Calibration,
                      icetray.I3Frame.DetectorStatus,
                      icetray.I3Frame.DAQ,
                      icetray.I3Frame.Physics])

    tray.Execute()
    tray.Finish()
    print("Wrote {}".format(outfile))

    # Read results and compare
    print("\nReading results...")
    results = []
    current_daq = None
    f = dataio.I3File(outfile)
    while f.more():
        frame = f.pop_frame()
        if frame.Stop == icetray.I3Frame.DAQ:
            current_daq = frame
            continue
        if frame.Stop != icetray.I3Frame.Physics:
            continue

        muon = get_primary_muon(frame)
        if muon is None and current_daq is not None:
            muon = get_primary_muon(current_daq)
        if muon is None:
            continue

        mc_dir = np.array([muon.dir.x, muon.dir.y, muon.dir.z])

        std_err = np.nan
        if 'PoleMuonLinefit' in frame:
            lf = frame['PoleMuonLinefit']
            std_dir = np.array([lf.dir.x, lf.dir.y, lf.dir.z])
            std_err = angular_diff_deg(mc_dir, std_dir)

        pivot_err = np.nan
        if COMBINED_LINEFIT_KEY in frame:
            lf2 = frame[COMBINED_LINEFIT_KEY]
            pivot_dir = np.array([lf2.dir.x, lf2.dir.y, lf2.dir.z])
            pivot_err = angular_diff_deg(mc_dir, pivot_dir)

        results.append(dict(std_err=std_err, pivot_err=pivot_err))

    f.close()

    import pandas as pd
    df = pd.DataFrame(results)
    df.to_csv(args.output, index=False)

    has = df['pivot_err'].notna() & df['std_err'].notna()
    improved = df.loc[has, 'pivot_err'] < df.loc[has, 'std_err']
    print("\n════════════════════════════════════════════════════════")
    print("Events with both recos:        {}".format(has.sum()))
    print("PoleMuonLinefit median error:   {:.2f} deg".format(df.loc[has,'std_err'].median()))
    print("DM-Ice Pivot LF median error:   {:.2f} deg".format(df.loc[has,'pivot_err'].median()))
    print("DM-Ice improves: {} / {} ({:.1f}%)".format(
        improved.sum(), has.sum(), 100*improved.mean()))
    print("════════════════════════════════════════════════════════")

    if args.plot and HAS_MPL:
        ic_err    = df['std_err'].dropna()
        pivot_err = df['pivot_err'].dropna()
        max_err   = max(ic_err.max(), pivot_err.max())
        bins = np.linspace(0, min(max_err * 1.05, 30), 46)

        fig, ax = plt.subplots(figsize=(9, 6))
        fig.suptitle('PoleMuonLinefit vs DM-Ice Pivot LineFit\n'
                     'Angular error vs MC truth  (n={} events)'.format(len(df)), fontsize=13)

        ax.hist(ic_err, bins=bins, histtype='stepfilled', alpha=0.5,
                color='steelblue', edgecolor='steelblue',
                label='PoleMuonLinefit  median={:.1f}°'.format(ic_err.median()))
        ax.axvline(ic_err.median(), color='navy', linewidth=2, linestyle='--')

        ax.hist(pivot_err, bins=bins, histtype='stepfilled', alpha=0.5,
                color='tomato', edgecolor='tomato',
                label='DM-Ice Pivot LineFit  median={:.1f}°'.format(pivot_err.median()))
        ax.axvline(pivot_err.median(), color='darkred', linewidth=2, linestyle='--')

        ax.set_xlabel('Angular error from MC truth (deg)', fontsize=12)
        ax.set_ylabel('Events', fontsize=12)
        ax.legend(fontsize=11)
        plt.tight_layout()
        plt.savefig(args.plot, dpi=150)
        print("Saved plot to {}".format(args.plot))


if __name__ == '__main__':
    main()
