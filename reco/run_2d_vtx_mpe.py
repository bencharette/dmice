#!/usr/bin/env python3
"""
run_2d_vtx_mpe.py

2D vertex-constrained MPE reconstruction for DM-Ice coincidence events.

The DM-Ice NaI hit fixes the track vertex POSITION at the crystal centre.
The vertex TIME t₀ is derived from the IceCube seed vertex:

    t₀(θ,φ) = seed.time + dot(dm_pos - seed.pos, d_hat(θ,φ)) / C_VAC

This ensures t₀ is always consistent with the IceCube pulse timing.
The DM-Ice timing measurement then enters as a Gaussian soft constraint:

    log L_NaI = -0.5 * ((t₀ - (t_DM - 280 ns)) / 81 ns)²

Only zenith θ and azimuth φ are free; the combined objective is:

    log L = log L_MPE(IC DOMs) + log L_NaI

Run on Cobalt:
  /cvmfs/icecube.opensciencegrid.org/py3-v4.3.0/RHEL_9_x86_64/metaprojects/icetray/v1.12.1/env-shell.sh \\
    python3 -u ~/dmice/run_2d_vtx_mpe.py [--year 2012]

Outputs:
  ~/dmice_work/output/vtx2d_mpe.csv
"""

import os, sys, math, argparse
import numpy as np
from scipy.optimize import minimize
from scipy.special import gammaln, gammainc

# ── Physical constants ────────────────────────────────────────────────────────

C_VAC     = 0.2998           # speed of light in vacuum [m/ns]
N_ICE     = 1.3195           # ice refractive index (group)
C_ICE     = C_VAC / N_ICE
THETA_C   = math.acos(1.0 / N_ICE)
TAN_C     = math.tan(THETA_C)

# Pandel parameters — homogeneous IceCube ice
P_LAMBDA  = 47.0    # effective scattering length [m]
P_TAU     = 557.0   # scattering time scale [ns]
P_ABS     = 98.0    # absorption length [m]
RHO       = 1.0/P_TAU + C_ICE/P_ABS

# DM-Ice NaI timing model (real data)
MU_NS     = 280.0
SIGMA_NS  = 81.0

# DM-Ice crystal positions [m, IceCube coords]
DMICE_POS = {
    "det1": np.array([ 31.25,  -72.93, -511.05]),
    "det2": np.array([-334.80, -424.50, -511.26]),
}

IC_STRINGS   = set(range(1, 87))
MUON_STREAMS = {'', 'in_ice', 'InIceSplit'}

# Fit quality cuts
T_RES_MIN  = -15.0   # ns
T_RES_MAX  = 3000.0  # ns
D_PERP_MAX = 400.0   # m
D_PERP_MIN = 0.5     # m
N_DOMS_MIN = 5

# Seed priority
SEED_KEYS = ["RealIterMPE", "IterMPE", "MPEFit", "PoleMuonLlhFit",
             "RealPivotLF", "PivotLineFit", "LineFit", "PoleMuonLinefit"]

PULSE_PRIORITY = ["SplitInIcePulses", "OnlineL2_CleanedMuonPulses",
                  "OfflinePulses", "SRTInIcePulses",
                  "ReextractedInIcePulses", "InIcePulses"]

DEBUG_N = 0   # set via --debug N; print frame key info for first N events


# ── Pandel MPE + NaI combined log-likelihood ──────────────────────────────────

def combined_log_llh(theta, phi,
                     vtx_pos, seed_pos, seed_time, dm_t_corrected,
                     dom_xyz, dom_t_first, dom_charge):
    """
    Combined MPE Pandel (IC DOMs) + NaI Gaussian (DM-Ice) log-likelihood.

    Vertex position fixed at vtx_pos (crystal).
    Vertex time derived from IceCube seed:
        t₀(θ,φ) = seed_time + dot(vtx_pos - seed_pos, d_hat) / C_VAC

    DM-Ice constraint:
        log L_NaI = -0.5 * ((t₀ - dm_t_corrected) / SIGMA_NS)²

    IceCube convention: zenith=0 → downgoing → dz = -cos(θ).
    """
    sin_t = math.sin(theta);  cos_t = math.cos(theta)
    sin_p = math.sin(phi);    cos_p = math.cos(phi)
    d_hat = np.array([sin_t*cos_p, sin_t*sin_p, -cos_t])

    # Vertex time: how long along this track from seed vertex to crystal
    s_crystal = float(np.dot(vtx_pos - seed_pos, d_hat))
    t0_ns = seed_time + s_crystal / C_VAC

    # NaI Gaussian penalty
    dt_nai = t0_ns - dm_t_corrected
    log_nai = -0.5 * (dt_nai / SIGMA_NS)**2

    # Displacement from crystal to each DOM
    r      = dom_xyz - vtx_pos
    s      = r @ d_hat
    d_sq   = np.maximum(0.0, np.einsum('ij,ij->i', r, r) - s**2)
    d_perp = np.sqrt(d_sq)

    # Geometric first-photon arrival (Cherenkov cone)
    t_geo = t0_ns + (s - d_perp / TAN_C) / C_VAC
    t_res = dom_t_first - t_geo

    # Quality cuts
    mask = ((t_res  > T_RES_MIN) & (t_res  < T_RES_MAX) &
            (d_perp > D_PERP_MIN) & (d_perp < D_PERP_MAX))
    if mask.sum() < N_DOMS_MIN:
        return -np.inf

    t_r = t_res[mask]
    d_p = d_perp[mask]
    Npe = np.maximum(dom_charge[mask], 1.0)

    xi    = np.maximum(d_p / P_LAMBDA, 0.01)
    log_f = (xi*np.log(RHO) + (xi-1.0)*np.log(np.maximum(t_r, 0.1))
             - RHO*t_r - d_p/P_ABS - gammaln(xi))
    log_F = (np.log(np.maximum(gammainc(xi, RHO*np.maximum(t_r, 0.1)), 1e-300))
             - d_p/P_ABS)

    log_mpe = float(np.sum(np.log(Npe) + log_f + (Npe-1.0)*log_F))
    return log_mpe + log_nai


def fit_2d(vtx_pos, seed_pos, seed_time, dm_t_corrected,
           dom_xyz, dom_t_first, dom_charge, seed_zen, seed_azi):
    """Minimise −log L over (θ, φ). Returns (zen, azi, llh, n_doms)."""

    llh_seed = combined_log_llh(seed_zen, seed_azi,
                                vtx_pos, seed_pos, seed_time, dm_t_corrected,
                                dom_xyz, dom_t_first, dom_charge)
    if not np.isfinite(llh_seed):
        return np.nan, np.nan, -np.inf, 0

    def neg_llh(angles):
        v = combined_log_llh(angles[0], angles[1],
                             vtx_pos, seed_pos, seed_time, dm_t_corrected,
                             dom_xyz, dom_t_first, dom_charge)
        return -v if np.isfinite(v) else 1e9

    result = minimize(neg_llh, x0=[seed_zen, seed_azi], method='Nelder-Mead',
                      options=dict(xatol=1e-4, fatol=0.1, maxiter=2000, adaptive=True))

    if result.fun > 1e8:
        return np.nan, np.nan, -np.inf, 0

    zen = float(result.x[0]) % math.pi
    azi = float(result.x[1]) % (2*math.pi)
    llh = -float(result.fun)

    # Count DOMs used at the final solution
    sin_z = math.sin(zen); cos_z = math.cos(zen)
    sin_a = math.sin(azi); cos_a = math.cos(azi)
    d_hat = np.array([sin_z*cos_a, sin_z*sin_a, -cos_z])
    s_crystal = float(np.dot(vtx_pos - seed_pos, d_hat))
    t0 = seed_time + s_crystal / C_VAC
    r  = dom_xyz - vtx_pos
    s  = r @ d_hat
    d_sq   = np.maximum(0.0, np.einsum('ij,ij->i', r, r) - s**2)
    d_perp = np.sqrt(d_sq)
    t_res  = dom_t_first - (t0 + (s - d_perp/TAN_C)/C_VAC)
    n_doms = int(((t_res > T_RES_MIN) & (t_res < T_RES_MAX) &
                  (d_perp > D_PERP_MIN) & (d_perp < D_PERP_MAX)).sum())

    return zen, azi, llh, n_doms


# ── IceTray module ────────────────────────────────────────────────────────────

from icecube import icetray, dataio, dataclasses


class DMIce2DVertexMPE(icetray.I3Module):
    """2D vertex-constrained MPE fit anchored at DM-Ice crystal."""

    def __init__(self, ctx):
        super().__init__(ctx)
        self.AddParameter("OutputKey", "Frame key for output I3Particle", "VtxMPE2D")
        self.AddParameter("PulseKey",  "IC pulse series key (default: auto)", "")
        self.AddOutBox("OutBox")

    def Configure(self):
        self.output_key  = self.GetParameter("OutputKey")
        self.pulse_key   = self.GetParameter("PulseKey")
        self.om_pos      = {}
        self.debug_count = 0

    def Geometry(self, frame):
        self._load_geo(frame)
        self.PushFrame(frame)

    def _load_geo(self, frame):
        if "I3Geometry" not in frame:
            return
        geo = frame["I3Geometry"]
        new_pos = {}
        all_strings = set()
        for omk, omg in geo.omgeo.items():
            all_strings.add(omk.string)
            if omk.string in IC_STRINGS:
                new_pos[(omk.string, omk.om)] = np.array(
                    [omg.position.x, omg.position.y, omg.position.z])
        if new_pos:
            self.om_pos = new_pos
            print(f"[Geometry] strings: {min(all_strings)}-{max(all_strings)} "
                  f"({len(all_strings)} total), IC DOMs loaded={len(self.om_pos)}")
        else:
            print(f"[Geometry] WARNING: no IC DOMs found in geometry! "
                  f"strings present: {sorted(all_strings)[:10]}")

    def _get_seed(self, frame):
        for k in SEED_KEYS:
            if k in frame:
                p = frame[k]
                if (hasattr(p, "fit_status") and
                        p.fit_status == dataclasses.I3Particle.FitStatus.OK):
                    return p, k
        return None, None

    def _get_pulses(self, frame, debug=False):
        keys = ([self.pulse_key] if self.pulse_key else []) + PULSE_PRIORITY
        if debug:
            # show all frame keys that look like pulse series
            all_keys = list(frame.keys())
            pulse_like = [k for k in all_keys if any(
                w in k for w in ("Pulse","pulse","Hit","hit"))]
            print(f"  [debug] frame keys (pulse-like): {pulse_like}")
            print(f"  [debug] trying pulse keys: {keys}")
            print(f"  [debug] om_pos has {len(self.om_pos)} IC DOMs")
        for k in keys:
            if k not in frame:
                if debug:
                    print(f"  [debug]   {k}: not in frame")
                continue
            try:
                pmap   = dataclasses.I3RecoPulseSeriesMap.from_frame(frame, k)
                result = {}
                for omk, pulses in pmap:
                    if omk.string not in IC_STRINGS or not pulses:
                        continue
                    gk = (omk.string, omk.om)
                    if gk not in self.om_pos:
                        continue
                    result[gk] = (min(p.time for p in pulses),
                                  sum(p.charge for p in pulses))
                if debug:
                    print(f"  [debug]   {k}: found {len(result)} IC DOMs")
                if len(result) >= N_DOMS_MIN:
                    return result, k
            except Exception as e:
                if debug:
                    print(f"  [debug]   {k}: exception {e}")
                continue
        return {}, ""

    def Physics(self, frame):
        # Lazy geometry init — handles files where G frames come mixed with P frames
        if not self.om_pos:
            self._load_geo(frame)

        hdr    = frame["I3EventHeader"]
        stream = getattr(hdr, "sub_event_stream", "")
        if stream not in MUON_STREAMS:
            self.PushFrame(frame); return

        if "DMIce_detection_time" not in frame:
            self.PushFrame(frame); return

        det_str = str(frame["DMIce_detector"]) if "DMIce_detector" in frame else "det1"
        det_key = "det1" if "det1" in det_str else "det2"
        vtx_pos = DMICE_POS[det_key]

        event_start_daq = hdr.start_time.utc_daq_time
        dm_t_ns         = (frame["DMIce_detection_time"].value - event_start_daq) * 0.1
        dm_t_corrected  = dm_t_ns - MU_NS   # what t₀ should equal for signal events

        seed, seed_key = self._get_seed(frame)
        if seed is None:
            self.PushFrame(frame); return

        seed_pos  = np.array([seed.pos.x, seed.pos.y, seed.pos.z])
        seed_time = float(seed.time)
        seed_zen  = float(seed.dir.zenith)
        seed_azi  = float(seed.dir.azimuth)

        do_debug = (self.debug_count < DEBUG_N)
        if do_debug:
            self.debug_count += 1
            print(f"\n[debug event {self.debug_count}/{DEBUG_N}] "
                  f"run={hdr.run_id} evt={hdr.event_id} stream='{stream}'")
        pulses, pulse_key_used = self._get_pulses(frame, debug=do_debug)
        if do_debug:
            print(f"  [debug] pulse key used: '{pulse_key_used}', ndoms={len(pulses)}")
        if len(pulses) < N_DOMS_MIN:
            self.PushFrame(frame); return

        dom_xyz     = np.array([self.om_pos[k] for k in pulses])
        dom_t_first = np.array([v[0]           for v in pulses.values()])
        dom_charge  = np.array([v[1]           for v in pulses.values()])

        zen, azi, llh, n_doms = fit_2d(
            vtx_pos, seed_pos, seed_time, dm_t_corrected,
            dom_xyz, dom_t_first, dom_charge, seed_zen, seed_azi)

        p_out = dataclasses.I3Particle()
        if np.isfinite(zen) and np.isfinite(azi):
            p_out.dir      = dataclasses.I3Direction(zen, azi)
            p_out.pos      = dataclasses.I3Position(*vtx_pos)
            p_out.time     = seed_time + float(np.dot(vtx_pos - seed_pos,
                             np.array([math.sin(zen)*math.cos(azi),
                                       math.sin(zen)*math.sin(azi),
                                       -math.cos(zen)]))) / C_VAC
            p_out.speed    = dataclasses.I3Constants.c
            p_out.shape    = dataclasses.I3Particle.ParticleShape.InfiniteTrack
            p_out.fit_status = dataclasses.I3Particle.FitStatus.OK
        else:
            p_out.fit_status = dataclasses.I3Particle.FitStatus.FailedToConverge

        # dm_dt = DM-Ice corrected time minus IceCube-derived t0 at the crystal,
        # evaluated at the FITTED direction (or seed if fit failed).
        # Near 0 → signal; large → accidental.
        if np.isfinite(zen) and np.isfinite(azi):
            fit_zen, fit_azi = zen, azi
        else:
            fit_zen, fit_azi = seed_zen, seed_azi
        d_fit = np.array([math.sin(fit_zen)*math.cos(fit_azi),
                          math.sin(fit_zen)*math.sin(fit_azi),
                          -math.cos(fit_zen)])
        t0_fit = seed_time + float(np.dot(vtx_pos - seed_pos, d_fit)) / C_VAC
        dm_dt  = dm_t_corrected - t0_fit

        frame[self.output_key]              = p_out
        frame[self.output_key + "_llh"]     = dataclasses.I3Double(llh)
        frame[self.output_key + "_ndoms"]   = icetray.I3Int(n_doms)
        frame[self.output_key + "_seed"]    = dataclasses.I3String(seed_key or "")
        frame[self.output_key + "_dm_dt"]   = dataclasses.I3Double(dm_dt)
        self.PushFrame(frame)


# ── Scorer ────────────────────────────────────────────────────────────────────

class ScorerModule(icetray.I3Module):
    def __init__(self, ctx):
        super().__init__(ctx)
        self.AddParameter("OutputCSV", "Path to output CSV", "")
        self.AddOutBox("OutBox")

    def Configure(self):
        self.csv_path = self.GetParameter("OutputCSV")
        self.records  = []
        self.seen     = set()

    def _zen(self, frame, key):
        if key in frame:
            p = frame[key]
            if (hasattr(p, "fit_status") and
                    p.fit_status == dataclasses.I3Particle.FitStatus.OK):
                return math.degrees(p.dir.zenith), math.degrees(p.dir.azimuth)
        return np.nan, np.nan

    def Physics(self, frame):
        hdr    = frame["I3EventHeader"]
        stream = getattr(hdr, "sub_event_stream", "")
        if stream not in MUON_STREAMS:
            self.PushFrame(frame); return

        dk = (hdr.run_id, hdr.event_id, stream)
        if dk in self.seen:
            self.PushFrame(frame); return
        self.seen.add(dk)

        vtx_zen, vtx_azi = self._zen(frame, "VtxMPE2D")
        llh    = frame["VtxMPE2D_llh"].value   if "VtxMPE2D_llh"   in frame else np.nan
        ndoms  = frame["VtxMPE2D_ndoms"].value  if "VtxMPE2D_ndoms" in frame else 0
        dm_dt  = frame["VtxMPE2D_dm_dt"].value  if "VtxMPE2D_dm_dt" in frame else np.nan
        seed_k = str(frame["VtxMPE2D_seed"])    if "VtxMPE2D_seed"  in frame else ""

        lf_zen,   lf_azi   = self._zen(frame, "LineFit")
        if np.isnan(lf_zen):
            lf_zen, lf_azi = self._zen(frame, "PoleMuonLinefit")
        mpe_zen,  mpe_azi  = self._zen(frame, "MPEFit")
        iter_zen, iter_azi = self._zen(frame, "RealIterMPE")
        if np.isnan(iter_zen):
            iter_zen, iter_azi = self._zen(frame, "IterMPE")

        dm_t_ns = np.nan
        if "DMIce_detection_time" in frame:
            dm_t_ns = (frame["DMIce_detection_time"].value
                       - hdr.start_time.utc_daq_time) * 0.1

        self.records.append(dict(
            year        = hdr.start_time.utc_year,
            run_id      = hdr.run_id,
            event_id    = hdr.event_id,
            vtx2d_zen   = vtx_zen,
            vtx2d_azi   = vtx_azi,
            vtx2d_llh   = llh,
            vtx2d_ndoms = ndoms,
            vtx2d_seed  = seed_k,
            vtx2d_dm_dt = dm_dt,   # t₀_seed - t₀_DM-Ice [ns]; near 0 for signal
            lf_zen      = lf_zen,
            mpe_zen     = mpe_zen,
            iter_zen    = iter_zen,
            dm_t_ns     = dm_t_ns,
        ))
        self.PushFrame(frame)

    def Finish(self):
        import pandas as pd
        df = pd.DataFrame(self.records)
        df.to_csv(self.csv_path, index=False)
        n_ok = df["vtx2d_zen"].notna().sum()
        print(f"\nDone: {len(df)} events, {n_ok} with valid VtxMPE2D fit")
        print(f"CSV: {self.csv_path}")
        if n_ok > 0:
            print(f"VtxMPE2D zenith median: {df['vtx2d_zen'].median():.1f}°")
            print(f"vtx2d_dm_dt median: {df['vtx2d_dm_dt'].median():.0f} ns "
                  f"(near 0 = signal, large = accidental)")


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--year",   type=int, default=None)
    parser.add_argument("--i3",     default="/data/user/bcharett/dmice_coincidences_2011_2022/"
                                            "all_dmice_coincidences_2011_2022_fixed.i3")
    parser.add_argument("--gcd",    default="/cvmfs/icecube.opensciencegrid.org/data/GCD/"
                                            "GeoCalibDetectorStatus_2013.56429_V1.i3.gz",
                        help="GCD file to prepend (provides I3Geometry)")
    parser.add_argument("--out",    default=os.path.expanduser(
                                            "~/dmice_work/output/vtx2d_mpe.csv"))
    parser.add_argument("--pulses", default="")
    parser.add_argument("--debug",  type=int, default=0,
                        help="Print frame/pulse diagnostics for first N events")
    args = parser.parse_args()

    global DEBUG_N
    DEBUG_N = args.debug

    from icecube.icetray import I3Tray

    print("=" * 60)
    print("DM-Ice 2D Vertex-Constrained MPE Fit")
    print("=" * 60)
    print(f"Input:  {args.i3}")
    print(f"GCD:    {args.gcd}")
    print(f"Output: {args.out}")
    print(f"Year:   {args.year or 'all'}")
    print(f"Pandel: λ={P_LAMBDA}m  τ={P_TAU}ns  λ_a={P_ABS}m")
    print(f"Vertex position fixed at DM-Ice crystal")
    print(f"Vertex time from IceCube seed + Gaussian NaI penalty (σ={SIGMA_NS}ns)")
    print()

    tray = I3Tray()
    file_list = ([args.gcd] if args.gcd else []) + [args.i3]
    tray.Add("I3Reader", FilenameList=file_list)

    if args.year:
        yr = args.year
        tray.Add(lambda f: (f.Stop != icetray.I3Frame.Physics or
                            f["I3EventHeader"].start_time.utc_year == yr))

    tray.Add(DMIce2DVertexMPE, OutputKey="VtxMPE2D", PulseKey=args.pulses)
    tray.Add(ScorerModule, OutputCSV=args.out)

    tray.Execute()
    tray.Finish()


if __name__ == "__main__":
    main()
