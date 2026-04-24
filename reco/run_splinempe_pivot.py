#!/usr/bin/env python3
"""
run_splinempe_pivot.py

Runs SplineMPE with two seeds on BLO simulation and compares using DM-Ice log L:
  A) Standard seed: LineFit on IC-only pulses
  B) Pivot seed:    LineFit anchored to DM-Ice transit time

Tray structure:
  NPZInjector → SplitPulses (IC / DM-Ice) → LineFit → PivotLF
              → SplineMPE(std seed) + SplineMPE(piv seed) → Score + CSV

Input:  ~/dmice_work/output/blo_dmice_targeted_det1det2_both_1000events_repacked.npz
Model:  ~/dmice_work/output/dmice_timing_model.npz
Output: ~/dmice_work/output/splinempe_pivot_comparison.{csv,png}

Run on Cobalt:
  /cvmfs/.../icetray/v1.12.1/env-shell.sh python3 ~/dmice/run_splinempe_pivot.py
"""

import os, csv, math
import numpy as np

# ── Paths ─────────────────────────────────────────────────────────────────────

NPZ_FILE   = os.path.expanduser(
    "~/dmice_work/output/muons_binned_5000ev_repacked_injected.npz")
MODEL_FILE = os.path.expanduser("~/dmice_work/output/dmice_timing_model.npz")
GEO_FILE   = os.path.expanduser(
    "~/dmice/BlueLightOrchestra.jl/resources/geofiles/icecube_with_dmice.geo")
OUT_BASE   = os.path.expanduser("~/dmice_work/output/splinempe_pivot_comparison")

SPLINE_PROB = ("/cvmfs/icecube.opensciencegrid.org/data/photon-tables/splines/"
               "InfBareMu_mie_prob_z20a10_V2.fits")
SPLINE_AMP  = ("/cvmfs/icecube.opensciencegrid.org/data/photon-tables/splines/"
               "InfBareMu_mie_abs_z20a10_V2.fits")

# ── Constants ─────────────────────────────────────────────────────────────────

Z_OFFSET = 1948.07     # BLO z → IceCube z  [m]
C_M_NS   = 0.2998      # speed of light in vacuum [m/ns]
N_ICE    = 1.3195
THETA_C  = math.acos(1.0 / N_ICE)

DMICE_OMKEYS = {0: (87, 1), 1: (88, 1)}
DMICE_POS_IC = {
    0: np.array([ 31.25,  -72.93, -511.05]),
    1: np.array([-334.80, -424.50, -511.26]),
}

# Frame keys
IC_PULSES  = "InIcePulses"       # IC-only (strings 1–86), used for SplineMPE
ALL_PULSES = "AllPulses"         # IC + DM-Ice, used for Pivot LF
DM_T_KEY   = "DMIce_t"          # DM-Ice first-photon time [ns] (I3Double)
DM_ID_KEY  = "DMIce_id"         # DM-Ice det id (I3Int, 0=det1, 1=det2)
LF_KEY       = "LineFit"
PIV_LF_KEY   = "PivotLineFit"
SMPE_STD     = "SplineMPE_Std"
SMPE_PIV     = "SplineMPE_Pivot"
MPE_STD      = "MPEFit_Std"
MPE_PIV      = "MPEFit_Pivot"     # pivot seeded from LineFit direction
PIV_MPE_KEY  = "PivotMPE_LF"     # pivot anchor computed using MPEFit direction
MPE_PIV2     = "MPEFit_Pivot2"   # MPEFit seeded from MPEFit-direction pivot
SPE_STD      = "SPEFit_Std"
SPE_PIV      = "SPEFit_Pivot"

# ── Load timing model ─────────────────────────────────────────────────────────

_m      = np.load(MODEL_FILE, allow_pickle=True)
MU_NS   = float(_m["mu_ns"])
SIGMA_NS= float(_m["sigma_ns"])
EPS0    = float(_m["efficiency"])
D_MAX   = float(_m["d_perp_max_m"])
print(f"Timing model: μ={MU_NS:+.1f}ns  σ={SIGMA_NS:.1f}ns  ε={EPS0:.3f}  d⊥<{D_MAX:.0f}m")

# ── IceTray imports ───────────────────────────────────────────────────────────

from icecube import icetray, dataclasses, dataio
from icecube import linefit, spline_reco, lilliput, gulliver, gulliver_modules
import icecube.lilliput.segments
from icecube.icetray import I3Units, I3Tray

# ── Load geometry ─────────────────────────────────────────────────────────────

def load_geo(path):
    doms = {}
    with open(path) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 5:
                continue
            try:
                x, y, z_dep = float(parts[0]), float(parts[1]), float(parts[2])
                s, dom = int(parts[3]), int(parts[4])
                doms[(s, dom)] = (x, y, z_dep + Z_OFFSET)
            except ValueError:
                continue
    return doms

geo_doms = load_geo(GEO_FILE)

geo_obj = dataclasses.I3Geometry()
for (s, dom), (px, py, pz) in geo_doms.items():
    omkey           = icetray.OMKey(s, dom)
    omgeo           = dataclasses.I3OMGeo()
    omgeo.position  = dataclasses.I3Position(px, py, pz)
    omgeo.omtype    = dataclasses.I3OMGeo.IceCube
    geo_obj.omgeo[omkey] = omgeo

print(f"Geometry loaded: {len(geo_doms)} DOMs (incl. DM-Ice strings 87/88)")

# ── Load NPZ ──────────────────────────────────────────────────────────────────

d  = np.load(NPZ_FILE, allow_pickle=True)
N  = len(d["energy_GeV"])
print(f"NPZ: {N} events from {os.path.basename(NPZ_FILE)}")

def load_ragged(key):
    if f"{key}_flat" in d:
        flat    = d[f"{key}_flat"]
        offsets = d[f"{key}_offsets"]
        return [flat[offsets[i]:offsets[i+1]] for i in range(N)]
    return list(d[key])

_dom_x      = load_ragged("dom_x")
_dom_y      = load_ragged("dom_y")
_dom_z      = load_ragged("dom_z")
_dom_t      = load_ragged("dom_t")
_dom_nhits  = load_ragged("dom_nhits")
_dom_string = load_ragged("dom_string")
_dom_sensor = load_ragged("dom_sensor")

# ── Helper: t_geometric and DM-Ice log L ─────────────────────────────────────

def t_geometric(track_pos, track_dir, t0, dm_pos):
    r     = np.asarray(dm_pos) - np.asarray(track_pos)
    d_hat = np.asarray(track_dir, dtype=float)
    d_hat = d_hat / np.linalg.norm(d_hat)
    s     = np.dot(r, d_hat)
    d_perp = math.sqrt(max(0.0, np.dot(r, r) - s**2))
    t_pca  = t0 + s / C_M_NS
    t_geo  = t_pca if d_perp < 0.01 else t_pca + d_perp / (C_M_NS * math.sin(THETA_C))
    return t_geo, d_perp


def dm_log_l(track_pos, track_dir, t0, dm_t, dm_pos, d_perp_override=None):
    t_geo, d_perp_geom = t_geometric(track_pos, track_dir, t0, dm_pos)
    d_perp = d_perp_override if d_perp_override is not None else d_perp_geom
    if d_perp > D_MAX:
        return float("nan"), dm_t - t_geo, d_perp
    dt = dm_t - t_geo
    ll = (-0.5 * ((dt - MU_NS) / SIGMA_NS)**2
          - math.log(SIGMA_NS * math.sqrt(2 * math.pi))
          + math.log(max(EPS0, 1e-9)))
    return ll, dt, d_perp


def ang_err_deg(d1, d2):
    dot = max(-1.0, min(1.0, float(
        np.dot(np.asarray(d1), np.asarray(d2)))))
    return math.degrees(math.acos(abs(dot)))


def pivot_linefit_ic(xs, ys, zs, ts, ws, dm_pos, dm_t_corrected, seed_dir):
    """Pivot LineFit anchored to μ-corrected DM-Ice transit time.

    dm_t_corrected: DM-Ice hit time with timing model μ already subtracted.
    Returns (direction unit vector, t0_pivot) or None on failure.
    t0_pivot is the event time at the DM-Ice position (for seeding MPEFit vertex).
    """
    ws = np.asarray(ws, dtype=float)
    W  = ws.sum()
    if W == 0:
        return None
    cx = np.dot(ws, xs) / W; cy = np.dot(ws, ys) / W; cz = np.dot(ws, zs) / W
    tb = np.dot(ws, ts) / W
    d_proj = ((dm_pos[0]-cx)*seed_dir[0] + (dm_pos[1]-cy)*seed_dir[1]
              + (dm_pos[2]-cz)*seed_dir[2])
    # Use μ-corrected DM-Ice time as the anchor (removes NaI scintillation bias)
    t_dm = dm_t_corrected
    dts  = ts - t_dm
    drxs = xs - dm_pos[0]; drys = ys - dm_pos[1]; drzs = zs - dm_pos[2]
    den  = np.dot(ws, dts**2)
    if den < 1e-10:
        return None
    vx = np.dot(ws * dts, drxs) / den
    vy = np.dot(ws * dts, drys) / den
    vz = np.dot(ws * dts, drzs) / den
    spd = math.sqrt(vx**2 + vy**2 + vz**2)
    if spd < 1e-6:
        return None
    # Disambiguate direction: ensure consistent with seed
    if vx*seed_dir[0] + vy*seed_dir[1] + vz*seed_dir[2] < 0:
        vx, vy, vz = -vx, -vy, -vz
    direction = np.array([vx, vy, vz]) / spd
    # t0 at DM-Ice position: track was at dm_pos at time dm_t_corrected
    # back-project to a standard reference using the fitted direction
    t0_pivot = dm_t_corrected
    return direction, t0_pivot


# ── Combined IC + DM-Ice likelihood (scipy) ───────────────────────────────────

from scipy.optimize import minimize as scipy_minimize
from scipy.special import gammaln

# Approximate uniform ice model (SPICEMie bulk values)
_PANDEL_LA  = 98.0    # absorption length [m]
_PANDEL_LS  = 30.0    # scattering length [m]
_JITTER_NS  = 15.0    # DOM timing jitter [ns]

def _pandel_log_spe(t_res, d_perp):
    """Log SPE Pandel PDF (conditional on a hit). Gamma distribution in t_res."""
    d = max(d_perp, 1.0)
    alpha = d / _PANDEL_LS          # Gamma shape
    beta  = C_M_NS / _PANDEL_LA    # Gamma rate [1/ns]
    if t_res < 0:
        # Gaussian jitter smearing for early photons
        return -0.5 * (t_res / _JITTER_NS)**2 - math.log(_JITTER_NS * math.sqrt(2*math.pi))
    return (alpha * math.log(beta) + (alpha - 1) * math.log(t_res + 1e-6)
            - beta * t_res - gammaln(alpha))


def _ic_log_l(zen, azi, t0, vertex, ic_pulse_list):
    """Sum of SPE Pandel log L over all IC pulses. vertex = (x,y,z) reference point."""
    sin_z, cos_z = math.sin(zen), math.cos(zen)
    sin_a, cos_a = math.sin(azi), math.cos(azi)
    dx = sin_z * cos_a
    dy = sin_z * sin_a
    dz = -cos_z   # IceCube: zenith=0 → downgoing → dz<0
    vx, vy, vz = vertex
    ll = 0.0
    for (px, py, pz, t_hit, charge) in ic_pulse_list:
        rx, ry, rz = px - vx, py - vy, pz - vz
        s = rx*dx + ry*dy + rz*dz
        d_perp2 = rx*rx + ry*ry + rz*rz - s*s
        d_perp = math.sqrt(max(0.0, d_perp2))
        t_geo = t0 + s / C_M_NS
        if d_perp > 0.01:
            t_geo += d_perp / (C_M_NS * math.sin(THETA_C))
        ll += charge * _pandel_log_spe(t_hit - t_geo, d_perp)
    return ll


def _neg_combined_ll(params, vertex, dm_pos, dm_t_corrected, ic_pulse_list, sigma):
    """Negative combined log L: IC Pandel + DM-Ice Gaussian. For scipy.minimize.

    vertex: track reference point (kept at seed position — NOT dm_pos)
    dm_pos: DM-Ice detector position
    dm_t_corrected: observed DM-Ice hit time minus μ_NaI (280 ns pivot correction)

    The DM-Ice expected time is computed from the track geometry:
        t_geo_DM = t0 + (r_DM - vertex)·d̂ / c
    This keeps the IC Pandel vertex geometry intact.
    """
    zen, azi, t0 = params
    if not (0.0 < zen < math.pi):
        return 1e9
    sin_z, cos_z = math.sin(zen), math.cos(zen)
    sin_a, cos_a = math.sin(azi), math.cos(azi)
    dx = sin_z * cos_a; dy = sin_z * sin_a; dz = -cos_z
    # Expected time muon reaches dm_pos along this track
    vx, vy, vz = vertex
    rx = dm_pos[0] - vx; ry = dm_pos[1] - vy; rz = dm_pos[2] - vz
    s_dm = rx*dx + ry*dy + rz*dz
    t_geo_dm = t0 + s_dm / C_M_NS
    ll_ic = _ic_log_l(zen, azi, t0, vertex, ic_pulse_list)
    ll_dm = -0.5 * ((dm_t_corrected - t_geo_dm) / sigma)**2
    return -(ll_ic + ll_dm)


MPE_DM_KEY    = "MPEFit_DM"       # MPEFit + DM-Ice combined likelihood (std seed)
SPE_DM_KEY    = "SPEFit_DM"       # SPEFit + DM-Ice combined likelihood (std seed)
PIV_SPAT_KEY  = "PivotLF_Spatial" # LineFit vertex at dm_pos, NO timing correction
MPE_SPAT_KEY  = "MPEFit_Spatial"  # MPEFit seeded from spatial-only pivot
SPE_SPAT_KEY  = "SPEFit_Spatial"  # SPEFit seeded from spatial-only pivot


class SpatialPivotModule(icetray.I3Module):
    """
    Spatial-only DM-Ice pivot: re-anchors track vertex to dm_pos using IC timing only.
    Does NOT use DM-Ice hit time — only the POSITION of the crystal (d_perp≈0).
    Seeds MPEFit and SPEFit to isolate the contribution of spatial vs temporal info.

    Applied to ALL events: every event in this dataset is a DM-Ice coincidence,
    meaning the muon physically traversed the NaI crystal. NaI scintillation only
    fires on direct ionisation, so d_perp ≡ 0 by event selection. The previous
    guard on reconstructed d_perp was wrong — it dropped events where LineFit
    was inaccurate, not events where the muon missed the crystal.
    """
    def __init__(self, ctx):
        super().__init__(ctx)

    def Configure(self): pass

    def Physics(self, frame):
        if DM_T_KEY not in frame or LF_KEY not in frame:
            self.PushFrame(frame)
            return

        dm_id  = frame[DM_ID_KEY].value
        dm_pos = DMICE_POS_IC[dm_id]

        lf = frame[LF_KEY]
        lf_dir = np.array([lf.dir.x, lf.dir.y, lf.dir.z])
        lf_pos = np.array([lf.pos.x, lf.pos.y, lf.pos.z])

        # NOTE: All events in this dataset are DM-Ice coincidences — the muon
        # physically passed through the NaI crystal (d_perp ≈ 0 by selection).
        # The scintillator only fires on direct ionisation, so a DM-Ice hit is
        # definitive proof the muon crossed the crystal. We therefore apply the
        # spatial anchor unconditionally; the old d_perp > D_MAX guard was
        # incorrectly cutting events where the *reconstructed* (not true) track
        # had large d_perp due to poor LineFit at low energies.

        # Re-reference vertex to dm_pos using IC-derived timing (no DM-Ice time used)
        r = np.array(dm_pos) - lf_pos
        s = float(np.dot(r, lf_dir))
        t0_at_dm = lf.time + s / C_M_NS  # IC-only t0 at dm_pos

        pp = dataclasses.I3Particle()
        pp.dir    = lf.dir                           # same direction as LineFit
        pp.pos    = dataclasses.I3Position(*dm_pos)  # vertex at DM-Ice
        pp.time   = float(t0_at_dm)
        pp.fit_status = dataclasses.I3Particle.OK
        frame[PIV_SPAT_KEY] = pp
        self.PushFrame(frame)


class DMCombinedFitModule(icetray.I3Module):
    """
    Combined IC Pandel + DM-Ice Gaussian likelihood minimisation.

    Fixes track vertex at dm_pos (d_perp≈0 by definition), then runs
    scipy Nelder-Mead over (zenith, azimuth, t0) to maximise:
        log L_IC_pandel(zen, azi, t0) + log L_DM_gaussian(t0)

    Produces two outputs:
      MPEFit_DM  — started from MPEFit(std) seed
      SPEFit_DM  — started from SPEFit(std) seed
    """
    def __init__(self, ctx):
        super().__init__(ctx)

    def Configure(self): pass

    def Physics(self, frame):
        if DM_T_KEY not in frame:
            self.PushFrame(frame)
            return

        dm_id  = frame[DM_ID_KEY].value
        dm_pos = tuple(DMICE_POS_IC[dm_id])
        dm_t_corrected = frame[DM_T_KEY].value - MU_NS

        # Build IC pulse list (charge-weighted; use first pulse per DOM hit)
        if IC_PULSES not in frame:
            self.PushFrame(frame)
            return
        pulses = []
        for omk, plist in frame[IC_PULSES]:
            if omk not in geo_obj.omgeo:
                continue
            pos = geo_obj.omgeo[omk].position
            for p in plist:
                pulses.append((pos.x, pos.y, pos.z, float(p.time), float(p.charge)))
        if len(pulses) < 4:
            self.PushFrame(frame)
            return

        def run_fit(seed_key, out_key):
            if seed_key not in frame:
                return
            seed = frame[seed_key]
            seed_pos = (seed.pos.x, seed.pos.y, seed.pos.z)

            # Keep vertex at seed position — moving it to dm_pos breaks IC Pandel
            # geometry for all IceCube DOMs. The DM-Ice term computes t_geo_DM
            # from the track direction and seed vertex internally.
            x0 = [seed.dir.zenith, seed.dir.azimuth, seed.time]
            result = scipy_minimize(
                _neg_combined_ll,
                x0,
                args=(seed_pos, dm_pos, dm_t_corrected, pulses, SIGMA_NS),
                method='Nelder-Mead',
                options={'maxiter': 800, 'xatol': 1e-4, 'fatol': 0.5},
            )
            zen_opt, azi_opt, t0_opt = result.x

            pp = dataclasses.I3Particle()
            pp.dir    = dataclasses.I3Direction(float(zen_opt % math.pi),
                                                 float(azi_opt % (2*math.pi)))
            pp.pos    = seed.pos   # vertex stays at seed position
            pp.time   = float(t0_opt)
            pp.fit_status = dataclasses.I3Particle.OK
            frame[out_key] = pp

        run_fit(MPE_STD, MPE_DM_KEY)
        run_fit(SPE_STD, SPE_DM_KEY)
        self.PushFrame(frame)


# ── IceTray modules ───────────────────────────────────────────────────────────

class NPZInjector(icetray.I3Module):
    """Inject NPZ events into IceTray as I3Frames."""
    def __init__(self, ctx):
        super().__init__(ctx)
        self.idx = 0

    def Configure(self): pass

    def Process(self):
        if self.idx == 0:
            gframe = icetray.I3Frame(icetray.I3Frame.Geometry)
            gframe["I3Geometry"] = geo_obj
            self.PushFrame(gframe)

        if self.idx >= N:
            self.RequestSuspension()
            return

        i = self.idx; self.idx += 1

        frame = icetray.I3Frame(icetray.I3Frame.Physics)

        hdr = dataclasses.I3EventHeader()
        hdr.run_id = 3000; hdr.event_id = i
        frame["I3EventHeader"] = hdr

        # MC truth (BLO convention: travel direction, downgoing dz<0)
        zen = float(d["zenith_rad"][i]); azi = float(d["azimuth_rad"][i])
        dx = math.sin(zen) * math.cos(azi)
        dy = math.sin(zen) * math.sin(azi)
        dz = math.cos(zen)    # <0 downgoing
        primary = dataclasses.I3Particle()
        primary.type       = dataclasses.I3Particle.MuMinus
        primary.shape      = dataclasses.I3Particle.InfiniteTrack
        primary.energy     = float(d["energy_GeV"][i]) * I3Units.GeV
        primary.dir        = dataclasses.I3Direction(dx, dy, dz)
        primary.fit_status = dataclasses.I3Particle.FitStatus.OK
        mc_tree = dataclasses.I3MCTree()
        mc_tree.add_primary(primary)
        frame["I3MCTree"]  = mc_tree
        frame["MCTruth"]   = primary

        # Pulses: all DOMs (IC + DM-Ice)
        dom_x_ev = np.array(_dom_x[i]); dom_y_ev = np.array(_dom_y[i])
        dom_z_ev = np.array(_dom_z[i]); dom_t_ev = np.array(_dom_t[i])
        dom_nh   = np.array(_dom_nhits[i])
        dom_str  = np.array(_dom_string[i], dtype=int)
        dom_sen  = np.array(_dom_sensor[i], dtype=int)

        all_pm = dataclasses.I3RecoPulseSeriesMap()
        ic_pm  = dataclasses.I3RecoPulseSeriesMap()
        dm_t   = None
        dm_id  = -1

        for j in range(len(dom_x_ev)):
            s_j = int(dom_str[j]); sen_j = int(dom_sen[j])
            omkey = icetray.OMKey(s_j, sen_j)
            ps    = dataclasses.I3RecoPulseSeries()
            p     = dataclasses.I3RecoPulse()
            p.time   = float(dom_t_ev[j])
            p.charge = float(dom_nh[j])
            ps.append(p)
            all_pm[omkey] = ps
            # DM-Ice DOMs: strings 87/88
            is_dm = False
            for det_id, (dm_s, dm_sen) in DMICE_OMKEYS.items():
                if s_j == dm_s and sen_j == dm_sen:
                    if dm_t is None or float(dom_t_ev[j]) < dm_t:
                        dm_t  = float(dom_t_ev[j])
                        dm_id = det_id
                    is_dm = True
                    break
            if not is_dm:
                ic_pm[omkey] = ps

        frame[ALL_PULSES] = all_pm
        frame[IC_PULSES]  = ic_pm

        # Use PPC DM-Ice hit if available; fall back to injected direct-ionization time
        if dm_t is None and "dm_t_injected_ns" in d:
            injected = float(d["dm_t_injected_ns"][i])
            if not math.isnan(injected):
                dm_t  = injected
                dm_id = int(d["target_det"][i]) if "target_det" in d else 0

        if dm_t is not None:
            frame[DM_T_KEY] = dataclasses.I3Double(dm_t)
            frame[DM_ID_KEY] = icetray.I3Int(dm_id)

        # Metadata
        if "target_det" in d:
            tgt_det = int(d["target_det"][i])
        elif "det_id" in d:
            tgt_det = 1 if "det2" in str(d["det_id"][i]) else 0
        else:
            tgt_det = 0
        frame["TargetDet"] = icetray.I3Int(tgt_det)
        if "bin_id" in d:
            frame["BinId"] = icetray.I3Int(int(d["bin_id"][i]))
        else:
            frame["BinId"] = icetray.I3Int(0)
        frame["NDoms"]     = icetray.I3Int(len(ic_pm))
        frame["NHits"]     = icetray.I3Int(int(np.sum(dom_nh[dom_str < 87])))

        self.PushFrame(frame)


class PivotLFModule(icetray.I3Module):
    """Compute Pivot LineFit anchored to DM-Ice hit time from IC-only pulses."""
    def __init__(self, ctx):
        super().__init__(ctx)

    def Configure(self): pass

    def Physics(self, frame):
        if LF_KEY not in frame or DM_T_KEY not in frame:
            self.PushFrame(frame); return

        lf    = frame[LF_KEY]
        dm_id = frame[DM_ID_KEY].value if DM_ID_KEY in frame else 0
        dm_pos = DMICE_POS_IC[dm_id]
        seed_dir = np.array([lf.dir.x, lf.dir.y, lf.dir.z])

        # Apply μ correction: remove timing model offset (NaI scintillation delay)
        dm_t_corrected = frame[DM_T_KEY].value - MU_NS

        # Collect IC-only pulses
        try:
            pulses = dataclasses.I3RecoPulseSeriesMap.from_frame(frame, IC_PULSES)
        except Exception:
            pulses = frame[IC_PULSES]
        geo = frame["I3Geometry"].omgeo
        xs, ys, zs, ts, ws = [], [], [], [], []
        for omk, plist in pulses:
            if omk not in geo:
                continue
            pos = geo[omk].position
            for p in plist:
                xs.append(pos.x); ys.append(pos.y); zs.append(pos.z)
                ts.append(p.time); ws.append(p.charge)

        if len(xs) < 4:
            self.PushFrame(frame); return

        xs = np.array(xs); ys = np.array(ys); zs = np.array(zs)
        ts = np.array(ts); ws = np.array(ws, dtype=float)

        result = pivot_linefit_ic(xs, ys, zs, ts, ws, dm_pos, dm_t_corrected, seed_dir)
        if result is None:
            self.PushFrame(frame); return
        piv_dir, t0_pivot = result

        pp = dataclasses.I3Particle()
        pp.dir        = dataclasses.I3Direction(float(piv_dir[0]), float(piv_dir[1]), float(piv_dir[2]))
        # Anchor vertex at DM-Ice position with DM-Ice-constrained t0
        pp.pos        = dataclasses.I3Position(float(dm_pos[0]), float(dm_pos[1]), float(dm_pos[2]))
        pp.time       = float(t0_pivot)
        pp.fit_status = dataclasses.I3Particle.FitStatus.OK
        frame[PIV_LF_KEY] = pp
        self.PushFrame(frame)


class MPEPivotModule(icetray.I3Module):
    """Pivot LineFit anchored to DM-Ice, using MPEFit(std) direction as seed."""
    def __init__(self, ctx):
        super().__init__(ctx)

    def Configure(self): pass

    def Physics(self, frame):
        if MPE_STD not in frame or DM_T_KEY not in frame:
            self.PushFrame(frame); return

        mpe = frame[MPE_STD]
        if mpe.fit_status != dataclasses.I3Particle.FitStatus.OK:
            self.PushFrame(frame); return

        seed_dir = np.array([mpe.dir.x, mpe.dir.y, mpe.dir.z])
        dm_id    = frame[DM_ID_KEY].value if DM_ID_KEY in frame else 0
        dm_pos   = DMICE_POS_IC[dm_id]
        dm_t_corrected = frame[DM_T_KEY].value - MU_NS

        try:
            pulses = dataclasses.I3RecoPulseSeriesMap.from_frame(frame, IC_PULSES)
        except Exception:
            pulses = frame[IC_PULSES]
        geo = frame["I3Geometry"].omgeo
        xs, ys, zs, ts, ws = [], [], [], [], []
        for omk, plist in pulses:
            if omk not in geo:
                continue
            pos = geo[omk].position
            for p in plist:
                xs.append(pos.x); ys.append(pos.y); zs.append(pos.z)
                ts.append(p.time); ws.append(p.charge)

        if len(xs) < 4:
            self.PushFrame(frame); return

        xs = np.array(xs); ys = np.array(ys); zs = np.array(zs)
        ts = np.array(ts); ws = np.array(ws, dtype=float)

        result = pivot_linefit_ic(xs, ys, zs, ts, ws, dm_pos, dm_t_corrected, seed_dir)
        if result is None:
            self.PushFrame(frame); return
        piv_dir, t0_pivot = result

        pp = dataclasses.I3Particle()
        pp.dir        = dataclasses.I3Direction(float(piv_dir[0]), float(piv_dir[1]), float(piv_dir[2]))
        pp.pos        = dataclasses.I3Position(float(dm_pos[0]), float(dm_pos[1]), float(dm_pos[2]))
        pp.time       = float(t0_pivot)
        pp.fit_status = dataclasses.I3Particle.FitStatus.OK
        frame[PIV_MPE_KEY] = pp
        self.PushFrame(frame)


# ── Scoring and CSV extraction ────────────────────────────────────────────────

rows = []

class Scorer(icetray.I3Module):
    def __init__(self, ctx):
        super().__init__(ctx)

    def Configure(self): pass

    def Physics(self, frame):
        if "MCTruth" not in frame:
            self.PushFrame(frame); return

        truth  = frame["MCTruth"]
        mc_dir = np.array([truth.dir.x, truth.dir.y, truth.dir.z])
        e_GeV  = truth.energy / I3Units.GeV

        def get_dir(key):
            if key not in frame:
                return None, None, None
            p = frame[key]
            if p.fit_status != dataclasses.I3Particle.FitStatus.OK:
                return None, None, None
            return (np.array([p.dir.x, p.dir.y, p.dir.z]),
                    np.array([p.pos.x, p.pos.y, p.pos.z]),
                    p.time)

        def ang(key):
            dr, _, _ = get_dir(key)
            return ang_err_deg(mc_dir, dr) if dr is not None else float("nan")

        def score(key, d_perp_override=None):
            dr, dp, dt0 = get_dir(key)
            if dr is None or DM_T_KEY not in frame:
                return float("nan"), float("nan"), float("nan")
            dm_t  = frame[DM_T_KEY].value
            dm_id = frame[DM_ID_KEY].value if DM_ID_KEY in frame else 0
            dm_pos = DMICE_POS_IC[dm_id]
            return dm_log_l(dp, dr, dt0, dm_t, dm_pos, d_perp_override)

        ll_mc,  dt_mc,  _ = score("MCTruth", d_perp_override=0.0)
        ll_lf,  dt_lf,  dp_lf  = score(LF_KEY)
        ll_plf, dt_plf, _       = score(PIV_LF_KEY, d_perp_override=0.0)
        ll_mpe,  _, dp_mpe = score(MPE_STD)
        ll_mpiv, _, _      = score(MPE_PIV,  d_perp_override=0.0)
        ll_mpiv2,_, _      = score(MPE_PIV2, d_perp_override=0.0)
        ll_spe,  _, _      = score(SPE_STD)
        ll_spiv, _, _      = score(SPE_PIV, d_perp_override=0.0)
        ll_mpe_dm,  _, _ = score(MPE_DM_KEY,   d_perp_override=0.0)
        ll_spe_dm,  _, _ = score(SPE_DM_KEY,   d_perp_override=0.0)
        ll_mpe_sp,  _, _ = score(MPE_SPAT_KEY, d_perp_override=0.0)
        ll_spe_sp,  _, _ = score(SPE_SPAT_KEY, d_perp_override=0.0)

        rows.append(dict(
            event_id           = frame["I3EventHeader"].event_id,
            bin_id             = frame["BinId"].value if "BinId" in frame else -1,
            mc_energy_GeV      = e_GeV,
            n_doms_ic          = frame["NDoms"].value if "NDoms" in frame else 0,
            has_dm_hit         = int(DM_T_KEY in frame),
            # Angular errors — progression story
            lf_ang_err         = ang(LF_KEY),           # LineFit: no DM-Ice
            piv_spat_ang_err   = ang(PIV_SPAT_KEY),     # LineFit: DM spatial only
            piv_lf_ang_err     = ang(PIV_LF_KEY),       # LineFit: DM spatial+time
            mpe_std_ang_err    = ang(MPE_STD),           # MPEFit: no DM-Ice
            mpe_spat_ang_err   = ang(MPE_SPAT_KEY),     # MPEFit: DM spatial seed
            mpe_piv_ang_err    = ang(MPE_PIV),           # MPEFit: DM spatial+time seed
            mpe_dm_ang_err     = ang(MPE_DM_KEY),        # MPEFit: DM combined LL (no re-seed)
            mpe_piv2_ang_err   = ang(MPE_PIV2),          # MPEFit: MPE-pivot seed
            spe_std_ang_err    = ang(SPE_STD),           # SPEFit: no DM-Ice
            spe_spat_ang_err   = ang(SPE_SPAT_KEY),     # SPEFit: DM spatial seed
            spe_piv_ang_err    = ang(SPE_PIV),           # SPEFit: DM spatial+time seed
            spe_dm_ang_err     = ang(SPE_DM_KEY),        # SPEFit: DM combined LL (no re-seed)
            # DM-Ice log L
            ll_mc              = ll_mc,
            ll_lf              = ll_lf,
            ll_plf             = ll_plf,
            ll_mpe_std         = ll_mpe,
            ll_mpe_piv         = ll_mpiv,
            ll_mpe_piv2        = ll_mpiv2,
            ll_mpe_dm          = ll_mpe_dm,
            ll_mpe_spat        = ll_mpe_sp,
            ll_spe_std         = ll_spe,
            ll_spe_piv         = ll_spiv,
            ll_spe_dm          = ll_spe_dm,
            ll_spe_spat        = ll_spe_sp,
            dp_lf_m            = dp_lf,
            dp_mpe_m           = dp_mpe,
        ))
        self.PushFrame(frame)


# ── Build tray ────────────────────────────────────────────────────────────────

tray = I3Tray()

tray.Add(NPZInjector, "NPZInjector")

# LineFit on IC-only pulses
tray.Add("I3LineFit",
    Name            = LF_KEY,
    InputRecoPulses = IC_PULSES,
    AmpWeightPower  = 1.0,
)

# Pivot LineFit (seeded from LineFit direction)
tray.Add(PivotLFModule, "PivotLF")

# SplineMPE — standard seed
tray.Add(spline_reco.SplineMPE, SMPE_STD,
    fitname           = SMPE_STD,
    PulsesName        = IC_PULSES,
    TrackSeedList     = [LF_KEY],
    BareMuTimingSpline    = SPLINE_PROB,
    BareMuAmplitudeSpline = SPLINE_AMP,
    configuration     = "default",
    If = lambda f: LF_KEY in f and len(f[IC_PULSES]) >= 4,
)

# SplineMPE — pivot seed
tray.Add(spline_reco.SplineMPE, SMPE_PIV,
    fitname           = SMPE_PIV,
    PulsesName        = IC_PULSES,
    TrackSeedList     = [PIV_LF_KEY],
    BareMuTimingSpline    = SPLINE_PROB,
    BareMuAmplitudeSpline = SPLINE_AMP,
    configuration     = "default",
    If = lambda f: PIV_LF_KEY in f and len(f[IC_PULSES]) >= 4,
)

# MPEFit — standard seed
tray.Add(icecube.lilliput.segments.I3SinglePandelFitter,
    fitname = MPE_STD,
    domllh  = "MPE",
    pulses  = IC_PULSES,
    seeds   = [LF_KEY],
    If      = lambda f: LF_KEY in f and len(f[IC_PULSES]) >= 4,
)

# MPEFit — LineFit-pivot seed
tray.Add(icecube.lilliput.segments.I3SinglePandelFitter,
    fitname = MPE_PIV,
    domllh  = "MPE",
    pulses  = IC_PULSES,
    seeds   = [PIV_LF_KEY],
    If      = lambda f: PIV_LF_KEY in f and len(f[IC_PULSES]) >= 4,
)

# Pivot anchored from MPEFit(std) direction — runs after MPEFit std
tray.Add(MPEPivotModule, "MPEPivot")

# MPEFit — MPEFit-pivot seed (second-stage pivot)
tray.Add(icecube.lilliput.segments.I3SinglePandelFitter,
    fitname = MPE_PIV2,
    domllh  = "MPE",
    pulses  = IC_PULSES,
    seeds   = [PIV_MPE_KEY],
    If      = lambda f: PIV_MPE_KEY in f and len(f[IC_PULSES]) >= 4,
)

# SPEFit — standard seed
tray.Add(icecube.lilliput.segments.I3SinglePandelFitter,
    fitname = SPE_STD,
    domllh  = "SPE1st",
    pulses  = IC_PULSES,
    seeds   = [LF_KEY],
    If      = lambda f: LF_KEY in f and len(f[IC_PULSES]) >= 4,
)

# SPEFit — pivot seed
tray.Add(icecube.lilliput.segments.I3SinglePandelFitter,
    fitname = SPE_PIV,
    domllh  = "SPE1st",
    pulses  = IC_PULSES,
    seeds   = [PIV_LF_KEY],
    If      = lambda f: PIV_LF_KEY in f and len(f[IC_PULSES]) >= 4,
)

# Spatial-only pivot (position anchor, no DM-Ice timing)
tray.Add(SpatialPivotModule, "SpatialPivot")

# MPEFit — spatial-only pivot seed (position info only, no timing)
tray.Add(icecube.lilliput.segments.I3SinglePandelFitter,
    fitname = MPE_SPAT_KEY,
    domllh  = "MPE",
    pulses  = IC_PULSES,
    seeds   = [PIV_SPAT_KEY],
    If      = lambda f: PIV_SPAT_KEY in f and len(f[IC_PULSES]) >= 4,
)

# SPEFit — spatial-only pivot seed
tray.Add(icecube.lilliput.segments.I3SinglePandelFitter,
    fitname = SPE_SPAT_KEY,
    domllh  = "SPE1st",
    pulses  = IC_PULSES,
    seeds   = [PIV_SPAT_KEY],
    If      = lambda f: PIV_SPAT_KEY in f and len(f[IC_PULSES]) >= 4,
)

# Combined IC Pandel + DM-Ice Gaussian — MPEFit and SPEFit seeds
tray.Add(DMCombinedFitModule, "DMCombinedFit")

tray.Add(Scorer, "Scorer")

tray.Execute()
tray.Finish()

# ── Write CSV ─────────────────────────────────────────────────────────────────

import pandas as pd
df = pd.DataFrame(rows)
csv_path = OUT_BASE + ".csv"
df.to_csv(csv_path, index=False)
print(f"\nDone: {len(rows)} events")
print(f"CSV: {csv_path}")

# ── Summary ───────────────────────────────────────────────────────────────────

METHODS = [
    ("lf_ang_err",       "LineFit"),
    ("piv_lf_ang_err",   "Pivot LineFit"),
    ("mpe_std_ang_err",  "MPEFit (std seed)"),
    ("mpe_piv_ang_err",  "MPEFit (LF-pivot seed)"),
    ("mpe_piv2_ang_err", "MPEFit (MPE-pivot seed)"),
    ("spe_std_ang_err",  "SPEFit (std seed)"),
    ("spe_piv_ang_err",  "SPEFit (piv seed)"),
]

has_dm = df[df.has_dm_hit == 1]
print(f"\nEvents with DM-Ice hit: {len(has_dm)}/{len(df)}")
print("\n── Overall ──")
for col, label in METHODS:
    vals = df[col].dropna()
    print(f"  {label:28s}: med={vals.median():.2f}°  n={len(vals)}")

print("\n── Per energy bin ──")
for bin_id in sorted(df.bin_id.dropna().unique()):
    sub = df[df.bin_id == bin_id]
    e_med = sub.mc_energy_GeV.median()
    print(f"\n  Bin {int(bin_id)} (median E={e_med:.0f} GeV, n={len(sub)}):")
    for col, label in METHODS:
        vals = sub[col].dropna()
        if len(vals):
            print(f"    {label:28s}: med={vals.median():.2f}°  n={len(vals)}")

# ── Plot ──────────────────────────────────────────────────────────────────────

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

abins = np.linspace(0, 15, 61)
style = [
    ("lf_ang_err",      "LineFit",          "steelblue",   "-"),
    ("piv_lf_ang_err",  "Pivot LF",         "darkorange",  "-"),
    ("mpe_std_ang_err", "MPEFit (std)",     "forestgreen", "--"),
    ("mpe_piv_ang_err", "MPEFit (pivot)",   "red",         "--"),
    ("spe_std_ang_err", "SPEFit (std)",     "purple",      ":"),
    ("spe_piv_ang_err", "SPEFit (pivot)",   "magenta",     ":"),
]

# ── Figure 1: overall + MPEFit vs SPEFit scatter ─────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(13, 5))

ax = axes[0]
for col, label, color, ls in style:
    vals = df[col].dropna()
    if len(vals) > 5:
        ax.hist(vals[vals <= 15], bins=abins, histtype="step", lw=2, ls=ls,
                label=f"{label} ({vals.median():.2f}°)", color=color, density=True)
ax.set_xlabel("Angular error (°)"); ax.set_ylabel("Normalised events / bin")
ax.set_title(f"All {len(df)} events"); ax.legend(fontsize=7); ax.grid(True, alpha=0.3)

ax = axes[1]
df["mpe_best"] = df[["mpe_std_ang_err", "mpe_piv_ang_err"]].min(axis=1)
df["spe_best"] = df[["spe_std_ang_err", "spe_piv_ang_err"]].min(axis=1)
sub = df[df.mpe_best.notna() & df.spe_best.notna()]
if len(sub):
    ax.scatter(sub.mpe_best.clip(0, 15), sub.spe_best.clip(0, 15),
               alpha=0.5, s=12, color="steelblue")
    ax.plot([0, 15], [0, 15], "k--", lw=1, alpha=0.5)
    n_spe = (sub.spe_best < sub.mpe_best).sum()
    ax.set_title(f"Best SPEFit vs best MPEFit\nSPE wins: {n_spe}/{len(sub)} events")
ax.set_xlabel("MPEFit best ang err (°)"); ax.set_ylabel("SPEFit best ang err (°)")
ax.grid(True, alpha=0.3)

fig.suptitle(
    f"Reco comparison — BLO sim ({len(df)} events)  "
    f"μ={MU_NS:+.0f}ns σ={SIGMA_NS:.0f}ns d⊥<{D_MAX:.0f}m",
    fontsize=10)
plt.tight_layout()
fig.savefig(OUT_BASE + ".png", dpi=150, bbox_inches="tight")
plt.close()
print(f"Plot: {OUT_BASE}.png")

# ── Figure 2: per energy bin ──────────────────────────────────────────────────
bins_list = sorted(df.bin_id.dropna().unique())
fig2, axes2 = plt.subplots(1, len(bins_list), figsize=(4*len(bins_list), 4), sharey=True)
if len(bins_list) == 1:
    axes2 = [axes2]
for ax2, bid in zip(axes2, bins_list):
    sub = df[df.bin_id == bid]
    e_med = sub.mc_energy_GeV.median()
    for col, label, color, ls in style:
        vals = sub[col].dropna()
        if len(vals) > 3:
            ax2.hist(vals[vals <= 15], bins=abins, histtype="step", lw=1.5, ls=ls,
                     color=color, density=True,
                     label=f"{label.split('(')[0].strip()} {vals.median():.2f}°")
    ax2.set_title(f"Bin {int(bid)}\n{e_med:.0f} GeV (n={len(sub)})", fontsize=9)
    ax2.set_xlabel("Angular error (°)"); ax2.legend(fontsize=6); ax2.grid(True, alpha=0.3)
axes2[0].set_ylabel("Normalised events / bin")
fig2.suptitle("Per-energy-bin angular error — MPEFit vs SPEFit vs pivot", fontsize=10)
plt.tight_layout()
fig2.savefig(OUT_BASE + "_per_bin.png", dpi=150, bbox_inches="tight")
plt.close()
print(f"Per-bin plot: {OUT_BASE}_per_bin.png")
