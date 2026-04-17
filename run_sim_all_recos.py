#!/usr/bin/env python3
"""
run_sim_all_recos.py

Runs LineFit, Pivot LineFit, MPEFit, and Pivot MPEFit on BLO 200-event
binned downgoing muon sim. Computes angular error vs MC truth for each.

Output CSV columns:
    mc_energy_GeV, n_doms, n_hits, bin_id,
    ic_lf_ang_err_deg, pivot_lf_ang_err_deg,
    mpe_ang_err_deg, pivot_mpe_ang_err_deg

Usage (on Cobalt with IceTray env):
    /cvmfs/.../env-shell.sh python3 run_sim_all_recos.py
"""

import os, sys, csv, math
import numpy as np
from scipy.optimize import minimize as scipy_minimize
from scipy.special import gammaln

_default_npz = "~/dmice_work/output/muons_binned_200ev_repacked.npz"
NPZ_FILE     = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else _default_npz)
GEO_FILE     = os.path.expanduser("~/dmice/BlueLightOrchestra.jl/resources/geofiles/icecube_with_dmice.geo")
_default_csv = "~/dmice_work/output/comparison/sim_all_recos.csv"
OUT_CSV      = os.path.expanduser(sys.argv[2] if len(sys.argv) > 2 else _default_csv)
# Optional: override detector for all events (e.g. det_center whose target_id=1 by legacy bug)
DET_OVERRIDE = sys.argv[3] if len(sys.argv) > 3 else None  # e.g. "det_center"
# Optional: chunk slicing for Condor parallelism  argv[4]=chunk_id  argv[5]=n_chunks
CHUNK_ID     = int(sys.argv[4]) if len(sys.argv) > 4 else None
N_CHUNKS     = int(sys.argv[5]) if len(sys.argv) > 5 else None

# DM-Ice positions in IceCube coordinates [m] (z_BLO + 1948.07)
DMICE_POS_IC = {
    "det1":       np.array([ 31.25,  -72.93, -511.05]),
    "det2":       np.array([-334.80, -424.50, -511.26]),
    "det_center": np.array([  0.0,     0.0,     0.0  ]),  # IceCube geometric center
}
Z_OFFSET = 1948.07   # BLO z → IceCube z
C_M_NS   = 0.2998
N_ICE    = 1.3195
THETA_C  = math.acos(1.0 / N_ICE)
MU_NS    = 280.0     # NaI mean scintillation delay [ns]
SIGMA_NS =  81.0     # NaI scintillation jitter [ns]

# Approximate uniform ice Pandel parameters (SPICEMie bulk)
_PANDEL_LA = 98.0    # absorption length [m]
_PANDEL_LS = 30.0    # scattering length [m]
_JITTER_NS = 15.0    # DOM timing jitter [ns]

os.makedirs(os.path.dirname(OUT_CSV), exist_ok=True)

from icecube import icetray, dataclasses, dataio, simclasses
from icecube import linefit, lilliput, gulliver, gulliver_modules
import icecube.lilliput.segments
from icecube.spline_reco import SplineMPE
from icecube.icetray import I3Units, I3Tray

SPLINE_TIMING_BARE  = "/cvmfs/icecube.opensciencegrid.org/data/photon-tables/splines/InfBareMu_mie_prob_z20a10_V2.fits"
SPLINE_AMP_BARE     = "/cvmfs/icecube.opensciencegrid.org/data/photon-tables/splines/InfBareMu_mie_abs_z20a10_V2.fits"
SPLINE_TIMING_STOCH = "/cvmfs/icecube.opensciencegrid.org/data/photon-tables/splines/InfHighEStoch_mie_prob_z20a10.fits"
SPLINE_AMP_STOCH    = "/cvmfs/icecube.opensciencegrid.org/data/photon-tables/splines/InfHighEStoch_mie_abs_z20a10.fits"

# ── Geometry ─────────────────────────────────────────────────────────────────

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

# ── Build I3Geometry frame ────────────────────────────────────────────────────

geo_obj = dataclasses.I3Geometry()
for (s, dom), (px, py, pz) in geo_doms.items():
    omkey = icetray.OMKey(s, dom)
    omgeo = dataclasses.I3OMGeo()
    omgeo.position = dataclasses.I3Position(px, py, pz)
    omgeo.omtype   = dataclasses.I3OMGeo.IceCube
    geo_obj.omgeo[omkey] = omgeo

# ── Load npz ─────────────────────────────────────────────────────────────────

d = np.load(NPZ_FILE, allow_pickle=True)
N_TOTAL = len(d["energy_GeV"])
if CHUNK_ID is not None and N_CHUNKS is not None:
    chunk_size = math.ceil(N_TOTAL / N_CHUNKS)
    EV_START   = CHUNK_ID * chunk_size
    EV_END     = min(EV_START + chunk_size, N_TOTAL)
else:
    EV_START, EV_END = 0, N_TOTAL
N = EV_END - EV_START
print(f"Loaded {N_TOTAL} events from {NPZ_FILE}  (processing {EV_START}:{EV_END})")

def load_ragged(key):
    if f"{key}_flat" in d:
        flat    = d[f"{key}_flat"]
        offsets = d[f"{key}_offsets"]
        return [flat[offsets[i]:offsets[i+1]] for i in range(EV_START, EV_END)]
    return d[key][EV_START:EV_END]

_dom_x      = load_ragged("dom_x")
_dom_y      = load_ragged("dom_y")
_dom_z      = load_ragged("dom_z")
_dom_t      = load_ragged("dom_t")
_dom_nhits  = load_ragged("dom_nhits")
_dom_string = load_ragged("dom_string")
_dom_sensor = load_ragged("dom_sensor")

# Scalar arrays — sliced to chunk range, indexed with local i
_zenith      = d["zenith_rad"][EV_START:EV_END]
_azimuth     = d["azimuth_rad"][EV_START:EV_END]
_energy      = d["energy_GeV"][EV_START:EV_END]
_bin_id      = d["bin_id"][EV_START:EV_END]
_target_det  = d["target_det"][EV_START:EV_END]
_n_hits      = d["n_hits"][EV_START:EV_END]
_n_doms      = d["n_doms"][EV_START:EV_END]
_dm_t        = d["dm_t_injected_ns"][EV_START:EV_END] if "dm_t_injected_ns" in d else None

# ── Pivot LineFit (Python) ────────────────────────────────────────────────────

def _wm(vals, ws):
    W = sum(ws)
    return sum(v * w for v, w in zip(vals, ws)) / W if W else 0.0

def pivot_linefit_ic(xs, ys, zs, ts, ws, dm_pos_ic, seed_dir, dm_t=None):
    cx, cy, cz = _wm(xs, ws), _wm(ys, ws), _wm(zs, ws)
    tb = _wm(ts, ws)
    d_proj = ((dm_pos_ic[0]-cx)*seed_dir[0] + (dm_pos_ic[1]-cy)*seed_dir[1]
              + (dm_pos_ic[2]-cz)*seed_dir[2])
    # Use actual μ-corrected DM-Ice hit time if available;
    # fall back to LineFit extrapolation for events without a DM-Ice hit.
    t_dm = dm_t if dm_t is not None else tb + d_proj / C_M_NS
    dts  = [t - t_dm for t in ts]
    drxs = [x - dm_pos_ic[0] for x in xs]
    drys = [y - dm_pos_ic[1] for y in ys]
    drzs = [z - dm_pos_ic[2] for z in zs]
    den  = sum(w * dt * dt for w, dt in zip(ws, dts))
    if not den:
        return None
    vx = sum(w * dt * dr for w, dt, dr in zip(ws, dts, drxs)) / den
    vy = sum(w * dt * dr for w, dt, dr in zip(ws, dts, drys)) / den
    vz = sum(w * dt * dr for w, dt, dr in zip(ws, dts, drzs)) / den
    spd = math.sqrt(vx*vx + vy*vy + vz*vz)
    return (vx/spd, vy/spd, vz/spd) if spd else None

def ang_err_deg(truth, reco):
    dot = max(-1.0, min(1.0, truth[0]*reco[0] + truth[1]*reco[1] + truth[2]*reco[2]))
    return math.degrees(math.acos(abs(dot)))

# ── Frame injection module ────────────────────────────────────────────────────

class NPZInjector(icetray.I3Module):
    def __init__(self, context):
        super().__init__(context)
        self.idx = 0

    def Configure(self): pass

    def Process(self):
        # Push geometry frame once
        if self.idx == 0:
            gframe = icetray.I3Frame(icetray.I3Frame.Geometry)
            gframe["I3Geometry"] = geo_obj
            self.PushFrame(gframe)

        if self.idx >= N:
            self.RequestSuspension()
            return

        i = self.idx          # local index into sliced arrays
        self.idx += 1

        frame = icetray.I3Frame(icetray.I3Frame.Physics)

        hdr = dataclasses.I3EventHeader()
        hdr.run_id   = 3000
        hdr.event_id = EV_START + i   # global event id for merge dedup
        frame["I3EventHeader"] = hdr

        # MC truth direction (BLO frame: zenith=arccos(dz), dz<0=downgoing)
        zen = float(_zenith[i])
        azi = float(_azimuth[i])
        dx_mc = math.sin(zen) * math.cos(azi)
        dy_mc = math.sin(zen) * math.sin(azi)
        dz_mc = math.cos(zen)   # < 0 for downgoing

        primary = dataclasses.I3Particle()
        primary.type          = dataclasses.I3Particle.MuMinus
        primary.shape         = dataclasses.I3Particle.InfiniteTrack
        primary.energy        = float(_energy[i]) * I3Units.GeV
        primary.dir           = dataclasses.I3Direction(dx_mc, dy_mc, dz_mc)
        primary.fit_status    = dataclasses.I3Particle.FitStatus.OK
        mc_tree = dataclasses.I3MCTree()
        mc_tree.add_primary(primary)
        frame["I3MCTree"] = mc_tree
        frame["MCTruth"]  = primary   # convenience shortcut

        # Pulses: one per DOM, charge=nhits, t=first-hit time
        # DOM positions in IceCube coords (z_BLO + Z_OFFSET)
        dom_x_ev = np.array(_dom_x[i])
        dom_y_ev = np.array(_dom_y[i])
        dom_z_ev = np.array(_dom_z[i])  # BLO coords
        dom_t_ev = np.array(_dom_t[i])
        dom_nh   = np.array(_dom_nhits[i])
        dom_str  = np.array(_dom_string[i], dtype=int)
        dom_sen  = np.array(_dom_sensor[i], dtype=int)

        pulse_map = dataclasses.I3RecoPulseSeriesMap()
        for j in range(len(dom_x_ev)):
            omkey = icetray.OMKey(int(dom_str[j]), int(dom_sen[j]))
            ps    = dataclasses.I3RecoPulseSeries()
            p     = dataclasses.I3RecoPulse()
            p.time   = float(dom_t_ev[j])
            p.charge = float(dom_nh[j])
            ps.append(p)
            pulse_map[omkey] = ps
        frame["InIcePulses"] = pulse_map

        # Store metadata for extraction
        frame["BinId"]     = icetray.I3Int(int(_bin_id[i]))
        frame["TargetDet"] = icetray.I3Int(int(_target_det[i]))
        frame["NHits"]     = icetray.I3Int(int(_n_hits[i]))
        frame["NDoms"]     = icetray.I3Int(int(_n_doms[i]))

        # DM-Ice hit time (analytically injected, NaI scintillation model)
        if _dm_t is not None:
            dm_t = float(_dm_t[i])
            if math.isfinite(dm_t):
                frame[DM_T_KEY] = dataclasses.I3Double(dm_t)

        self.PushFrame(frame)

# ── Pivot LineFit + Pivot MPEFit modules ─────────────────────────────────────

PIVOT_LF_KEY     = "PivotLineFit"
SPE_KEY          = "SPEFit"
PIVOT_SPE_KEY    = "PivotSPEFit"
PIVOT_MPE_KEY    = "PivotMPEFit"
COMBINED_SPE_KEY = "CombinedSPEFit"   # IC Pandel + DM-Ice Gaussian, pivot seed
COMBINED_MPE_KEY = "CombinedMPEFit"   # IC Pandel + DM-Ice Gaussian, pivot seed
ITER_MPE_KEY      = "IterMPE"           # Iterative Pandel SPE (3 iterations)
ITER_PIVOT_LF_KEY = "IterPivotLineFit"  # Pivot LineFit seeded from IterMPE direction
ITER_PIVOT_MPE_KEY= "IterPivotMPEFit"   # MPEFit seeded from IterPivotLineFit
SPLINE_STD_KEY    = "SplineMPE_std"     # SplineMPE seeded from LineFit
SPLINE_PIV_KEY    = "SplineMPE_piv"     # SplineMPE seeded from PivotLineFit
SPLINE_ITER_KEY   = "SplineMPE_iter"    # SplineMPE seeded from IterPivotLineFit
DM_T_KEY         = "DMIce_t"          # DM-Ice corrected hit time in frame

def compute_pivot_lf(frame):
    if "LineFit" not in frame or "InIcePulses" not in frame:
        return

    lf     = frame["LineFit"]
    lf_dir = (lf.dir.x, lf.dir.y, lf.dir.z)

    tgt_id = frame["TargetDet"].value if "TargetDet" in frame else 0
    if DET_OVERRIDE:
        dm_pos = DMICE_POS_IC[DET_OVERRIDE]
    elif tgt_id == 0:
        dm_pos = DMICE_POS_IC["det1"]
    elif tgt_id == 2:
        dm_pos = DMICE_POS_IC["det_center"]
    else:
        dm_pos = DMICE_POS_IC["det2"]

    try:
        pulses = dataclasses.I3RecoPulseSeriesMap.from_frame(frame, "InIcePulses")
    except Exception:
        return

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
        return

    # Use LineFit-extrapolated crossing time at DM-Ice (dm_t=None).
    # Ablation tests showed that subtracting MU_NS (NaI scintillation mean) here
    # makes Pivot LineFit worse in simulation because BLO/PPC generates Cherenkov
    # photon arrival times for DM-Ice, not NaI scintillation times.
    piv = pivot_linefit_ic(xs, ys, zs, ts, ws, dm_pos, lf_dir, dm_t=None)
    if piv is None:
        return

    # Anchor MPEFit seed vertex at DM-Ice position.
    # Use raw DM-Ice hit time (no μ correction) to set seed t0.
    # Step-2 ablation: this is neutral on Pivot LineFit, slightly improves MPEFit seed.
    dm_t_raw = frame[DM_T_KEY].value if DM_T_KEY in frame else None
    lf_particle = frame["LineFit"]
    pp = dataclasses.I3Particle()
    pp.dir = dataclasses.I3Direction(piv[0], piv[1], piv[2])
    if dm_t_raw is not None:
        s = ((dm_pos[0] - lf_particle.pos.x) * piv[0] +
             (dm_pos[1] - lf_particle.pos.y) * piv[1] +
             (dm_pos[2] - lf_particle.pos.z) * piv[2])
        pp.pos  = dataclasses.I3Position(dm_pos[0], dm_pos[1], dm_pos[2])
        pp.time = dm_t_raw - s / C_M_NS
    else:
        pp.pos  = lf_particle.pos
        pp.time = lf_particle.time
    pp.fit_status = dataclasses.I3Particle.FitStatus.OK
    frame[PIVOT_LF_KEY] = pp

# ── IterMPE-seeded pivot ──────────────────────────────────────────────────────

def compute_iter_pivot_lf(frame):
    """Pivot LineFit anchored using IterMPE direction instead of LineFit."""
    if ITER_MPE_KEY not in frame or DM_T_KEY not in frame:
        return
    iter_p = frame[ITER_MPE_KEY]
    if iter_p.fit_status != dataclasses.I3Particle.FitStatus.OK:
        return

    dm_id  = frame["TargetDet"].value
    if DET_OVERRIDE:
        dm_key = DET_OVERRIDE
    elif dm_id == 0:
        dm_key = "det1"
    elif dm_id == 2:
        dm_key = "det_center"
    else:
        dm_key = "det2"
    dm_pos = DMICE_POS_IC[dm_key]

    dm_t_raw = frame[DM_T_KEY].value
    piv = (iter_p.dir.x, iter_p.dir.y, iter_p.dir.z)

    pp = dataclasses.I3Particle()
    pp.dir = dataclasses.I3Direction(*piv)

    if DM_T_KEY in frame:
        s = ((dm_pos[0] - iter_p.pos.x) * piv[0] +
             (dm_pos[1] - iter_p.pos.y) * piv[1] +
             (dm_pos[2] - iter_p.pos.z) * piv[2])
        pp.pos  = dataclasses.I3Position(dm_pos[0], dm_pos[1], dm_pos[2])
        pp.time = dm_t_raw - s / C_M_NS
    else:
        pp.pos  = iter_p.pos
        pp.time = iter_p.time
    pp.fit_status = dataclasses.I3Particle.FitStatus.OK
    frame[ITER_PIVOT_LF_KEY] = pp


# ── Extraction ────────────────────────────────────────────────────────────────

rows = []

def extract(frame):
    if "MCTruth" not in frame:
        return

    truth = frame["MCTruth"]
    mc_dir = (truth.dir.x, truth.dir.y, truth.dir.z)
    e_GeV  = truth.energy / I3Units.GeV

    def ang(key):
        try:
            p = frame[key]
            if p.fit_status != dataclasses.I3Particle.FitStatus.OK:
                return float("nan")
            return ang_err_deg(mc_dir, (p.dir.x, p.dir.y, p.dir.z))
        except Exception:
            return float("nan")

    rows.append(dict(
        mc_energy_GeV            = e_GeV,
        n_doms                   = frame["NDoms"].value if "NDoms" in frame else 0,
        n_hits                   = frame["NHits"].value if "NHits" in frame else 0,
        bin_id                   = frame["BinId"].value if "BinId" in frame else -1,
        ic_lf_ang_err_deg        = ang("LineFit"),
        pivot_lf_ang_err_deg     = ang(PIVOT_LF_KEY),
        spe_ang_err_deg          = ang(SPE_KEY),
        pivot_spe_ang_err_deg    = ang(PIVOT_SPE_KEY),
        mpe_ang_err_deg          = ang("MPEFit"),
        pivot_mpe_ang_err_deg    = ang(PIVOT_MPE_KEY),
        combined_spe_ang_err_deg  = ang(COMBINED_SPE_KEY),
        combined_mpe_ang_err_deg  = ang(COMBINED_MPE_KEY),
        iter_mpe_ang_err_deg       = ang(ITER_MPE_KEY),
        iter_pivot_lf_ang_err_deg  = ang(ITER_PIVOT_LF_KEY),
        iter_pivot_mpe_ang_err_deg = ang(ITER_PIVOT_MPE_KEY),
        spline_std_ang_err_deg     = ang(SPLINE_STD_KEY),
        spline_piv_ang_err_deg     = ang(SPLINE_PIV_KEY),
        spline_iter_ang_err_deg    = ang(SPLINE_ITER_KEY),
    ))

# ── Combined IC Pandel + DM-Ice Gaussian likelihood ──────────────────────────

def _pandel_log_spe(t_res, d_perp):
    d = max(d_perp, 1.0)
    alpha = d / _PANDEL_LS
    beta  = C_M_NS / _PANDEL_LA
    if t_res < 0:
        return -0.5 * (t_res / _JITTER_NS)**2 - math.log(_JITTER_NS * math.sqrt(2*math.pi))
    return (alpha * math.log(beta) + (alpha - 1) * math.log(t_res + 1e-6)
            - beta * t_res - gammaln(alpha))


def _ic_log_l(zen, azi, t0, vertex, ic_pulse_list):
    sin_z, cos_z = math.sin(zen), math.cos(zen)
    sin_a, cos_a = math.sin(azi), math.cos(azi)
    dx = sin_z * cos_a; dy = sin_z * sin_a; dz = -cos_z
    vx, vy, vz = vertex
    ll = 0.0
    for (px, py, pz, t_hit, charge) in ic_pulse_list:
        rx, ry, rz = px - vx, py - vy, pz - vz
        s = rx*dx + ry*dy + rz*dz
        d_perp = math.sqrt(max(0.0, rx*rx + ry*ry + rz*rz - s*s))
        t_geo  = t0 + s / C_M_NS
        if d_perp > 0.01:
            t_geo += d_perp / (C_M_NS * math.sin(THETA_C))
        ll += charge * _pandel_log_spe(t_hit - t_geo, d_perp)
    return ll


def _neg_combined_ll(params, vertex, dm_pos, dm_t_corrected, ic_pulse_list):
    zen, azi, t0 = params
    if not (0.0 < zen < math.pi):
        return 1e9
    sin_z, cos_z = math.sin(zen), math.cos(zen)
    sin_a, cos_a = math.sin(azi), math.cos(azi)
    dx = sin_z * cos_a; dy = sin_z * sin_a; dz = -cos_z
    vx, vy, vz = vertex
    rx = dm_pos[0]-vx; ry = dm_pos[1]-vy; rz = dm_pos[2]-vz
    s_dm     = rx*dx + ry*dy + rz*dz
    t_geo_dm = t0 + s_dm / C_M_NS
    ll_ic = _ic_log_l(zen, azi, t0, vertex, ic_pulse_list)
    ll_dm = -0.5 * ((dm_t_corrected - t_geo_dm) / SIGMA_NS)**2
    return -(ll_ic + ll_dm)


class DMCombinedFitModule(icetray.I3Module):
    """
    Combined IC Pandel + DM-Ice Gaussian likelihood.
    Seeded from PivotLineFit; optimises (zenith, azimuth, t0) with scipy Nelder-Mead.
    DM-Ice constraint: Gaussian(mu=0, sigma=SIGMA_NS) on residual after subtracting MU_NS.
    Produces CombinedSPEFit and CombinedMPEFit (identical optimizer, different labels
    reflecting that SPE/MPE IceTray fits are replaced by the combined objective).
    """
    def __init__(self, ctx):
        super().__init__(ctx)

    def Configure(self): pass

    def Physics(self, frame):
        if DM_T_KEY not in frame or PIVOT_LF_KEY not in frame or "InIcePulses" not in frame:
            self.PushFrame(frame)
            return

        dm_id  = frame["TargetDet"].value
        if DET_OVERRIDE:
            dm_key = DET_OVERRIDE
        elif dm_id == 0:
            dm_key = "det1"
        elif dm_id == 2:
            dm_key = "det_center"
        else:
            dm_key = "det2"
        dm_pos = tuple(DMICE_POS_IC[dm_key])
        dm_t_corrected = frame[DM_T_KEY].value - MU_NS

        pulses = []
        for omk, plist in frame["InIcePulses"]:
            if omk not in geo_obj.omgeo:
                continue
            pos = geo_obj.omgeo[omk].position
            for p in plist:
                pulses.append((pos.x, pos.y, pos.z, float(p.time), float(p.charge)))

        if len(pulses) < 4:
            self.PushFrame(frame)
            return

        seed = frame[PIVOT_LF_KEY]
        vertex = (seed.pos.x, seed.pos.y, seed.pos.z)

        # Compute t0 consistent with seed direction + DM-Ice hit time:
        # the track reaches dm_pos at dm_t_corrected, so t0 = dm_t_corrected - s/c
        # where s = (dm_pos - vertex) · d̂
        sin_z = math.sin(seed.dir.zenith); cos_z = math.cos(seed.dir.zenith)
        sin_a = math.sin(seed.dir.azimuth); cos_a = math.cos(seed.dir.azimuth)
        dx0 = sin_z*cos_a; dy0 = sin_z*sin_a; dz0 = -cos_z
        vx, vy, vz = vertex
        s_seed = ((dm_pos[0]-vx)*dx0 + (dm_pos[1]-vy)*dy0 + (dm_pos[2]-vz)*dz0)
        t0_init = dm_t_corrected - s_seed / C_M_NS

        x0 = [seed.dir.zenith, seed.dir.azimuth, t0_init]

        result = scipy_minimize(
            _neg_combined_ll, x0,
            args=(vertex, dm_pos, dm_t_corrected, pulses),
            method='Nelder-Mead',
            options={'maxiter': 800, 'xatol': 1e-4, 'fatol': 0.5},
        )
        zen_opt, azi_opt, t0_opt = result.x

        pp = dataclasses.I3Particle()
        pp.dir        = dataclasses.I3Direction(float(zen_opt % math.pi),
                                                float(azi_opt % (2*math.pi)))
        pp.pos        = seed.pos
        pp.time       = float(t0_opt)
        pp.fit_status = dataclasses.I3Particle.FitStatus.OK
        frame[COMBINED_SPE_KEY] = pp
        frame[COMBINED_MPE_KEY] = pp
        self.PushFrame(frame)


# ── Build tray ────────────────────────────────────────────────────────────────

tray = I3Tray()

tray.Add(NPZInjector)

tray.Add("I3LineFit",
    Name            = "LineFit",
    InputRecoPulses = "InIcePulses",
    AmpWeightPower  = 1.0,
)

tray.Add(compute_pivot_lf, Streams=[icetray.I3Frame.Physics])

tray.Add(icecube.lilliput.segments.I3SinglePandelFitter,
    fitname = SPE_KEY,
    domllh  = "SPE1st",
    pulses  = "InIcePulses",
    seeds   = [PIVOT_LF_KEY],
    If      = lambda f: PIVOT_LF_KEY in f,
)

tray.Add(icecube.lilliput.segments.I3SinglePandelFitter,
    fitname = PIVOT_SPE_KEY,
    domllh  = "SPE1st",
    pulses  = "InIcePulses",
    seeds   = [PIVOT_LF_KEY],
    If      = lambda f: PIVOT_LF_KEY in f,
)

tray.Add(icecube.lilliput.segments.I3SinglePandelFitter,
    fitname = "MPEFit",
    domllh  = "MPE",
    pulses  = "InIcePulses",
    seeds   = [PIVOT_LF_KEY],
    If      = lambda f: PIVOT_LF_KEY in f,
)

tray.Add(icecube.lilliput.segments.I3SinglePandelFitter,
    fitname = PIVOT_MPE_KEY,
    domllh  = "MPE",
    pulses  = "InIcePulses",
    seeds   = [PIVOT_LF_KEY],
    If      = lambda f: PIVOT_LF_KEY in f,
)

tray.Add(DMCombinedFitModule)

# IterativePandelFit: 2 iterations, SPE1st (more stable than MPE for iterative)
# then use result as seed for a final MPEFit pass
tray.Add(icecube.lilliput.segments.I3IterativePandelFitter,
    fitname      = ITER_MPE_KEY,
    domllh       = "SPE1st",
    pulses       = "InIcePulses",
    seeds        = ["LineFit"],
    n_iterations = 3,
    If           = lambda f: "LineFit" in f,
)

# Pivot LineFit seeded from IterMPE direction
tray.Add(compute_iter_pivot_lf, Streams=[icetray.I3Frame.Physics])

# MPEFit seeded from IterPivotLineFit
tray.Add(icecube.lilliput.segments.I3SinglePandelFitter,
    fitname = ITER_PIVOT_MPE_KEY,
    domllh  = "MPE",
    pulses  = "InIcePulses",
    seeds   = [ITER_PIVOT_LF_KEY],
    If      = lambda f: ITER_PIVOT_LF_KEY in f,
)

# SplineMPE — three seeds: standard LineFit, PivotLineFit, IterPivotLineFit
SPLINE_COMMON = dict(
    configuration        = "recommended",
    PulsesName           = "InIcePulses",
    EnergyEstimators     = ["MCTruth"],  # MC truth energy (sim only)
    BareMuTimingSpline   = SPLINE_TIMING_BARE,
    BareMuAmplitudeSpline= SPLINE_AMP_BARE,
    StochTimingSpline    = SPLINE_TIMING_STOCH,
    StochAmplitudeSpline = SPLINE_AMP_STOCH,
)

tray.Add(SplineMPE, fitname=SPLINE_STD_KEY,
    TrackSeedList=["LineFit"],
    If=lambda f: "LineFit" in f,
    **SPLINE_COMMON)

tray.Add(SplineMPE, fitname=SPLINE_PIV_KEY,
    TrackSeedList=[PIVOT_LF_KEY],
    If=lambda f: PIVOT_LF_KEY in f,
    **SPLINE_COMMON)

tray.Add(SplineMPE, fitname=SPLINE_ITER_KEY,
    TrackSeedList=[ITER_PIVOT_LF_KEY],
    If=lambda f: ITER_PIVOT_LF_KEY in f,
    **SPLINE_COMMON)

tray.Add(extract, Streams=[icetray.I3Frame.Physics])

tray.Execute()
tray.Finish()

# ── Write CSV ─────────────────────────────────────────────────────────────────

with open(OUT_CSV, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader()
    w.writerows(rows)

def valid(key): return sum(1 for r in rows if not math.isnan(r[key]))
def median_ang(key):
    vals = [r[key] for r in rows if not math.isnan(r[key])]
    return float(np.median(vals)) if vals else float("nan")

print(f"Done: {len(rows)} events")
print(f"  {'Fit':<22}  {'Valid':>5}  {'Median ang err':>14}")
for key, label in [
    ("ic_lf_ang_err_deg",        "LineFit"),
    ("pivot_lf_ang_err_deg",     "Pivot LineFit"),
    ("spe_ang_err_deg",          "SPEFit"),
    ("pivot_spe_ang_err_deg",    "Pivot SPEFit"),
    ("mpe_ang_err_deg",          "MPEFit"),
    ("pivot_mpe_ang_err_deg",    "Pivot MPEFit"),
    ("combined_spe_ang_err_deg",  "Combined SPEFit"),
    ("combined_mpe_ang_err_deg",  "Combined MPEFit"),
    ("iter_mpe_ang_err_deg",      "IterMPE"),
    ("iter_pivot_lf_ang_err_deg", "IterPivot LineFit"),
    ("iter_pivot_mpe_ang_err_deg","IterPivot MPEFit"),
    ("spline_std_ang_err_deg",    "SplineMPE (LineFit seed)"),
    ("spline_piv_ang_err_deg",    "SplineMPE (Pivot seed)"),
    ("spline_iter_ang_err_deg",   "SplineMPE (IterPivot seed)"),
]:
    print(f"  {label:<22}  {valid(key):>5}  {median_ang(key):>13.2f}°")
print(f"Saved: {OUT_CSV}")
