#!/usr/bin/env python3
"""
run_splinempe_nai_lambda.py

Hyperparameter scan over NaI likelihood weight λ for SplineMPE reconstruction.

Track model: the muon track passes through the DM-Ice detector position r_DM
at some time t0 (the transit time). We optimize over (zenith, azimuth, t0)
using a combined likelihood:

  log L_combined = log L_Pandel(zen, azi, t0; vertex=r_DM, IC pulses)
                 + λ · log G(T_meas; t0 + μ, σ)

where log G is the NaI timing Gaussian (μ, σ from the fitted timing model).
The NaI term is a soft constraint on t0 — at λ=0 the transit time is free
(determined only by IC pulses), at large λ it is pinned to T_meas - μ.

This fixes the bad-residual problem: Pandel evaluated from r_DM at the
NaI-measured transit time gives physically grounded residuals (no time
convention mismatch with SplineMPE).

Dataset-level loss functions stored per λ:
  loss_ang_mean : mean angular error
  loss_huber    : mean(huber(Δψ, δ=0.5°))
  loss_t0_res   : mean |t0_opt - dm_t_corrected| [ns] (timing residual)

Output: ~/dmice_work/output/nai_lambda_scan.csv  (one row per λ per bin)
        ~/dmice_work/output/nai_lambda_scan.png

Usage (Cobalt, IceTray env):
  env-shell.sh python3 ~/dmice/run_splinempe_nai_lambda.py
  env-shell.sh python3 ~/dmice/run_splinempe_nai_lambda.py --true-time
"""

import os, csv, math, argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("--npz", default=os.path.expanduser(
    "~/dmice_work/output/muons_binned_5000ev_repacked_injected.npz"))
parser.add_argument("--model", default=os.path.expanduser(
    "~/dmice_work/output/dmice_timing_model.npz"))
parser.add_argument("--out", default=os.path.expanduser(
    "~/dmice_work/output/nai_lambda_scan.csv"))
parser.add_argument("--lam-max", type=float, default=100.0)
parser.add_argument("--n-lam", type=int, default=11,
    help="Number of λ values (log-spaced between 0.05 and lam-max)")
parser.add_argument("--true-time", action="store_true",
    help="Use MC true DM-Ice transit time (upper bound on improvement)")
parser.add_argument("--huber-delta", type=float, default=0.5,
    help="δ for Huber loss [deg]")
args = parser.parse_args()

# ── Paths ─────────────────────────────────────────────────────────────────────

GEO_FILE    = os.path.expanduser(
    "~/dmice/BlueLightOrchestra.jl/resources/geofiles/icecube_with_dmice.geo")
SPLINE_PROB = ("/cvmfs/icecube.opensciencegrid.org/data/photon-tables/splines/"
               "InfBareMu_mie_prob_z20a10_V2.fits")
SPLINE_AMP  = ("/cvmfs/icecube.opensciencegrid.org/data/photon-tables/splines/"
               "InfBareMu_mie_abs_z20a10_V2.fits")

Z_OFFSET = 1948.07
C_M_NS   = 0.2998
N_ICE    = 1.3195
THETA_C  = math.acos(1.0 / N_ICE)

DMICE_OMKEYS = {0: (87, 1), 1: (88, 1)}
DMICE_POS_IC = {
    0: np.array([ 31.25,  -72.93, -511.05]),
    1: np.array([-334.80, -424.50, -511.26]),
}

# ── IceTray imports ───────────────────────────────────────────────────────────

from icecube import icetray, dataclasses, spline_reco, gulliver, linefit
from icecube import lilliput
import icecube.lilliput.segments
from icecube.icetray import I3Units, I3Tray
from scipy.optimize import minimize as scipy_minimize
from scipy.special import gammaln

# ── Timing model ──────────────────────────────────────────────────────────────

_m       = np.load(args.model, allow_pickle=True)
MU_NS    = float(_m["mu_ns"])
SIGMA_NS = float(_m["sigma_ns"])
print(f"Timing model: μ={MU_NS:+.1f} ns  σ={SIGMA_NS:.1f} ns")

# ── Geometry ──────────────────────────────────────────────────────────────────

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
geo_obj  = dataclasses.I3Geometry()
for (s, dom), (px, py, pz) in geo_doms.items():
    omkey = icetray.OMKey(s, dom)
    omgeo = dataclasses.I3OMGeo()
    omgeo.position = dataclasses.I3Position(px, py, pz)
    omgeo.omtype   = dataclasses.I3OMGeo.IceCube
    geo_obj.omgeo[omkey] = omgeo

# ── NPZ ───────────────────────────────────────────────────────────────────────

d = np.load(args.npz, allow_pickle=True)
N = len(d["energy_GeV"])

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

# ── IC likelihood: try spline tables, fall back to Pandel ────────────────────

# Pandel bulk ice parameters (SPICEMie bulk values)
_LA  = 98.0   # absorption length [m]
_LS  = 30.0   # scattering length [m]
_JIT = 15.0   # timing jitter [ns]

def _pandel_ic_log_l(zen, azi, t0, vertex, pulses_arr):
    """Vectorized Pandel log-likelihood. pulses_arr is (N,5) numpy array."""
    sz, cz = math.sin(zen), math.cos(zen)
    sa, ca = math.sin(azi), math.cos(azi)
    dx, dy, dz = sz*ca, sz*sa, -cz
    vx, vy, vz = vertex

    rx = pulses_arr[:,0] - vx
    ry = pulses_arr[:,1] - vy
    rz = pulses_arr[:,2] - vz
    s      = rx*dx + ry*dy + rz*dz
    d_perp = np.sqrt(np.maximum(0.0, rx**2 + ry**2 + rz**2 - s**2))
    d_perp = np.maximum(d_perp, 1.0)

    t_geo = t0 + s / C_M_NS + np.where(d_perp > 0.01,
                                        d_perp / (C_M_NS * math.sin(THETA_C)), 0.0)
    t_res = pulses_arr[:,3] - t_geo
    charge = pulses_arr[:,4]

    alpha = d_perp / _LS
    beta  = C_M_NS / _LA

    ll_pos = (alpha * np.log(beta)
              + (alpha - 1.0) * np.log(np.maximum(t_res, 1e-6))
              - beta * t_res - gammaln(alpha))
    ll_neg = (-0.5 * (t_res / _JIT)**2
              - math.log(_JIT * math.sqrt(2.0 * math.pi)))

    return float(np.sum(charge * np.where(t_res >= 0, ll_pos, ll_neg)))


class SplineLikelihoodWrapper:
    """
    Wraps I3SplineRecoLikelihood for use in scipy minimization.

    Requires I3PhotoSplineService for the photon tables.
    Falls back to Pandel automatically if anything fails.
    """

    def __init__(self):
        self.using_spline = False
        self._service = None
        self._current_frame = None
        try:
            from icecube import photonics_service
            # Load spline tables into a PhotoSplineService
            phot_svc = photonics_service.I3PhotoSplineService(
                SPLINE_AMP, SPLINE_PROB, "")

            svc = spline_reco.I3SplineRecoLikelihood()
            svc.Pulses = "InIcePulses"
            svc.PhotonicsService = phot_svc
            svc.SetGeometry(geo_obj)

            # Smoke-test with a real particle hypothesis on a dummy frame
            dummy = icetray.I3Frame(icetray.I3Frame.Physics)
            dummy["InIcePulses"] = dataclasses.I3RecoPulseSeriesMap()
            dummy["I3Geometry"] = geo_obj
            svc.SetEvent(dummy)
            p_test = dataclasses.I3Particle()
            p_test.dir  = dataclasses.I3Direction(1.0, 0.0)
            p_test.pos  = dataclasses.I3Position(0., 0., 0.)
            p_test.time = 0.0
            p_test.fit_status = dataclasses.I3Particle.FitStatus.OK
            _ = svc.GetLogLikelihood(gulliver.I3EventHypothesis(p_test))

            self._service = svc
            self.using_spline = True
            print("IC likelihood: I3SplineRecoLikelihood (spline tables)")
        except Exception as e:
            print(f"Spline likelihood not usable ({type(e).__name__}: {e}); "
                  "falling back to Pandel approximation")

    def set_event(self, frame):
        self._current_frame = frame
        if self.using_spline:
            try:
                self._service.SetEvent(frame)
            except Exception:
                self.using_spline = False

    def log_l(self, zen, azi, t0, vertex, pulses):
        if self.using_spline and self._current_frame is not None:
            try:
                p = dataclasses.I3Particle()
                p.dir  = dataclasses.I3Direction(float(zen), float(azi))
                p.pos  = dataclasses.I3Position(*[float(v) for v in vertex])
                p.time = float(t0)
                p.fit_status = dataclasses.I3Particle.FitStatus.OK
                hyp = gulliver.I3EventHypothesis(p)
                return float(self._service.GetLogLikelihood(hyp))
            except Exception:
                self.using_spline = False
        return _pandel_ic_log_l(zen, azi, t0, vertex,
                               pulses if isinstance(pulses, np.ndarray)
                               else np.array(pulses, dtype=np.float64))


spline_wrapper = SplineLikelihoodWrapper()
# Force Pandel: I3SplineRecoLikelihood.GetLogLikelihood() returns incorrect values
# (likely NaN or constant) when called from scipy Nelder-Mead outside Gulliver.
# Pandel is internally consistent across all λ values and fast enough for the scan.
spline_wrapper.using_spline = False
print("IC likelihood overridden to: Pandel (spline service unreliable in scipy context)")

# ── NaI term ──────────────────────────────────────────────────────────────────

def _nai_log_l(zen, azi, t0, vertex, dm_pos, dm_t_corrected, lam):
    sz, cz = math.sin(zen), math.cos(zen)
    sa, ca = math.sin(azi), math.cos(azi)
    dx = sz*ca; dy = sz*sa; dz = -cz
    vx, vy, vz = vertex
    s_dm = (dm_pos[0]-vx)*dx + (dm_pos[1]-vy)*dy + (dm_pos[2]-vz)*dz
    t_geo_dm = t0 + s_dm / C_M_NS
    return lam * (-0.5 * ((dm_t_corrected - t_geo_dm) / SIGMA_NS)**2)


def neg_combined(params, vertex, dm_pos, dm_t_corrected, pulses, lam):
    zen, azi, t0 = params
    if not (0.0 < zen < math.pi):
        return 1e9
    ll = spline_wrapper.log_l(zen, azi, t0, vertex, pulses)
    if lam > 0:
        ll += _nai_log_l(zen, azi, t0, vertex, dm_pos, dm_t_corrected, lam)
    return -ll

# ── Loss functions ────────────────────────────────────────────────────────────

def huber(x, delta):
    x = np.abs(x)
    return np.where(x <= delta, x**2, 2*delta*x - delta**2)


def compute_losses(ang_errs, t0_residuals_ns, delta=args.huber_delta):
    ang  = np.array(ang_errs)
    t0r  = np.array(t0_residuals_ns)
    return dict(
        n              = len(ang),
        ang_median     = float(np.median(ang)),
        ang_mean       = float(np.mean(ang)),
        t0_res_median  = float(np.median(np.abs(t0r))),
        t0_res_mean    = float(np.mean(np.abs(t0r))),
        loss_ang_mean  = float(np.mean(ang)),
        loss_huber     = float(np.mean(huber(ang, delta))),
        loss_t0_res    = float(np.mean(np.abs(t0r))),
    )

# ── Helper ────────────────────────────────────────────────────────────────────

def ang_err_deg(d1, d2):
    dot = max(-1.0, min(1.0, float(np.dot(np.asarray(d1), np.asarray(d2)))))
    return math.degrees(math.acos(abs(dot)))


# ── λ grid ────────────────────────────────────────────────────────────────────

lam_grid = [0.0] + list(np.logspace(-1, math.log10(args.lam_max), args.n_lam))
print(f"λ grid ({len(lam_grid)} values): {[f'{l:.3g}' for l in lam_grid]}")

# ── Frame store: holds the per-event I3Frame for spline service ───────────────

_frame_store = {}   # ev_id -> I3Frame (IC pulses only)
event_cache  = []   # list of dicts

# ── IceTray keys ─────────────────────────────────────────────────────────────

IC_PULSES = "InIcePulses"
LF_KEY    = "LineFit"
SMPE_KEY  = "SplineMPE_Std"
DM_T_KEY  = "DMIce_t"
DM_ID_KEY = "DMIce_id"

# ── IceTray modules ───────────────────────────────────────────────────────────

class NPZInjector(icetray.I3Module):
    def __init__(self, ctx):
        super().__init__(ctx)
        self.idx = 0

    def Configure(self): pass

    def Process(self):
        if self.idx == 0:
            gf = icetray.I3Frame(icetray.I3Frame.Geometry)
            gf["I3Geometry"] = geo_obj
            self.PushFrame(gf)

        if self.idx >= N:
            self.RequestSuspension(); return

        i = self.idx; self.idx += 1
        frame = icetray.I3Frame(icetray.I3Frame.Physics)

        hdr = dataclasses.I3EventHeader()
        hdr.run_id = 1; hdr.event_id = i
        frame["I3EventHeader"] = hdr

        zen = float(d["zenith_rad"][i]); azi = float(d["azimuth_rad"][i])
        dx = math.sin(zen)*math.cos(azi)
        dy = math.sin(zen)*math.sin(azi)
        dz = math.cos(zen)
        mc = dataclasses.I3Particle()
        mc.type  = dataclasses.I3Particle.MuMinus
        mc.dir   = dataclasses.I3Direction(dx, dy, dz)
        mc.energy = float(d["energy_GeV"][i]) * I3Units.GeV
        mc.fit_status = dataclasses.I3Particle.FitStatus.OK
        mct = dataclasses.I3MCTree(); mct.add_primary(mc)
        frame["I3MCTree"] = mct; frame["MCTruth"] = mc

        dom_x_ev = np.array(_dom_x[i]); dom_y_ev = np.array(_dom_y[i])
        dom_z_ev = np.array(_dom_z[i]); dom_t_ev = np.array(_dom_t[i])
        dom_nh   = np.array(_dom_nhits[i])
        dom_str  = np.array(_dom_string[i], dtype=int)
        dom_sen  = np.array(_dom_sensor[i], dtype=int)

        ic_pm = dataclasses.I3RecoPulseSeriesMap()
        dm_t  = None; dm_id = -1

        for j in range(len(dom_x_ev)):
            s_j = int(dom_str[j]); sen_j = int(dom_sen[j])
            omkey = icetray.OMKey(s_j, sen_j)
            ps = dataclasses.I3RecoPulseSeries()
            p  = dataclasses.I3RecoPulse()
            p.time   = float(dom_t_ev[j])
            p.charge = float(dom_nh[j])
            ps.append(p)
            is_dm = False
            for det_id, (dm_s, dm_sen) in DMICE_OMKEYS.items():
                if s_j == dm_s and sen_j == dm_sen:
                    if dm_t is None or float(dom_t_ev[j]) < dm_t:
                        dm_t = float(dom_t_ev[j]); dm_id = det_id
                    is_dm = True; break
            if not is_dm:
                ic_pm[omkey] = ps

        if dm_t is None and "dm_t_injected_ns" in d:
            inj = float(d["dm_t_injected_ns"][i])
            if not math.isnan(inj):
                dm_t  = inj
                dm_id = int(d["target_det"][i]) if "target_det" in d else 0

        frame[IC_PULSES] = ic_pm
        if dm_t is not None:
            frame[DM_T_KEY]  = dataclasses.I3Double(dm_t)
            frame[DM_ID_KEY] = icetray.I3Int(dm_id)

        frame["BinId"] = icetray.I3Int(int(d["bin_id"][i]) if "bin_id" in d else -1)
        self.PushFrame(frame)


class CacheExtractor(icetray.I3Module):
    """Cache SplineMPE result and event frame for the λ scan."""
    def __init__(self, ctx):
        super().__init__(ctx)

    def Configure(self): pass

    def Physics(self, frame):
        if "MCTruth" not in frame or DM_T_KEY not in frame:
            self.PushFrame(frame); return

        seed_key = SMPE_KEY if SMPE_KEY in frame else LF_KEY
        if seed_key not in frame:
            self.PushFrame(frame); return
        seed = frame[seed_key]
        if seed.fit_status != dataclasses.I3Particle.FitStatus.OK:
            self.PushFrame(frame); return
        if len(frame[IC_PULSES]) < 4:
            self.PushFrame(frame); return

        mc      = frame["MCTruth"]
        mc_dir  = (mc.dir.x, mc.dir.y, mc.dir.z)
        ev_id   = frame["I3EventHeader"].event_id
        dm_id   = frame[DM_ID_KEY].value
        dm_pos  = tuple(DMICE_POS_IC[dm_id])
        dm_t_obs = frame[DM_T_KEY].value

        # Build pulse array (N,5): x,y,z,t,charge — numpy for vectorized Pandel
        pulse_list = []
        geo = frame["I3Geometry"].omgeo
        for omk, plist in frame[IC_PULSES]:
            if omk not in geo: continue
            pos = geo[omk].position
            for p in plist:
                pulse_list.append((pos.x, pos.y, pos.z, float(p.time), float(p.charge)))
        pulses = np.array(pulse_list, dtype=np.float64) if pulse_list else np.zeros((0,5))

        # True transit time for --true-time mode
        true_dm_t = None
        if args.true_time and "dm_t_injected_ns" in d:
            inj = float(d["dm_t_injected_ns"][ev_id])
            if not math.isnan(inj):
                true_dm_t = inj

        smpe_dir = (seed.dir.x, seed.dir.y, seed.dir.z)

        # Store the frame for I3SplineRecoLikelihood.SetEvent().
        # SetEvent() requires I3Geometry to be in the frame (not just set via SetGeometry).
        if "I3Geometry" not in frame:
            frame["I3Geometry"] = geo_obj
        _frame_store[ev_id] = frame

        event_cache.append(dict(
            ev_id          = ev_id,
            bin_id         = frame["BinId"].value,
            energy_GeV     = mc.energy / I3Units.GeV,
            mc_dir         = mc_dir,
            seed_key       = seed_key,
            seed_zen       = seed.dir.zenith,
            seed_azi       = seed.dir.azimuth,
            seed_t0        = seed.time,
            seed_pos       = (seed.pos.x, seed.pos.y, seed.pos.z),
            smpe_dir       = smpe_dir,
            pulses         = pulses,
            dm_pos         = dm_pos,
            dm_t_corrected = dm_t_obs - MU_NS,
            true_dm_t      = true_dm_t,
            smpe_ang_err   = ang_err_deg(mc_dir, smpe_dir)
                             if seed_key == SMPE_KEY else float("nan"),
        ))
        self.PushFrame(frame)


# ── Run IceTray once ──────────────────────────────────────────────────────────

print(f"\nRunning SplineMPE (once) on {N} events...")
tray = I3Tray()
tray.Add(NPZInjector)
tray.Add("I3LineFit",
    Name            = LF_KEY,
    InputRecoPulses = IC_PULSES,
    AmpWeightPower  = 1.0,
)
tray.Add(spline_reco.SplineMPE, SMPE_KEY,
    fitname               = SMPE_KEY,
    PulsesName            = IC_PULSES,
    TrackSeedList         = [LF_KEY],
    BareMuTimingSpline    = SPLINE_PROB,
    BareMuAmplitudeSpline = SPLINE_AMP,
    configuration         = "default",
    If = lambda f: LF_KEY in f and len(f[IC_PULSES]) >= 4,
)
tray.Add(CacheExtractor)
tray.Execute()
tray.Finish()

print(f"Cached {len(event_cache)} events  (IC llh: pandel)")

# ── λ scan ────────────────────────────────────────────────────────────────────

all_rows = []   # one dict per (λ, bin_id)

_n_diag = 3   # print diagnostic for first N events of first λ

for lam in lam_grid:
    bin_results = {}   # bin_id -> lists of ang_err, t0_residual, energies
    _diag_count = 0

    for ev in event_cache:
        dm_pos  = np.array(ev["dm_pos"])
        dm_t_c  = (ev["true_dm_t"] if (args.true_time and ev["true_dm_t"] is not None)
                   else ev["dm_t_corrected"])

        if math.isnan(ev["smpe_ang_err"]):
            continue

        # Anchor track vertex at DM-Ice position — this makes Pandel residuals
        # physically meaningful (muon was at r_DM at time t0).
        # Optimize (zen, azi, t0); NaI soft-constrains t0 toward dm_t_c.
        vertex = dm_pos
        x0 = [ev["seed_zen"], ev["seed_azi"], dm_t_c]
        seed_val = neg_combined(x0, vertex, dm_pos, dm_t_c, ev["pulses"], lam)

        res = scipy_minimize(
            neg_combined, x0,
            args=(vertex, dm_pos, dm_t_c, ev["pulses"], lam),
            method="Nelder-Mead",
            options={"maxiter": 600, "xatol": 1e-3, "fatol": 0.5},
        )

        # Safety: if optimizer diverged or made things worse, fall back to SplineMPE seed
        fell_back = False
        if math.isnan(res.fun) or res.fun > seed_val + 1.0:
            zen_o, azi_o, t0_o = ev["seed_zen"], ev["seed_azi"], dm_t_c
            fell_back = True
        else:
            zen_o, azi_o, t0_o = res.x[0], res.x[1], res.x[2]

        sz, cz = math.sin(zen_o), math.cos(zen_o)
        sa, ca = math.sin(azi_o), math.cos(azi_o)
        reco_dir = (sz*ca, sz*sa, -cz)
        ang_e  = ang_err_deg(ev["mc_dir"], reco_dir)
        t0_res = t0_o - dm_t_c   # signed t0 residual [ns]

        if _diag_count < _n_diag:
            _diag_count += 1
            print(f"    [diag λ={lam:.2f} ev{ev['ev_id']}] "
                  f"smpe_err={ev['smpe_ang_err']:.2f}° → "
                  f"{'seed' if fell_back else 'opt'}={ang_e:.2f}°  "
                  f"t0_res={t0_res:+.1f}ns  "
                  f"seed_obj={seed_val:.1f}  res_obj={res.fun:.1f}  "
                  f"nfev={res.nfev}")

        b = ev["bin_id"]
        if b not in bin_results:
            bin_results[b] = {"ang": [], "t0r": [], "energies": []}
        bin_results[b]["ang"].append(ang_e)
        bin_results[b]["t0r"].append(t0_res)
        bin_results[b]["energies"].append(ev["energy_GeV"])

    # Compute losses per bin and overall
    all_ang = []; all_t0r = []
    for b, br in sorted(bin_results.items()):
        losses = compute_losses(br["ang"], br["t0r"])
        losses["lam"]    = lam
        losses["bin_id"] = b
        losses["med_energy_GeV"] = float(np.median(br["energies"]))
        losses["true_time"] = int(args.true_time)
        all_rows.append(losses)
        all_ang.extend(br["ang"]); all_t0r.extend(br["t0r"])

    # Overall row (bin_id = -1)
    if all_ang:
        losses = compute_losses(all_ang, all_t0r)
        losses["lam"]    = lam
        losses["bin_id"] = -1
        losses["med_energy_GeV"] = float("nan")
        losses["true_time"] = int(args.true_time)
        all_rows.append(losses)

    # Print summary
    if all_ang:
        print(f"  λ={lam:6.2f}  ang_med={np.median(all_ang):.3f}°  "
              f"t0_res_med={np.median(np.abs(all_t0r)):.1f}ns  "
              f"huber={np.mean(huber(np.array(all_ang), args.huber_delta)):.4f}  "
              f"n={len(all_ang)}")

# ── Save CSV ──────────────────────────────────────────────────────────────────

os.makedirs(os.path.dirname(args.out), exist_ok=True)
fieldnames = ["lam", "bin_id", "med_energy_GeV", "true_time", "n",
              "ang_median", "ang_mean", "t0_res_median", "t0_res_mean",
              "loss_ang_mean", "loss_huber", "loss_t0_res"]
with open(args.out, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
    w.writeheader(); w.writerows(all_rows)
print(f"\nSaved: {args.out}  ({len(all_rows)} rows)")

# ── Plot ──────────────────────────────────────────────────────────────────────

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

overall = [r for r in all_rows if r["bin_id"] == -1]
lams_o  = [r["lam"]            for r in overall]
ameds   = [r["ang_median"]     for r in overall]
amns    = [r["ang_mean"]       for r in overall]
t0meds  = [r["t0_res_median"]  for r in overall]
l_huber = [r["loss_huber"]     for r in overall]
l_ang   = [r["loss_ang_mean"]  for r in overall]

fig, axes = plt.subplots(3, 1, figsize=(8, 11), sharex=True)

ax = axes[0]
ax.plot(lams_o, ameds, "o-",  color="steelblue",  label="median Δψ")
ax.plot(lams_o, amns,  "s--", color="steelblue", alpha=0.6, label="mean Δψ")
ax.axhline(ameds[0], color="steelblue", lw=0.8, ls=":", alpha=0.5)
ax.set_ylabel("Angular error [deg]")
ax.legend(fontsize=9); ax.grid(True, alpha=0.3)
ax.set_title(f"SplineMPE+NaI (vertex=r_DM) — λ hyperparameter scan"
             + (" [true time]" if args.true_time else ""))

ax = axes[1]
ax.plot(lams_o, t0meds, "o-", color="darkorange", label="median |t0 residual| [ns]")
ax.axhline(0, color="darkorange", lw=0.8, ls=":", alpha=0.5)
ax.set_ylabel("|t0 − dm_t_corrected| [ns]")
ax.legend(fontsize=9); ax.grid(True, alpha=0.3)

ax = axes[2]
ax.plot(lams_o, l_ang,   "o-",  color="tomato",       label="loss: mean Δψ")
ax.plot(lams_o, l_huber, "s--", color="mediumpurple",  label=f"loss: Huber(δ={args.huber_delta}°)")
best_ang   = lams_o[int(np.argmin(l_ang))]
best_huber = lams_o[int(np.argmin(l_huber))]
ax.axvline(best_ang,   color="tomato",       ls="--", lw=1.2,
           label=f"opt λ (mean Δψ) = {best_ang:.2f}")
ax.axvline(best_huber, color="mediumpurple", ls="--", lw=1.2,
           label=f"opt λ (Huber) = {best_huber:.2f}")
ax.set_ylabel("Loss"); ax.set_xlabel("λ (NaI t0 weight)")
ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

# Log x-axis if range spans > 1 decade
pos_lams = [l for l in lams_o if l > 0]
if pos_lams and max(pos_lams) / min(pos_lams) > 10:
    for a in axes:
        a.set_xscale("symlog", linthresh=0.1)

plt.tight_layout()
out_png = args.out.replace(".csv", ".png")
fig.savefig(out_png, dpi=150, bbox_inches="tight")
plt.close()
print(f"Plot: {out_png}")
print(f"\nOptimal λ (mean Δψ)   = {best_ang:.3f}")
print(f"Optimal λ (Huber loss) = {best_huber:.3f}")
