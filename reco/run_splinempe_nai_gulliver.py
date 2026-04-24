#!/usr/bin/env python3
"""
run_splinempe_nai_gulliver.py

SplineMPE + NaI combined likelihood fit using Gulliver — the correct approach.

Implements:
    log L_combined = log L_SplineMPE(track) + λ · log G(T_meas; t_pred + μ, σ)

The NaI term is added as a Python subclass of I3EventLogLikelihood.
Gulliver's I3SimpleFitter optimises the combined likelihood using the
actual spline photon tables (not Pandel), so the IC term is correct.

SplineMPE runs once per event (standard fit). For each λ in the grid, a
second I3SimpleFitter pass refines the direction starting from the
SplineMPE seed with the combined likelihood.

Output:
    ~/dmice_work/output/nai_gulliver_scan.csv
    ~/dmice_work/output/nai_gulliver_scan.png

Usage (Cobalt):
    /cvmfs/.../env-shell.sh python3 -u ~/dmice/run_splinempe_nai_gulliver.py
    /cvmfs/.../env-shell.sh python3 -u ~/dmice/run_splinempe_nai_gulliver.py --true-time
"""

import os, csv, math, argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("--npz", default=os.path.expanduser(
    "~/dmice_work/output/muons_binned_5000ev_repacked_injected.npz"))
parser.add_argument("--model", default=os.path.expanduser(
    "~/dmice_work/output/dmice_timing_model.npz"))
parser.add_argument("--out", default=os.path.expanduser(
    "~/dmice_work/output/nai_gulliver_scan.csv"))
parser.add_argument("--true-time", action="store_true",
    help="Use MC true DM-Ice transit time (upper bound)")
parser.add_argument("--huber-delta", type=float, default=0.5)
args = parser.parse_args()

# ── Paths ─────────────────────────────────────────────────────────────────────

GEO_FILE    = os.path.expanduser(
    "~/dmice/BlueLightOrchestra.jl/resources/geofiles/icecube_with_dmice.geo")
SPLINE_PROB = ("/cvmfs/icecube.opensciencegrid.org/data/photon-tables/splines/"
               "InfBareMu_mie_prob_z20a10_V2.fits")
SPLINE_AMP  = ("/cvmfs/icecube.opensciencegrid.org/data/photon-tables/splines/"
               "InfBareMu_mie_abs_z20a10_V2.fits")

Z_OFFSET     = 1948.07
C_M_NS       = 0.2998

DMICE_OMKEYS = {0: (87, 1), 1: (88, 1)}
DMICE_POS_IC = {
    0: np.array([ 31.25,  -72.93, -511.05]),
    1: np.array([-334.80, -424.50, -511.26]),
}

IC_PULSES = "InIcePulses"
LF_KEY    = "LineFit"
SMPE_KEY  = "SplineMPE_Std"
DM_T_KEY  = "DMIce_t"
DM_ID_KEY = "DMIce_id"
DMTT_KEY  = "DMIce_t_true"   # true MC transit time (--true-time mode)

# λ grid: 0 = SplineMPE baseline, rest = NaI-weighted refinements
LAM_GRID = [0.0, 0.5, 1.0, 2.0, 5.0, 10.0, 25.0, 50.0, 100.0]

# ── IceTray ───────────────────────────────────────────────────────────────────

from icecube import (icetray, dataclasses, spline_reco, gulliver,
                     linefit, lilliput, photonics_service)
from icecube.icetray import I3Units, I3Tray

# ── Timing model ──────────────────────────────────────────────────────────────

_m       = np.load(args.model, allow_pickle=True)
MU_NS    = float(_m["mu_ns"])
SIGMA_NS = float(_m["sigma_ns"])
print(f"Timing model: μ={MU_NS:+.1f} ns  σ={SIGMA_NS:.1f} ns")

# ── Geometry (IceTray object + lookup dict) ───────────────────────────────────

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
    omgeo.position  = dataclasses.I3Position(px, py, pz)
    omgeo.omtype    = dataclasses.I3OMGeo.IceCube
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

# ── Shared spline photon service ──────────────────────────────────────────────

print("Loading spline tables...")
photon_svc = photonics_service.I3PhotoSplineService(SPLINE_AMP, SPLINE_PROB, "")
print("Spline tables loaded.")

# ── NaI Combined Likelihood ───────────────────────────────────────────────────

class NaICombinedLikelihood(gulliver.I3EventLogLikelihood):
    """
    I3EventLogLikelihood subclass that adds a NaI Gaussian penalty to the
    spline IC likelihood.

        log L = log L_SplineMPE(track) + λ · log G(T_meas; t_pred + μ, σ)

    The spline service is shared across λ instances (safe for sequential
    Gulliver execution). The NaI timing and DM-Ice position are read from
    the frame in SetEvent.
    """

    def __init__(self, name, spline_llh, lam, mu_ns, sigma_ns, true_time=False):
        super().__init__()
        self._name          = name
        self._spline_llh    = spline_llh
        self._lam           = lam
        self._mu_ns         = mu_ns
        self._sigma_ns      = sigma_ns
        self._true_time     = true_time
        self._dm_pos        = None   # set per-event in SetEvent
        self._dm_t_corr     = None   # dm_t_obs - mu_ns
        self._multiplicity  = 0

    def SetGeometry(self, geo):
        self._spline_llh.SetGeometry(geo)

    def SetEvent(self, frame):
        self._spline_llh.SetEvent(frame)
        self._multiplicity = self._spline_llh.GetMultiplicity()

        # DM-Ice position and corrected time
        if DM_T_KEY in frame and DM_ID_KEY in frame:
            dm_id = frame[DM_ID_KEY].value
            self._dm_pos = DMICE_POS_IC[dm_id]

            if self._true_time and DMTT_KEY in frame:
                self._dm_t_corr = frame[DMTT_KEY].value
            else:
                dm_t_obs = frame[DM_T_KEY].value
                self._dm_t_corr = dm_t_obs - self._mu_ns
        else:
            self._dm_pos    = None
            self._dm_t_corr = None

    def GetLogLikelihood(self, hyp):
        ll_ic = self._spline_llh.GetLogLikelihood(hyp)

        if self._lam == 0.0 or self._dm_t_corr is None or self._dm_pos is None:
            return ll_ic

        # NaI Gaussian term: penalise mismatch between predicted and observed
        # transit time at the DM-Ice position.
        p   = hyp.particle
        zen = p.dir.zenith
        azi = p.dir.azimuth
        dx  = math.sin(zen) * math.cos(azi)
        dy  = math.sin(zen) * math.sin(azi)
        dz  = -math.cos(zen)
        vx, vy, vz = p.pos.x, p.pos.y, p.pos.z
        s_dm = ((self._dm_pos[0] - vx) * dx
                + (self._dm_pos[1] - vy) * dy
                + (self._dm_pos[2] - vz) * dz)
        t_geo_dm = p.time + s_dm / C_M_NS
        ll_nai = self._lam * (
            -0.5 * ((self._dm_t_corr - t_geo_dm) / self._sigma_ns) ** 2)

        return ll_ic + ll_nai

    def GetMultiplicity(self):
        return self._multiplicity

    def GetName(self):
        return self._name

    def HasGradient(self):
        return False

# ── Create one spline likelihood (shared) and per-λ combined likelihoods ──────

# One I3SplineRecoLikelihood (uses the shared photon service)
spline_llh = spline_reco.I3SplineRecoLikelihood()
spline_llh.Pulses            = IC_PULSES
spline_llh.PhotonicsService  = photon_svc
spline_llh.LlhChoice         = "MPE"   # consistent with SplineMPE default
# SetGeometry will be called by Gulliver when the G frame passes through

nai_llhs = {}   # lam -> NaICombinedLikelihood instance
for lam in LAM_GRID:
    tag  = f"{lam:.4g}".replace(".", "p")
    name = f"NaILLH_{tag}"
    nai_llhs[lam] = NaICombinedLikelihood(
        name=name, spline_llh=spline_llh,
        lam=lam, mu_ns=MU_NS, sigma_ns=SIGMA_NS,
        true_time=args.true_time)

print(f"λ grid ({len(LAM_GRID)} values): {LAM_GRID}")

# ── Loss functions ────────────────────────────────────────────────────────────

def huber(x, delta):
    x = np.abs(x)
    return np.where(x <= delta, x**2, 2*delta*x - delta**2)

def ang_err_deg(d1, d2):
    dot = max(-1.0, min(1.0, float(np.dot(np.asarray(d1), np.asarray(d2)))))
    return math.degrees(math.acos(abs(dot)))

def d_perp_to_point(track_pos, track_dir, point):
    r  = np.asarray(point) - np.asarray(track_pos)
    dh = np.asarray(track_dir, dtype=float)
    dh = dh / np.linalg.norm(dh)
    return float(np.linalg.norm(r - np.dot(r, dh) * dh))

def compute_losses(ang_errs, d_perps, delta=args.huber_delta):
    ang = np.array(ang_errs)
    dp  = np.array(d_perps)
    return dict(
        n             = len(ang),
        ang_median    = float(np.median(ang)),
        ang_mean      = float(np.mean(ang)),
        dp_median     = float(np.median(dp)),
        dp_mean       = float(np.mean(dp)),
        loss_ang_mean = float(np.mean(ang)),
        loss_dp_mean  = float(np.mean(dp)),
        loss_combined = float(np.mean(ang) + np.mean(dp) / 100.0),
        loss_huber    = float(np.mean(huber(ang, delta))),
    )

# ── Scorer module ─────────────────────────────────────────────────────────────

# {lam -> {bin_id -> {ang: [], dp: [], energies: []}}}
results = {lam: {} for lam in LAM_GRID}

class Scorer(icetray.I3Module):
    def __init__(self, ctx):
        super().__init__(ctx)

    def Configure(self): pass

    def Physics(self, frame):
        if "MCTruth" not in frame or DM_T_KEY not in frame:
            self.PushFrame(frame); return

        mc     = frame["MCTruth"]
        mc_dir = (mc.dir.x, mc.dir.y, mc.dir.z)
        dm_id  = frame[DM_ID_KEY].value
        dm_pos = DMICE_POS_IC[dm_id]
        ev_id  = frame["I3EventHeader"].event_id
        ene    = mc.energy / I3Units.GeV
        bin_id = frame["BinId"].value

        for lam in LAM_GRID:
            tag     = f"{lam:.4g}".replace(".", "p")
            fit_key = f"NaiFit_{tag}" if lam > 0.0 else SMPE_KEY

            if fit_key not in frame:
                continue
            fit = frame[fit_key]
            if fit.fit_status != dataclasses.I3Particle.FitStatus.OK:
                continue

            reco_dir = (fit.dir.x, fit.dir.y, fit.dir.z)
            ang_e = ang_err_deg(mc_dir, reco_dir)
            dp    = d_perp_to_point(
                (fit.pos.x, fit.pos.y, fit.pos.z), reco_dir, dm_pos)

            r = results[lam]
            if bin_id not in r:
                r[bin_id] = {"ang": [], "dp": [], "energies": []}
            r[bin_id]["ang"].append(ang_e)
            r[bin_id]["dp"].append(dp)
            r[bin_id]["energies"].append(ene)

        self.PushFrame(frame)

# ── NPZInjector ───────────────────────────────────────────────────────────────

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
        dx  = math.sin(zen)*math.cos(azi)
        dy  = math.sin(zen)*math.sin(azi)
        dz  = math.cos(zen)
        mc  = dataclasses.I3Particle()
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
            s_j   = int(dom_str[j]); sen_j = int(dom_sen[j])
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
                        dm_t  = float(dom_t_ev[j]); dm_id = det_id
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

            # True transit time for --true-time mode
            if args.true_time and "dm_t_injected_ns" in d:
                inj = float(d["dm_t_injected_ns"][i])
                if not math.isnan(inj):
                    frame[DMTT_KEY] = dataclasses.I3Double(inj)

        frame["BinId"] = icetray.I3Int(int(d["bin_id"][i]) if "bin_id" in d else -1)
        self.PushFrame(frame)

# ── Build tray ────────────────────────────────────────────────────────────────

print(f"\nBuilding tray for {N} events...")

tray = I3Tray()

# Minimizer: Minuit2 simplex — fast for direction refinement
minimizer = lilliput.I3GulliverMinuit2(
    "NaIMinimizer", tolerance=0.001, max_iterations=1000)
tray.context["NaIMinimizer"] = minimizer

# Parametrization: Zen + Azi + T free; vertex position fixed.
# StepX/Y/Z=0 keeps the vertex pinned at the SplineMPE position.
tray.Add("I3SimpleParametrizationFactory", "NaIParam",
    StepX           = 0.0,
    StepY           = 0.0,
    StepZ           = 0.0,
    StepZenith      = 0.1 * I3Units.radian,
    StepAzimuth     = 0.2 * I3Units.radian,
    StepT           = 50.0 * I3Units.ns,
    BoundsZenith    = [0.0, math.pi],
    BoundsAzimuth   = [0.0, 2.0 * math.pi],
)

# Seed: SplineMPE result; TNone = trust SplineMPE's vertex time
tray.Add("I3BasicSeedServiceFactory", "NaISeed",
    InputReadout  = IC_PULSES,
    FirstGuesses  = [SMPE_KEY],
    TimeShiftType = "TNone",
)

# Modules
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

# For each λ > 0: register combined likelihood + run I3SimpleFitter
_has_dm = lambda f: (SMPE_KEY in f
                     and f[SMPE_KEY].fit_status == dataclasses.I3Particle.FitStatus.OK
                     and DM_T_KEY in f
                     and len(f[IC_PULSES]) >= 4)

for lam in LAM_GRID:
    if lam == 0.0:
        continue   # baseline = SplineMPE itself, no extra fitter needed

    llh      = nai_llhs[lam]
    llh_name = llh.GetName()
    tag      = f"{lam:.4g}".replace(".", "p")
    fit_name = f"NaiFit_{tag}"

    # Register likelihood in tray context so I3SimpleFitter can find it
    tray.context[llh_name] = llh

    tray.Add("I3SimpleFitter", fit_name,
        SeedService    = "NaISeed",
        Parametrization = "NaIParam",
        LogLikelihood  = llh_name,
        Minimizer      = "NaIMinimizer",
        OutputName     = fit_name,
        If             = _has_dm,
    )

tray.Add(Scorer)
tray.Execute()
tray.Finish()

print("IceTray done.\n")

# ── Aggregate results ─────────────────────────────────────────────────────────

all_rows = []
for lam in LAM_GRID:
    all_ang = []; all_dp = []
    for bin_id, res in sorted(results[lam].items()):
        if not res["ang"]:
            continue
        losses = compute_losses(res["ang"], res["dp"])
        losses["lam"]            = lam
        losses["bin_id"]         = bin_id
        losses["med_energy_GeV"] = float(np.median(res["energies"]))
        losses["true_time"]      = int(args.true_time)
        all_rows.append(losses)
        all_ang.extend(res["ang"]); all_dp.extend(res["dp"])

    if all_ang:
        losses = compute_losses(all_ang, all_dp)
        losses["lam"]            = lam
        losses["bin_id"]         = -1
        losses["med_energy_GeV"] = float("nan")
        losses["true_time"]      = int(args.true_time)
        all_rows.append(losses)
        delta = args.huber_delta
        print(f"  λ={lam:6.2f}  ang_med={np.median(all_ang):.3f}°  "
              f"d⊥_med={np.median(all_dp):.1f}m  "
              f"huber={np.mean(huber(np.array(all_ang), delta)):.4f}  "
              f"n={len(all_ang)}")

# ── Save CSV ──────────────────────────────────────────────────────────────────

os.makedirs(os.path.dirname(args.out), exist_ok=True)
fieldnames = ["lam", "bin_id", "med_energy_GeV", "true_time", "n",
              "ang_median", "ang_mean", "dp_median", "dp_mean",
              "loss_ang_mean", "loss_dp_mean", "loss_combined", "loss_huber"]
with open(args.out, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
    w.writeheader(); w.writerows(all_rows)
print(f"\nSaved: {args.out}  ({len(all_rows)} rows)")

# ── Plot ──────────────────────────────────────────────────────────────────────

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

overall  = [r for r in all_rows if r["bin_id"] == -1]
lams_o   = [r["lam"]           for r in overall]
ameds    = [r["ang_median"]     for r in overall]
amns     = [r["ang_mean"]       for r in overall]
dpmds    = [r["dp_median"]      for r in overall]
l_huber  = [r["loss_huber"]     for r in overall]
l_comb   = [r["loss_combined"]  for r in overall]

bins_all = sorted({r["bin_id"] for r in all_rows if r["bin_id"] >= 0})
bin_colors = {0: "steelblue", 1: "darkorange", 2: "mediumseagreen",
              3: "orchid",    4: "tomato"}

fig, axes = plt.subplots(3, 1, figsize=(8, 11), sharex=True)

ax = axes[0]
ax.plot(lams_o, ameds, "o-",  color="steelblue", label="median Δψ")
ax.plot(lams_o, amns,  "s--", color="steelblue", alpha=0.6, label="mean Δψ")
ax.axhline(ameds[0], color="steelblue", lw=0.8, ls=":", alpha=0.5,
           label=f"SplineMPE baseline = {ameds[0]:.2f}°")
ax.set_ylabel("Angular error [deg]")
ax.legend(fontsize=8); ax.grid(True, alpha=0.3)
ax.set_title("SplineMPE + λ·NaI (Gulliver)"
             + (" [true time]" if args.true_time else ""))

ax = axes[1]
for b in bins_all:
    brows = [r for r in all_rows if r["bin_id"] == b]
    bls   = [r["lam"]        for r in brows]
    bmed  = [r["ang_median"] for r in brows]
    ene   = brows[0]["med_energy_GeV"] / 1e3 if brows else 0
    ax.plot(bls, bmed, "o-", color=bin_colors.get(b, "gray"),
            label=f"bin {b}  {ene:.2f} TeV")
ax.set_ylabel("Median Δψ per bin [deg]")
ax.legend(fontsize=7, ncol=2); ax.grid(True, alpha=0.3)

ax = axes[2]
ax.plot(lams_o, l_comb,  "o-",  color="tomato",       label="loss: mean(Δψ)+mean(d⊥)/100")
ax.plot(lams_o, l_huber, "s--", color="mediumpurple",  label=f"loss: Huber(δ={args.huber_delta}°)")
if l_comb:
    best_comb  = lams_o[int(np.argmin(l_comb))]
    best_huber = lams_o[int(np.argmin(l_huber))]
    ax.axvline(best_comb,  color="tomato",      ls="--", lw=1.2,
               label=f"opt λ (combined) = {best_comb:.2f}")
    ax.axvline(best_huber, color="mediumpurple", ls="--", lw=1.2,
               label=f"opt λ (Huber) = {best_huber:.2f}")
ax.set_ylabel("Loss"); ax.set_xlabel("λ (NaI weight)")
ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

if len(lams_o) > 2:
    pos_l = [l for l in lams_o if l > 0]
    if pos_l and max(pos_l) / min(pos_l) > 10:
        for a in axes:
            a.set_xscale("symlog", linthresh=0.5)

plt.tight_layout()
out_png = args.out.replace(".csv", ".png")
fig.savefig(out_png, dpi=150, bbox_inches="tight")
plt.close()
print(f"Plot: {out_png}")
if l_comb:
    print(f"\nOptimal λ (combined loss) = {best_comb:.3f}")
    print(f"Optimal λ (Huber loss)    = {best_huber:.3f}")
