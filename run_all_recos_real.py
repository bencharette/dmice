#!/usr/bin/env python3
"""
run_all_recos_real.py

Run all reconstruction methods on real DM-Ice coincidence data (2012-2019).

Methods:
  LineFit (std)        — already in frame, re-extracted
  Pivot LineFit        — DM-Ice-anchored, using calibrated μ correction
  MPEFit (std seed)    — Pandel MPE from LineFit seed
  MPEFit (pivot seed)  — Pandel MPE from Pivot LF seed
  SPEFit (std seed)    — Pandel SPE1st from LineFit seed
  SPEFit (pivot seed)  — Pandel SPE1st from Pivot LF seed

Calibrated timing model: μ=+280 ns, σ=81 ns (from real 2012 coincidences).

Output CSV: ~/dmice_work/output/real_all_recos.csv

Run on NPX via condor (see run_all_recos_real.sub) or on cobalt:
  /cvmfs/.../env-shell.sh python3 ~/dmice/run_all_recos_real.py
"""

import os, math
import numpy as np

# ── Constants ─────────────────────────────────────────────────────────────────
C_M_NS  = 0.2998
N_ICE   = 1.3195
THETA_C = math.acos(1.0 / N_ICE)

# DM-Ice positions [m] in IceCube coordinates
DMICE_POS = {
    "det1": np.array([ 31.25,  -72.93, -511.05]),
    "det2": np.array([-334.80, -424.50, -511.26]),
}

# Calibrated timing model (fit to real 2012 coincidence data)
MU_NS    =  280.0   # ns  — NaI scintillation + cable delays
SIGMA_NS =   81.0   # ns
D_MAX    =   15.0   # m   — max d_perp to use DM-Ice constraint

# Strings 1-86 are IceCube proper; 87-88 are DM-Ice
IC_STRINGS = set(range(1, 87))

# ── Paths ─────────────────────────────────────────────────────────────────────
SPICEMIE_DIR    = "/cvmfs/icecube.opensciencegrid.org/data/photon-tables/SPICEMie"
SPICEMIE_DRIVER = os.path.join(SPICEMIE_DIR, "driverfiles")
GCD_FILE  = "/cvmfs/icecube.opensciencegrid.org/data/GCD/GeoCalibDetectorStatus_2013.56429_V1.i3.gz"
IN_FILE   = "/data/user/bcharett/dmice_coincidences_2011_2022/all_dmice_coincidences_2011_2022_fixed.i3"
OUT_DIR   = os.path.expanduser("~/dmice_work/output")
OUT_CSV   = os.path.join(OUT_DIR, "real_all_recos.csv")
os.makedirs(OUT_DIR, exist_ok=True)

# ── Pulse key priority ────────────────────────────────────────────────────────
PULSE_PRIORITY = ["OfflinePulses", "SRTInIcePulses", "ReextractedInIcePulses", "InIcePulses"]
IC_PULSES      = "RealICPulses"
PIV_LF_KEY     = "RealPivotLF"       # pivot seeded from LineFit direction
PIV_MPE_KEY    = "RealPivotMPE_LF"  # pivot seeded from MPEFit(std) direction
MPE_STD_KEY    = "RealMPE_Std"
MPE_PIV_KEY    = "RealMPE_Piv"      # MPEFit seeded from LineFit-pivot
MPE_PIV2_KEY   = "RealMPE_Piv2"     # MPEFit seeded from MPEFit-pivot
SPE_STD_KEY    = "RealSPE_Std"
SPE_PIV_KEY    = "RealSPE_Piv"
TRUNC_E_KEY    = "RealTruncatedEnergy"
PHOT_SERVICE   = "PhotonicsServiceMu"

# ── Geometry helpers ──────────────────────────────────────────────────────────
def t_geometric(track_pos, track_dir, t0, dm_pos):
    r      = np.asarray(dm_pos) - np.asarray(track_pos)
    d_hat  = np.asarray(track_dir) / np.linalg.norm(track_dir)
    s      = float(np.dot(r, d_hat))
    d_perp = math.sqrt(max(0.0, float(np.dot(r, r)) - s**2))
    t_pca  = t0 + s / C_M_NS
    t_geo  = t_pca + (d_perp / (C_M_NS * math.sin(THETA_C)) if d_perp > 0.01 else 0.0)
    return t_geo, d_perp


def pivot_linefit_dm(xs, ys, zs, ts, ws, dm_pos, dm_t_corrected, seed_dir):
    """LineFit anchored to μ-corrected DM-Ice hit time."""
    if len(xs) < 4:
        return None, None
    W  = sum(ws)
    cx = sum(x*w for x,w in zip(xs,ws)) / W
    cy = sum(y*w for y,w in zip(ys,ws)) / W
    cz = sum(z*w for z,w in zip(zs,ws)) / W
    tb = sum(t*w for t,w in zip(ts,ws)) / W

    # Project DM-Ice pos onto seed direction to get transit time anchor
    d  = np.asarray(seed_dir, dtype=float)
    d  = d / np.linalg.norm(d)
    r  = np.array([dm_pos[0]-cx, dm_pos[1]-cy, dm_pos[2]-cz])
    t_dm = tb + float(np.dot(r, d)) / C_M_NS

    # Shift time origin to DM-Ice anchor
    dts  = [t - dm_t_corrected for t in ts]
    drxs = [x - dm_pos[0] for x in xs]
    drys = [y - dm_pos[1] for y in ys]
    drzs = [z - dm_pos[2] for z in zs]

    den  = sum(w*dt*dt for w,dt in zip(ws,dts))
    if not den:
        return None, None
    vx = sum(w*dt*dr for w,dt,dr in zip(ws,dts,drxs)) / den
    vy = sum(w*dt*dr for w,dt,dr in zip(ws,dts,drys)) / den
    vz = sum(w*dt*dr for w,dt,dr in zip(ws,dts,drzs)) / den
    spd = math.sqrt(vx**2 + vy**2 + vz**2)
    if not spd:
        return None, None

    # Disambiguate direction: ensure consistent with seed
    if vx*seed_dir[0] + vy*seed_dir[1] + vz*seed_dir[2] < 0:
        vx, vy, vz = -vx, -vy, -vz
    direction = (vx/spd, vy/spd, vz/spd)
    # t0 at DM-Ice: back-project along new direction to event reference
    d_new = np.array(direction)
    s_new = float(np.dot(np.array(dm_pos) - np.array([cx,cy,cz]), d_new))
    t0_pivot = dm_t_corrected - s_new / C_M_NS
    return direction, t0_pivot


# ── IceTray modules ───────────────────────────────────────────────────────────
from icecube import icetray, dataio, dataclasses, recclasses
import icecube.lilliput.segments
from icecube.icetray import I3Tray

rows     = []
seen     = set()   # deduplication: (run_id, event_id)
geo_omgeo = {}     # filled at G-frame

class SetupModule(icetray.I3Module):
    """Split IC-only pulses; extract DM-Ice hit time; compute Pivot LF."""
    def __init__(self, ctx):
        super().__init__(ctx)

    def Configure(self): pass

    def Geometry(self, frame):
        global geo_omgeo
        geo = frame["I3Geometry"]
        geo_omgeo = {k: v for k, v in geo.omgeo.items()}
        self.PushFrame(frame)

    def Physics(self, frame):
        hdr = frame["I3EventHeader"]
        # Include sub_event_stream so NullSplit and InIceSplit are not
        # collapsed together (2016+ data has multiple P-frames per DAQ event)
        uid = (hdr.run_id, hdr.event_id,
               getattr(hdr, 'sub_event_stream', ''))

        # ── Skip non-muon sub-event streams ──────────────────────────────
        stream = getattr(hdr, 'sub_event_stream', '')
        if stream not in ('', 'in_ice', 'InIceSplit'):
            return   # NullSplit, IceTopSplit, SLOPSplit etc. — no muon reco

        # ── Deduplicate ───────────────────────────────────────────────────
        if uid in seen:
            return   # drop duplicate
        seen.add(uid)

        # ── Find best pulse key ────────────────────────────────────────────
        pulse_key = None
        for k in PULSE_PRIORITY:
            if k in frame:
                pulse_key = k
                break
        if pulse_key is None:
            self.PushFrame(frame)
            return

        # ── Split IC-only pulses (strings 1-86) ───────────────────────────
        try:
            all_pulses = dataclasses.I3RecoPulseSeriesMap.from_frame(frame, pulse_key)
        except Exception:
            self.PushFrame(frame)
            return

        ic_map = dataclasses.I3RecoPulseSeriesMap()
        for omk, plist in all_pulses:
            if omk.string in IC_STRINGS:
                ic_map[omk] = plist

        if len(ic_map) < 4:
            self.PushFrame(frame)
            return

        frame[IC_PULSES] = ic_map

        # ── DM-Ice hit time (event-local ns) ──────────────────────────────
        if "DMIce_detection_time" not in frame or "DMIce_detector" not in frame:
            self.PushFrame(frame)
            return

        event_daq = hdr.start_time.utc_daq_time
        dm_t_ns   = (frame["DMIce_detection_time"].value - event_daq) * 0.1
        dm_t_corrected = dm_t_ns - MU_NS

        det_raw = str(frame["DMIce_detector"])
        det_key = "det1" if "det1" in det_raw else "det2"
        dm_pos  = DMICE_POS[det_key]

        # ── Get LineFit for seed (fall back to PoleMuonLinefit for 2017+) ──
        lf_key = None
        for k in ("LineFit", "PoleMuonLinefit"):
            if k in frame:
                lf_key = k
                break
        if lf_key is None:
            self.PushFrame(frame)
            return
        lf = frame[lf_key]
        seed_dir = (lf.dir.x, lf.dir.y, lf.dir.z)

        # ── Check d_perp — only anchor pivot if within D_MAX ──────────────
        lf_pos = np.array([lf.pos.x, lf.pos.y, lf.pos.z])
        _, d_perp = t_geometric(lf_pos, seed_dir, lf.time, dm_pos)
        frame["RealDPerp"] = dataclasses.I3Double(d_perp)
        frame["RealDetKey"] = dataclasses.I3String(det_key)
        frame["RealDMtNs"]  = dataclasses.I3Double(dm_t_ns)

        if d_perp > D_MAX:
            self.PushFrame(frame)
            return

        # ── Extract DOM positions for pivot ───────────────────────────────
        xs, ys, zs, ts, ws = [], [], [], [], []
        for omk, plist in ic_map:
            if omk not in geo_omgeo:
                continue
            pos = geo_omgeo[omk].position
            for p in plist:
                xs.append(pos.x); ys.append(pos.y); zs.append(pos.z)
                ts.append(p.time); ws.append(p.charge)

        if len(xs) < 4:
            self.PushFrame(frame)
            return

        direction, t0_pivot = pivot_linefit_dm(
            xs, ys, zs, ts, ws, dm_pos, dm_t_corrected, seed_dir)
        if direction is None:
            self.PushFrame(frame)
            return

        pp = dataclasses.I3Particle()
        pp.dir    = dataclasses.I3Direction(*direction)
        pp.pos    = dataclasses.I3Position(*dm_pos)
        pp.time   = float(t0_pivot)
        pp.fit_status = dataclasses.I3Particle.FitStatus.OK
        frame[PIV_LF_KEY] = pp

        # ── MPEFit-seeded pivot (computed after MPEFit std runs) ───────────
        # Stored as a flag; actual computation done in a second module after
        # MPEFit(std) has run. Store DM-Ice anchor info for reuse.
        frame["RealDMpos"]   = dataclasses.I3Position(*dm_pos)
        frame["RealDMtCorr"] = dataclasses.I3Double(dm_t_corrected)

        self.PushFrame(frame)


class MPEPivotModule(icetray.I3Module):
    """Compute pivot anchored from MPEFit(std) direction instead of LineFit."""
    def __init__(self, ctx):
        super().__init__(ctx)

    def Configure(self): pass

    def Physics(self, frame):
        # Only run if we have MPEFit(std) and the DM-Ice anchor info
        if (MPE_STD_KEY not in frame or "RealDMtCorr" not in frame
                or "RealDMpos" not in frame or IC_PULSES not in frame):
            self.PushFrame(frame)
            return

        mpe = frame[MPE_STD_KEY]
        if mpe.fit_status != dataclasses.I3Particle.FitStatus.OK:
            self.PushFrame(frame)
            return

        seed_dir = (mpe.dir.x, mpe.dir.y, mpe.dir.z)
        dm_pos   = np.array([frame["RealDMpos"].x,
                             frame["RealDMpos"].y,
                             frame["RealDMpos"].z])
        dm_t_corrected = frame["RealDMtCorr"].value

        # Extract IC pulse positions
        xs, ys, zs, ts, ws = [], [], [], [], []
        try:
            ic_map = dataclasses.I3RecoPulseSeriesMap.from_frame(frame, IC_PULSES)
        except Exception:
            self.PushFrame(frame)
            return
        for omk, plist in ic_map:
            if omk not in geo_omgeo:
                continue
            pos = geo_omgeo[omk].position
            for p in plist:
                xs.append(pos.x); ys.append(pos.y); zs.append(pos.z)
                ts.append(p.time); ws.append(p.charge)

        if len(xs) < 4:
            self.PushFrame(frame)
            return

        direction, t0_pivot = pivot_linefit_dm(
            xs, ys, zs, ts, ws, dm_pos, dm_t_corrected, seed_dir)
        if direction is None:
            self.PushFrame(frame)
            return

        pp = dataclasses.I3Particle()
        pp.dir    = dataclasses.I3Direction(*direction)
        pp.pos    = dataclasses.I3Position(*dm_pos)
        pp.time   = float(t0_pivot)
        pp.fit_status = dataclasses.I3Particle.FitStatus.OK
        frame[PIV_MPE_KEY] = pp
        self.PushFrame(frame)


class ScorerModule(icetray.I3Module):
    """Extract results from all reco keys into rows list."""
    def __init__(self, ctx):
        super().__init__(ctx)

    def Configure(self): pass

    def Physics(self, frame):
        if IC_PULSES not in frame:
            self.PushFrame(frame)
            return

        hdr     = frame["I3EventHeader"]
        year    = hdr.start_time.utc_year
        det_key = frame["RealDetKey"].value if "RealDetKey" in frame else "unknown"
        d_perp  = frame["RealDPerp"].value  if "RealDPerp"  in frame else float("nan")
        dm_t_ns = frame["RealDMtNs"].value  if "RealDMtNs"  in frame else float("nan")

        try:
            pulse_map = dataclasses.I3RecoPulseSeriesMap.from_frame(frame, IC_PULSES)
            n_doms = len(pulse_map)
            n_hits = sum(len(v) for v in pulse_map.values())
        except Exception:
            n_doms = 0
            n_hits = 0

        # Energy — try TruncatedEnergy first, fall back to existing keys
        energy_GeV = float("nan")
        for ekey in [TRUNC_E_KEY + "_ORIG_Muon",
                     "MPEFitTruncatedEnergy_SPICEMie_ORIG_Muon",
                     "TruncatedEnergy_SPICEMie_ORIG_Muon"]:
            try:
                e = frame[ekey].energy
                if e > 0:
                    energy_GeV = e
                    break
            except Exception:
                pass

        def get_zen_azi(key):
            if key not in frame:
                return float("nan"), float("nan")
            p = frame[key]
            if p.fit_status != dataclasses.I3Particle.FitStatus.OK:
                return float("nan"), float("nan")
            return math.degrees(p.dir.zenith), math.degrees(p.dir.azimuth)

        lf_zen,       lf_azi       = get_zen_azi("LineFit")
        piv_lf_zen,   piv_lf_azi   = get_zen_azi(PIV_LF_KEY)
        piv_mpe_zen,  piv_mpe_azi  = get_zen_azi(PIV_MPE_KEY)
        mpe_std_zen,  mpe_std_azi  = get_zen_azi(MPE_STD_KEY)
        mpe_piv_zen,  mpe_piv_azi  = get_zen_azi(MPE_PIV_KEY)
        mpe_piv2_zen, mpe_piv2_azi = get_zen_azi(MPE_PIV2_KEY)
        spe_std_zen,  spe_std_azi  = get_zen_azi(SPE_STD_KEY)
        spe_piv_zen,  spe_piv_azi  = get_zen_azi(SPE_PIV_KEY)

        rows.append(dict(
            year=year,
            run_id=hdr.run_id,
            event_id=hdr.event_id,
            detector=det_key,
            d_perp_m=round(d_perp, 2),
            dm_t_ns=round(dm_t_ns, 1),
            n_doms_ic=n_doms,
            n_hits_ic=n_hits,
            energy_GeV=round(energy_GeV, 2),
            lf_zen=round(lf_zen, 4),            lf_azi=round(lf_azi, 4),
            piv_lf_zen=round(piv_lf_zen, 4),    piv_lf_azi=round(piv_lf_azi, 4),
            piv_mpe_zen=round(piv_mpe_zen, 4),  piv_mpe_azi=round(piv_mpe_azi, 4),
            mpe_std_zen=round(mpe_std_zen, 4),  mpe_std_azi=round(mpe_std_azi, 4),
            mpe_piv_zen=round(mpe_piv_zen, 4),  mpe_piv_azi=round(mpe_piv_azi, 4),
            mpe_piv2_zen=round(mpe_piv2_zen,4), mpe_piv2_azi=round(mpe_piv2_azi,4),
            spe_std_zen=round(spe_std_zen, 4),  spe_std_azi=round(spe_std_azi, 4),
            spe_piv_zen=round(spe_piv_zen, 4),  spe_piv_azi=round(spe_piv_azi, 4),
        ))
        self.PushFrame(frame)


# ── Build tray ────────────────────────────────────────────────────────────────
from icecube import truncated_energy

tray = I3Tray()

tray.AddService("I3PhotonicsServiceFactory", PHOT_SERVICE,
    PhotonicsTopLevelDirectory = SPICEMIE_DIR,
    DriverFileDirectory        = SPICEMIE_DRIVER,
    PhotonicsLevel2DriverFile  = "mu_photorec.list",
    PhotonicsTableSelection    = 2,
    ServiceName                = PHOT_SERVICE,
)

tray.AddModule("I3Reader", FilenameList=[GCD_FILE, IN_FILE])

tray.AddModule(SetupModule, "Setup")

# MPEFit — std seed (LineFit)  [must run before MPEPivotModule]
tray.Add(icecube.lilliput.segments.I3SinglePandelFitter,
    fitname = MPE_STD_KEY,
    domllh  = "MPE",
    pulses  = IC_PULSES,
    seeds   = ["LineFit"],
    If      = lambda f: IC_PULSES in f and "LineFit" in f,
)

# MPEFit — pivot seed
tray.Add(icecube.lilliput.segments.I3SinglePandelFitter,
    fitname = MPE_PIV_KEY,
    domllh  = "MPE",
    pulses  = IC_PULSES,
    seeds   = [PIV_LF_KEY],
    If      = lambda f: IC_PULSES in f and PIV_LF_KEY in f,
)

# Pivot anchored from MPEFit(std) direction — runs after MPEFit std
tray.AddModule(MPEPivotModule, "MPEPivot")

# MPEFit — MPE-pivot seed (second-stage pivot)
tray.Add(icecube.lilliput.segments.I3SinglePandelFitter,
    fitname = MPE_PIV2_KEY,
    domllh  = "MPE",
    pulses  = IC_PULSES,
    seeds   = [PIV_MPE_KEY],
    If      = lambda f: IC_PULSES in f and PIV_MPE_KEY in f,
)

# SPEFit — std seed (LineFit)
tray.Add(icecube.lilliput.segments.I3SinglePandelFitter,
    fitname = SPE_STD_KEY,
    domllh  = "SPE1st",
    pulses  = IC_PULSES,
    seeds   = ["LineFit"],
    If      = lambda f: IC_PULSES in f and "LineFit" in f,
)

# SPEFit — pivot seed
tray.Add(icecube.lilliput.segments.I3SinglePandelFitter,
    fitname = SPE_PIV_KEY,
    domllh  = "SPE1st",
    pulses  = IC_PULSES,
    seeds   = [PIV_LF_KEY],
    If      = lambda f: IC_PULSES in f and PIV_LF_KEY in f,
)

# TruncatedEnergy — seeded from MPEFit pivot (best reco we have)
tray.AddModule("I3TruncatedEnergy",
    RecoPulsesName         = IC_PULSES,
    RecoParticleName       = MPE_PIV_KEY,
    ResultParticleName     = TRUNC_E_KEY,
    I3PhotonicsServiceName = PHOT_SERVICE,
    UseRDE                 = True,
    If = lambda f: IC_PULSES in f and MPE_PIV_KEY in f,
)

tray.AddModule(ScorerModule, "Scorer")

tray.Execute()
tray.Finish()

# ── Write CSV ─────────────────────────────────────────────────────────────────
import pandas as pd
df = pd.DataFrame(rows)
df.to_csv(OUT_CSV, index=False)

print(f"\nDone: {len(df)} events (deduplicated)")
print(f"  With DM-Ice pivot (d⊥ < {D_MAX}m): {df.piv_lf_zen.notna().sum()}")
print(f"  Years: {sorted(df.year.unique())}")
print(f"  Per year: {df.groupby('year').size().to_dict()}")
print(f"CSV: {OUT_CSV}")
