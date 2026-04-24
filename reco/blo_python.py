"""
blo_python.py — Python port of BlueLightOrchestra.jl

Replicates the three main BLO functions:
    propagate(particle, dist_km)  -> list of energy-loss ParticleStates
    run_ppc(particle, losses)     -> list of PhotonHit namedtuples
    process_hits(hits)            -> numpy arrays of per-DOM hit data

Dependencies:
    proposal  (pip install proposal, already on cobalt)
    numpy
    PPC CPU binary  (~/dmice_work/ppc_cpu/ppc)

Usage:
    from blo_python import ParticleState, propagate, run_ppc, process_hits, load_detector, smt8_trigger

    p = ParticleState(energy_GeV=2000, pos_m=[0,0,-1300], dir=[0,0,-1], pid=13, time_ns=0)
    losses = propagate(p, dist_km=1.5)
    hits   = run_ppc(p, losses)
    doms   = process_hits(hits)
    print(doms['x'], doms['z'], doms['t'], doms['nhits'])
    triggered, t_trig = smt8_trigger(doms)
    print("SMT8 triggered:", triggered, "at t =", t_trig, "ns")
"""

import os
import subprocess
import tempfile
import math
import random
import string
import numpy as np
from dataclasses import dataclass, field
from typing import List

# ── Paths ─────────────────────────────────────────────────────────────────────

_THIS_DIR    = os.path.dirname(os.path.abspath(__file__))
_BLO_DIR     = os.path.join(_THIS_DIR, "BlueLightOrchestra.jl")
_RESOURCE_DIR = os.path.join(_BLO_DIR, "resources")
_GEO_FILE    = os.environ.get(
    "BLO_GEO_FILE",
    os.path.join(_RESOURCE_DIR, "geofiles", "icecube_with_dmice.geo")
)
_PPC_TABLES  = os.environ.get(
    "BLO_PPC_TABLES",
    os.path.join(_RESOURCE_DIR, "PPC_tables", "south_pole")
)

# PPC binary: override via BLO_PPC_EXE env var, else default cobalt path
PPC_EXE = os.environ.get(
    "BLO_PPC_EXE",
    os.path.expanduser("~/dmice_work/ppc_cpu/ppc")
)

# PPC applies a fixed z offset to convert IceCube coords to its internal frame
PPC_MAGIC_Z = 1948.07  # metres

# ── PROPOSAL particle map ─────────────────────────────────────────────────────

_PDG_TO_PROPOSAL = {
    13:   "MuMinus",
    -13:  "MuPlus",
    11:   "EMinus",
    -11:  "EPlus",
    22:   "Gamma",
    15:   "TauMinus",
    -15:  "TauPlus",
}

_PROPOSAL_TO_F2K = {
    "BremsstrahlungLoss":        "brems",
    "IonizationLoss":            "delta",
    "EpairProductionLoss":       "epair",
    "PhotoNuclearLoss":          "hadr",
    "MupairProductionLoss":      "mupair",
    "WeakInteractionLoss":       "hadr",
    "ContinuousLoss":            "delta",
}

# ── Data types ────────────────────────────────────────────────────────────────

@dataclass
class ParticleState:
    energy_GeV: float
    pos_m:      List[float]   # [x, y, z] in metres
    dir:        List[float]   # unit vector [dx, dy, dz]
    pid:        int           # PDG ID (13 = muon-)
    time_ns:    float = 0.0

@dataclass
class PhotonHit:
    string_id:  int
    sensor_id:  int
    time_ns:    float
    source:     str
    wavelength_nm: float = 0.0
    om_zen:     float = 0.0
    om_az:      float = 0.0
    photon_zen: float = 0.0
    photon_az:  float = 0.0

# ── Geometry ──────────────────────────────────────────────────────────────────

def load_detector(geo_file=_GEO_FILE):
    """
    Parse a BLO .geo file and return arrays of DOM positions and IDs.
    Returns dict with keys: x, y, z (m), string_id, sensor_id.
    """
    rows = []
    with open(geo_file) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 5:
                continue
            try:
                rows.append([float(p) for p in parts])
            except ValueError:
                continue  # skip header lines
    arr = np.array(rows)
    return {
        "x":         arr[:, 0],
        "y":         arr[:, 1],
        "z":         arr[:, 2],
        "string_id": arr[:, 3].astype(int),
        "sensor_id": arr[:, 4].astype(int),
    }

# ── PROPOSAL propagation ──────────────────────────────────────────────────────

def _build_propagator(pid, E_cut_MeV=500.0, v_cut=0.1, cont_rand=True):
    import proposal as pp

    pdef_cls = getattr(pp.particle, _PDG_TO_PROPOSAL[pid] + "Def")
    pdef     = pdef_cls()
    medium   = pp.medium.Ice()
    cuts     = pp.EnergyCutSettings(E_cut_MeV, v_cut, cont_rand)

    xsec = pp.crosssection.make_std_crosssection(pdef, medium, cuts, True)

    coll = pp.PropagationUtilityCollection()
    coll.displacement = pp.make_displacement(xsec, True)
    coll.interaction  = pp.make_interaction(xsec, True)
    coll.time         = pp.make_time(xsec, pdef, True)
    utility  = pp.PropagationUtility(collection=coll)

    geometry = pp.geometry.Sphere(pp.Cartesian3D(0, 0, 0), 1e20)
    density  = pp.density_distribution.density_homogeneous(medium.mass_density)

    return pp.Propagator(pdef, [(geometry, utility, density)])


def propagate(particle: ParticleState, dist_km: float) -> List[ParticleState]:
    """
    Propagate particle through ice using PROPOSAL.
    Returns list of energy-loss ParticleStates (stochastic + discretised continuous).
    """
    import proposal as pp

    prop = _build_propagator(particle.pid)

    init = pp.particle.ParticleState()
    init.energy    = particle.energy_GeV * 1e3          # MeV
    init.position  = pp.Cartesian3D(*(np.array(particle.pos_m) * 100))  # cm
    init.direction = pp.Cartesian3D(*particle.dir)
    init.time      = particle.time_ns * 1e-9             # s

    result = prop.propagate(init, dist_km * 1e5)        # cm

    losses = []

    # stochastic losses
    for l in result.stochastic_losses():
        type_name = type(l).__name__
        losses.append(ParticleState(
            energy_GeV = l.energy * 1e-3,
            pos_m      = [l.position.x * 0.01, l.position.y * 0.01, l.position.z * 0.01],
            dir        = [l.direction.x, l.direction.y, l.direction.z],
            pid        = -1,   # interaction type stored in _f2k_type below
            time_ns    = l.time * 1e9,
        ))
        losses[-1]._f2k_type  = _PROPOSAL_TO_F2K.get(type_name, "hadr")

    # continuous losses: discretise at 1 m steps
    dx_m = 1.0
    for l in result.continuous_losses():
        x0 = np.array([l.start_position.x, l.start_position.y, l.start_position.z]) * 0.01
        x1 = np.array([l.end_position.x,   l.end_position.y,   l.end_position.z])   * 0.01
        d  = np.array([l.direction_initial.x, l.direction_initial.y, l.direction_initial.z])
        dist_m = np.linalg.norm(x1 - x0)
        n_pts  = max(1, int(math.ceil(dist_m / dx_m)))
        dE_GeV = (l.energy * 1e-3) / n_pts
        dt_ns  = (dist_m / n_pts) / 0.2998  # ns per metre in ice
        for i in range(1, n_pts + 1):
            pos = x0 + i * (dist_m / (n_pts + 1)) * d
            p = ParticleState(
                energy_GeV = dE_GeV,
                pos_m      = pos.tolist(),
                dir        = d.tolist(),
                pid        = -1,
                time_ns    = l.time_initial * 1e9 + i * dt_ns,
            )
            p._f2k_type = "delta"
            losses.append(p)

    return losses

# ── PPC interface ─────────────────────────────────────────────────────────────

def _random_serial():
    return "0x" + "".join(random.choices("0123456789abcdef", k=12))

def _random_mac():
    return "".join(random.choices(string.ascii_uppercase + string.digits, k=8))

def _dir_to_zen_azi(d):
    dx, dy, dz = d
    zen = math.atan2(math.hypot(dx, dy), dz)
    azi = math.atan2(dy, dx)
    return zen, azi


def _write_geo_f2k(tmp_dir, det):
    lines = []
    for x, y, z, sid, oid in zip(det["x"], det["y"], det["z"], det["string_id"], det["sensor_id"]):
        mac    = _random_mac()
        serial = _random_serial()
        lines.append(f"{mac}\t{serial}\t{x}\t{y}\t{z}\t{sid}\t{oid}")
    with open(os.path.join(tmp_dir, "geo-f2k"), "w") as f:
        f.write("\n".join(lines) + "\n")


def _write_event_f2k(fname, particle: ParticleState, losses: List[ParticleState]):
    zen, azi = _dir_to_zen_azi(particle.dir)
    x, y, z  = particle.pos_m

    lines = [f"EM 0 1 0 0 0 0 "]
    lines.append(
        f"MC E {particle.energy_GeV} "
        f"x {x} y {y} z {z + PPC_MAGIC_Z} "
        f"theta {zen} phi {azi}"
    )
    for loss in losses:
        lzen, lazi = _dir_to_zen_azi(loss.dir)
        lx, ly, lz = loss.pos_m
        f2k_type = getattr(loss, "_f2k_type", "delta")
        dt = loss.time_ns - particle.time_ns
        lines.append(
            f"TR 0 0 {f2k_type} "
            f"{lx} {ly} {lz + PPC_MAGIC_Z} "
            f"{lzen} {lazi} 0 "
            f"{loss.energy_GeV} {dt}"
        )
    lines.append("EE")
    with open(fname, "w") as f:
        f.write("\n".join(lines) + "\n")


def _read_hits_f2k(fname) -> List[PhotonHit]:
    hits = []
    src  = ""
    with open(fname) as f:
        for line in f:
            line = line.strip()
            if line.startswith("TR"):
                parts = line.split()
                src = parts[3] if len(parts) > 3 else ""
            elif line.startswith("HIT"):
                parts = line.split()
                vals  = [float(v) for v in parts[1:]]
                hit = PhotonHit(
                    string_id     = int(vals[0]),
                    sensor_id     = int(vals[1]),
                    time_ns       = vals[2],
                    source        = src,
                    wavelength_nm = vals[3] if len(vals) > 3 else 0.0,
                    om_zen        = vals[4] if len(vals) > 4 else 0.0,
                    om_az         = vals[5] if len(vals) > 5 else 0.0,
                    photon_zen    = vals[6] if len(vals) > 6 else 0.0,
                    photon_az     = vals[7] if len(vals) > 7 else 0.0,
                )
                hits.append(hit)
    return hits


def run_ppc(particle: ParticleState, losses: List[ParticleState],
            suppress_error=True, geo_file=_GEO_FILE) -> List[PhotonHit]:
    """
    Run the PPC CPU binary to simulate photon hits.
    Returns list of PhotonHit objects.
    """
    det = load_detector(geo_file)

    with tempfile.TemporaryDirectory(prefix="blo_ppc_") as tmp_dir:
        # copy ice tables
        import shutil
        for fname in os.listdir(_PPC_TABLES):
            shutil.copy(os.path.join(_PPC_TABLES, fname), os.path.join(tmp_dir, fname))

        _write_geo_f2k(tmp_dir, det)

        f2k_in  = os.path.join(tmp_dir, "event_in.f2k")
        f2k_out = os.path.join(tmp_dir, "event_out.f2k")
        _write_event_f2k(f2k_in, particle, losses)

        env = os.environ.copy()
        env["PPCTABLESDIR"] = tmp_dir

        with open(f2k_in)  as fin, \
             open(f2k_out, "w") as fout:
            stderr_dest = subprocess.DEVNULL if suppress_error else None
            subprocess.run(
                [PPC_EXE, "0"],
                stdin=fin, stdout=fout, stderr=stderr_dest,
                env=env, check=True
            )

        return _read_hits_f2k(f2k_out)


# ── SMT8 trigger ─────────────────────────────────────────────────────────────

# Parameters from IceCube instrumentation paper (arXiv:1612.05093), Table 8
SMT8_N_HLC      = 8        # required DOM multiplicity
SMT8_WINDOW_NS  = 5000.0   # 5 µs sliding window
SMT8_IC_STRINGS = (1, 86)  # IceCube in-ice string ID range (inclusive)


def _hlc_mask(doms, lc_window_ns=1000.0, lc_neighbor_range=2):
    """
    Return a boolean mask of DOMs that satisfy the Hard Local Coincidence (HLC)
    condition: a DOM counts as HLC if at least one string-neighbor within
    ±lc_neighbor_range positions on the same string also has a first hit within
    ±lc_window_ns of its own first hit.

    Parameters
    ----------
    doms              : dict from process_hits()
    lc_window_ns      : coincidence window (default ±1000 ns, standard IceCube LC)
    lc_neighbor_range : how many adjacent DOMs up/down the string to check (default 2,
                        i.e. nearest and second-nearest neighbours)
    """
    string_ids = doms["string_id"]
    sensor_ids = doms["sensor_id"]
    t          = doms["t"]
    n          = len(t)
    hlc        = np.zeros(n, dtype=bool)

    # Build a lookup: (string_id, sensor_id) -> first-hit time
    hit_map = {(string_ids[i], sensor_ids[i]): t[i] for i in range(n)}

    for i in range(n):
        sid = string_ids[i]
        oid = sensor_ids[i]
        t_i = t[i]
        for d in range(1, lc_neighbor_range + 1):
            for neighbor_oid in (oid - d, oid + d):
                t_nb = hit_map.get((sid, neighbor_oid))
                if t_nb is not None and abs(t_nb - t_i) <= lc_window_ns:
                    hlc[i] = True
                    break
            if hlc[i]:
                break

    return hlc


def smt8_trigger(doms, n_hits=SMT8_N_HLC, window_ns=SMT8_WINDOW_NS,
                 ic_only=True, lc_window_ns=1000.0, lc_neighbor_range=2):
    """
    Apply the IceCube Simple Multiplicity Trigger (SMT8) to simulated DOM hits.

    Requires n_hits HLC (Hard Local Coincidence) hits within a window_ns sliding
    window.  A DOM is HLC if a string-neighbor (within ±lc_neighbor_range positions
    on the same string) also has a hit within ±lc_window_ns.  This matches the real
    IceCube trigger which counts only neighbor-coincident hits toward the multiplicity.

    Parameters
    ----------
    doms              : dict returned by process_hits()
    n_hits            : multiplicity threshold (default 8)
    window_ns         : sliding trigger window in ns (default 5000 = 5 µs)
    ic_only           : if True, restrict to IceCube in-ice strings (IDs 1–86)
    lc_window_ns      : LC coincidence window in ns (default ±1000 ns)
    lc_neighbor_range : neighbor range on string for LC check (default ±2,
                        nearest and second-nearest neighbours)

    Returns
    -------
    triggered  : bool
    t_trigger  : float or None — start time (ns) of the first window that fires
    """
    string_ids = doms["string_id"]
    t          = doms["t"]

    if len(t) == 0:
        return False, None

    if ic_only:
        lo, hi = SMT8_IC_STRINGS
        ic_mask = (string_ids >= lo) & (string_ids <= hi)
        # Restrict doms dict for the LC check so neighbours are only IC DOMs
        doms_ic = {k: v[ic_mask] for k, v in doms.items()}
    else:
        doms_ic = doms

    hlc = _hlc_mask(doms_ic, lc_window_ns=lc_window_ns,
                    lc_neighbor_range=lc_neighbor_range)
    t_hlc = doms_ic["t"][hlc]

    if len(t_hlc) < n_hits:
        return False, None

    t_sorted = np.sort(t_hlc)

    # Sliding window: fires when n_hits consecutive (by time) HLC hits span ≤ window_ns
    for i in range(len(t_sorted) - n_hits + 1):
        if t_sorted[i + n_hits - 1] - t_sorted[i] <= window_ns:
            return True, float(t_sorted[i])

    return False, None


# ── Hit processing ────────────────────────────────────────────────────────────

def process_hits(hits: List[PhotonHit], geo_file=_GEO_FILE) -> dict:
    """
    Aggregate per-photon hits into per-DOM summary.
    Returns dict: x, y, z (m), t (ns, first hit), nhits, string_id, sensor_id.
    """
    det = load_detector(geo_file)

    # group by (string_id, sensor_id)
    dom_data = {}
    for h in hits:
        key = (h.string_id, h.sensor_id)
        if key not in dom_data:
            dom_data[key] = []
        dom_data[key].append(h.time_ns)

    if not dom_data:
        empty = np.array([])
        return {k: empty for k in ("x", "y", "z", "t", "nhits", "string_id", "sensor_id")}

    keys      = list(dom_data.keys())
    str_ids   = np.array([k[0] for k in keys])
    sen_ids   = np.array([k[1] for k in keys])
    t_first   = np.array([min(dom_data[k]) for k in keys])
    nhits     = np.array([len(dom_data[k]) for k in keys])

    # look up positions from geometry
    x_out = np.zeros(len(keys))
    y_out = np.zeros(len(keys))
    z_out = np.zeros(len(keys))
    for i, (sid, oid) in enumerate(keys):
        mask = (det["string_id"] == sid) & (det["sensor_id"] == oid)
        idx  = np.where(mask)[0]
        if len(idx):
            x_out[i] = det["x"][idx[0]]
            y_out[i] = det["y"][idx[0]]
            z_out[i] = det["z"][idx[0]]

    return {
        "x":         x_out,
        "y":         y_out,
        "z":         z_out,
        "t":         t_first,
        "nhits":     nhits,
        "string_id": str_ids,
        "sensor_id": sen_ids,
    }
