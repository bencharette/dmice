"""
steamshovel_artists.py — DMice LineFit Visualization Artists

Two Steamshovel PyArtist classes for real-time visual comparison of:
  - ICLineFitArtist:    IC-only LineFit (analytic, charge-weighted centroid pivot)
  - PivotLineFitArtist: DM-Ice Pivot LineFit (DM-Ice detector as fixed reference)

Written in pure Python stdlib (no numpy) for compatibility with steamshovel's
embedded Python 3.9 interpreter.

Loaded automatically via ~/.steamshovel/startup.py.
To add to the scene from the Python console:
    window.gl.scenario.add('ICLineFitArtist')
    window.gl.scenario.add('PivotLineFitArtist')
"""

import math

from icecube.shovelart import (
    PyArtist, PyQColor, PyQFont, RangeSetting,
    vec3d, ConstantVec3d,
)
from icecube import dataclasses, icetray
from icecube.icetray import I3Units

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

C_M_NS = 0.2998  # speed of light in m/ns

DMICE_POS = {
    "det1": ( 31.25,   -72.93,  -511.05),
    "det2": (-334.80, -424.50,  -511.26),
}

PULSE_KEYS = [
    "SRTInIcePulses", "SplitInIcePulses", "InIcePulses",
    "TWCMuonPulseSeriesReco", "OfflinePulses",
    "OnlineL2_CleanedMuonPulses", "SplitUncleanedInIcePulses",
    "UncleanedInIcePulses",
]

# ---------------------------------------------------------------------------
# Pure-Python vector helpers
# ---------------------------------------------------------------------------

def _dot(a, b):
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]

def _sub(a, b):
    return (a[0]-b[0], a[1]-b[1], a[2]-b[2])

def _add(a, b):
    return (a[0]+b[0], a[1]+b[1], a[2]+b[2])

def _scale(s, a):
    return (s*a[0], s*a[1], s*a[2])

def _norm(a):
    return math.sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2])

def _unit(a):
    n = _norm(a)
    if n == 0:
        return (0.0, 0.0, 0.0)
    return (a[0]/n, a[1]/n, a[2]/n)

# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

def _get_dom_positions(frame):
    """Build omkey -> (x,y,z) dict from I3Geometry in frame."""
    if "I3Geometry" not in frame:
        return {}
    geo = frame["I3Geometry"]
    dom_pos = {}
    for omkey, omgeo in geo.omgeo.items():
        p = omgeo.position
        dom_pos[omkey] = (p.x, p.y, p.z)
    return dom_pos


def _extract_hits(frame, dom_pos):
    """
    Extract charge-weighted pulses from the first available pulse key.
    Returns (xs, ys, zs, ts, ws) lists, or None if no pulses found.
    """
    pulse_map = None
    for key in PULSE_KEYS:
        if key in frame:
            raw = frame[key]
            pulse_map = raw.apply(frame) if hasattr(raw, "apply") else raw
            break
    if pulse_map is None:
        return None

    xs, ys, zs, ts, ws = [], [], [], [], []
    for omkey, pulses in pulse_map:
        if omkey not in dom_pos:
            continue
        px, py, pz = dom_pos[omkey]
        for pulse in pulses:
            charge = pulse.charge if pulse.charge > 0 else 1.0
            xs.append(px); ys.append(py); zs.append(pz)
            ts.append(pulse.time); ws.append(charge)

    if not xs:
        return None
    return xs, ys, zs, ts, ws


def _weighted_mean(vals, ws):
    W = sum(ws)
    if W == 0:
        return 0.0
    return sum(v*w for v, w in zip(vals, ws)) / W


def _run_ic_linefit(xs, ys, zs, ts, ws):
    """
    Standard analytic IC-only LineFit (centre-of-gravity pivot).
    Returns dict with dx,dy,dz, speed_m_ns, cx,cy,cz, or None.
    """
    W = sum(ws)
    if W == 0:
        return None
    cx = _weighted_mean(xs, ws)
    cy = _weighted_mean(ys, ws)
    cz = _weighted_mean(zs, ws)
    t_bar = _weighted_mean(ts, ws)

    dts = [t - t_bar for t in ts]
    denom = sum(w * dt*dt for w, dt in zip(ws, dts))
    if denom == 0:
        return None

    vx = sum(w*dt*(x-cx) for w,dt,x in zip(ws,dts,xs)) / denom
    vy = sum(w*dt*(y-cy) for w,dt,y in zip(ws,dts,ys)) / denom
    vz = sum(w*dt*(z-cz) for w,dt,z in zip(ws,dts,zs)) / denom

    speed = math.sqrt(vx*vx + vy*vy + vz*vz)
    if speed == 0:
        return None
    return dict(dx=vx/speed, dy=vy/speed, dz=vz/speed,
                speed_m_ns=speed, cx=cx, cy=cy, cz=cz)


def _compute_dmice_hit_time(xs, ys, zs, ts, ws, dm_pos, mc_dir):
    """
    Expected DM-Ice transit time from the MC-truth direction.
    Projects (dm_pos - charge_centroid) onto mc_dir to get travel offset.
    """
    cx = _weighted_mean(xs, ws)
    cy = _weighted_mean(ys, ws)
    cz = _weighted_mean(zs, ws)
    t_bar = _weighted_mean(ts, ws)
    d_proj = _dot(_sub(dm_pos, (cx, cy, cz)), mc_dir)
    return t_bar + d_proj / C_M_NS


def _run_pivot_linefit(xs, ys, zs, ts, ws, dm_pos, t_dm_ns):
    """
    DM-Ice Pivot LineFit: DM-Ice detector is the fixed reference in space+time.
      dt_i = t_i - t_dm,   dr_i = r_i - r_dm
      v = sum(w_i dt_i dr_i) / sum(w_i dt_i^2)
    """
    dts  = [t - t_dm_ns for t in ts]
    drxs = [x - dm_pos[0] for x in xs]
    drys = [y - dm_pos[1] for y in ys]
    drzs = [z - dm_pos[2] for z in zs]

    denom = sum(w*dt*dt for w,dt in zip(ws,dts))
    if denom == 0:
        return None

    vx = sum(w*dt*dr for w,dt,dr in zip(ws,dts,drxs)) / denom
    vy = sum(w*dt*dr for w,dt,dr in zip(ws,dts,drys)) / denom
    vz = sum(w*dt*dr for w,dt,dr in zip(ws,dts,drzs)) / denom

    speed = math.sqrt(vx*vx + vy*vy + vz*vz)
    if speed == 0:
        return None
    return dict(dx=vx/speed, dy=vy/speed, dz=vz/speed, speed_m_ns=speed)


def _get_mc_dir(frame):
    """
    Get primary muon travel direction from I3MCTree.
    BLO convention: momentum direction stored directly in primary.dir (no pi-flip).
    """
    for key in ("I3MCTree", "I3MCTree_preMuonProp"):
        if key in frame:
            tree = frame[key]
            primaries = tree.primaries
            if primaries:
                p = primaries[0]
                return (p.dir.x, p.dir.y, p.dir.z)
    return None


def _dir_to_zen_azi(dx, dy, dz):
    """
    Travel direction -> IceCube source zenith/azimuth (degrees).
    zen = arccos(-dz)  (source direction, anti-momentum convention).
    """
    zen = math.degrees(math.acos(max(-1.0, min(1.0, -dz))))
    azi = math.degrees(math.atan2(dy, dx)) % 360.0
    return zen, azi


def _add_track_arrow(output, centroid, direction, half_len, color, line_width):
    """Draw a directional arrow centred at centroid along direction."""
    start = _sub(centroid, _scale(half_len, direction))
    end   = _add(centroid, _scale(half_len, direction))
    a = output.addArrow(
        ConstantVec3d(vec3d(start[0], start[1], start[2])),
        ConstantVec3d(vec3d(end[0],   end[1],   end[2])),
        20.0 * I3Units.degree,
        60.0 * I3Units.m,
    )
    a.setColor(color)
    a.setLineWidth(line_width)
    return a


# ---------------------------------------------------------------------------
# Artist 1 — IC-only LineFit
# ---------------------------------------------------------------------------

class ICLineFitArtist(PyArtist):
    """
    IC-only LineFit track reconstruction visualized in Steamshovel.

    Reads InIcePulses, computes the analytic LineFit, and draws the
    reconstructed track direction as a blue arrow through the
    charge-weighted centroid of the hit DOMs.
    """

    numRequiredKeys = 0

    def __init__(self):
        PyArtist.__init__(self)
        self.defineSettings({
            "color":         PyQColor(70, 130, 180, 230),
            "track_length":  RangeSetting(200.0, 5000.0, 480, 2000.0),
            "line_width":    RangeSetting(0.5, 8.0, 150, 2.0),
            "show_centroid": True,
        })

    def description(self):
        return "IC-only LineFit"

    def create(self, frame, output):
        dom_pos = _get_dom_positions(frame)
        if not dom_pos:
            return

        hits = _extract_hits(frame, dom_pos)
        if hits is None:
            return
        xs, ys, zs, ts, ws = hits

        if len(xs) < 4:
            return

        fit = _run_ic_linefit(xs, ys, zs, ts, ws)
        if fit is None:
            return

        color    = self.setting("color")
        half_len = self.setting("track_length") / 2.0
        lw       = self.setting("line_width")

        centroid  = (fit["cx"], fit["cy"], fit["cz"])
        direction = (fit["dx"], fit["dy"], fit["dz"])

        _add_track_arrow(output, centroid, direction, half_len, color, lw)

        if self.setting("show_centroid"):
            s = output.addSphere(14.0, vec3d(centroid[0], centroid[1], centroid[2]))
            s.setColor(color)

        zen, azi = _dir_to_zen_azi(fit["dx"], fit["dy"], fit["dz"])
        font = PyQFont()
        font.pointSize = 10
        output.addTextOverlay(
            "IC LineFit   zen={:.1f}deg  azi={:.1f}deg  v={:.3f}c".format(
                zen, azi, fit["speed_m_ns"] / C_M_NS),
            font,
        )


# ---------------------------------------------------------------------------
# Artist 2 — DM-Ice Pivot LineFit
# ---------------------------------------------------------------------------

class PivotLineFitArtist(PyArtist):
    """
    DM-Ice Pivot LineFit track reconstruction visualized in Steamshovel.

    Uses the DM-Ice detector as a fixed space-time pivot point. The
    expected DM-Ice transit time is computed from MC-truth direction
    (I3MCTree, BLO convention: momentum stored directly, no pi-flip).

    Both arrows share the IC charge-weighted centroid as anchor, so the
    angular gap between them directly shows the improvement from DM-Ice.
    """

    numRequiredKeys = 0

    def __init__(self):
        PyArtist.__init__(self)
        self.defineSettings({
            "color":        PyQColor(210, 60, 60, 230),
            "track_length": RangeSetting(200.0, 5000.0, 480, 2000.0),
            "line_width":   RangeSetting(0.5, 8.0, 150, 2.0),
            "show_pivot":   True,
            "pivot_size":   RangeSetting(5.0, 100.0, 190, 30.0),
            "detector":     "auto",
        })

    def description(self):
        return "DM-Ice Pivot LineFit"

    def create(self, frame, output):
        dom_pos = _get_dom_positions(frame)
        if not dom_pos:
            return

        hits = _extract_hits(frame, dom_pos)
        if hits is None:
            return
        xs, ys, zs, ts, ws = hits

        if len(xs) < 4:
            return

        mc_dir = _get_mc_dir(frame)
        if mc_dir is None:
            return

        # Choose DM-Ice detector
        det_setting = self.setting("detector")
        if det_setting in ("det1", "det2"):
            dm_name = det_setting
        elif "BLO_DetId" in frame:
            det_tag = str(frame["BLO_DetId"].value)
            dm_name = det_tag if det_tag in DMICE_POS else "det1"
        else:
            def ca(dpos):
                cx = _weighted_mean(xs, ws)
                cy = _weighted_mean(ys, ws)
                cz = _weighted_mean(zs, ws)
                dp = _sub(dpos, (cx, cy, cz))
                proj = _dot(dp, mc_dir)
                perp = _sub(dp, _scale(proj, mc_dir))
                return _norm(perp)
            dm_name = "det1" if ca(DMICE_POS["det1"]) <= ca(DMICE_POS["det2"]) else "det2"

        dm_pos = DMICE_POS[dm_name]

        t_dm_ns = _compute_dmice_hit_time(xs, ys, zs, ts, ws, dm_pos, mc_dir)
        fit = _run_pivot_linefit(xs, ys, zs, ts, ws, dm_pos, t_dm_ns)
        if fit is None:
            return

        color    = self.setting("color")
        half_len = self.setting("track_length") / 2.0
        lw       = self.setting("line_width")

        cx = _weighted_mean(xs, ws)
        cy = _weighted_mean(ys, ws)
        cz = _weighted_mean(zs, ws)
        centroid  = (cx, cy, cz)
        direction = (fit["dx"], fit["dy"], fit["dz"])

        _add_track_arrow(output, centroid, direction, half_len, color, lw)

        if self.setting("show_pivot"):
            s = output.addSphere(
                self.setting("pivot_size"),
                vec3d(dm_pos[0], dm_pos[1], dm_pos[2])
            )
            s.setColor(color)

        zen, azi = _dir_to_zen_azi(fit["dx"], fit["dy"], fit["dz"])
        font = PyQFont()
        font.pointSize = 10
        output.addTextOverlay(
            "Pivot LineFit  [{}]   zen={:.1f}deg  azi={:.1f}deg  v={:.3f}c".format(
                dm_name, zen, azi, fit["speed_m_ns"] / C_M_NS),
            font,
        )


# ---------------------------------------------------------------------------
# Auto-register when exec'd from startup.py
# ---------------------------------------------------------------------------
try:
    window.gl.scenario.registerArtist(ICLineFitArtist)
    window.gl.scenario.registerArtist(PivotLineFitArtist)
    print('[DMice] ICLineFitArtist + PivotLineFitArtist registered.')
except NameError:
    pass
