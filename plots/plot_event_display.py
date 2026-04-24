#!/usr/bin/env python3
"""
plot_event_display.py — Event display plots for BLO-simulated muon events.

Reproduces the BLO README-style two-panel figure:
  Left:  x–z scatter of all DOMs (gray) + hit DOMs coloured by first-hit time,
         sized by nhits, with muon track (magenta dashed) and injection point (X).
  Right: Ice absorption profile vs depth from the SPICE icemodel.

Usage:
    python plot_event_display.py [--npz PATH] [--geo PATH] [--icemodel PATH] [--out DIR]

Defaults:
    npz      ~/dmice/output/blo_muons_200hits_rerun.npz
    geo      ~/dmice/BLO/icecube_with_dmice.geo
    icemodel ~/dmice/BlueLightOrchestra.jl/resources/PPC_tables/south_pole/icemodel.dat (if present)
             else ~/dmice/BLO/resources/PPC_tables/south_pole/icemodel.dat
"""

import os
import sys
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

# ── Paths ─────────────────────────────────────────────────────────────────────

def _find(candidates):
    for p in candidates:
        ep = os.path.expanduser(p)
        if os.path.exists(ep):
            return ep
    return None

DEFAULT_NPZ = _find([
    "~/dmice/output/blo_muons_200hits_rerun.npz",
    "~/dmice_work/output/blo_muons_200hits_rerun.npz",
])
DEFAULT_GEO = _find([
    "~/dmice/BLO/icecube_with_dmice.geo",
    "~/dmice/handoff/icecube_with_dmice.geo",
    "~/dmice/BlueLightOrchestra.jl/resources/geofiles/icecube_with_dmice.geo",
])
DEFAULT_ICE = _find([
    "~/dmice/BlueLightOrchestra.jl/resources/PPC_tables/south_pole/icemodel.dat",
    "~/dmice/BLO/resources/PPC_tables/south_pole/icemodel.dat",
])

# ── DM-Ice detector positions in BLO coords [m] ──────────────────────────────

DMICE_POS = {
    "det1": np.array([ 31.25,  -72.93, -2459.12]),
    "det2": np.array([-334.80, -424.50, -2459.33]),
}

C_M_NS = 0.2998   # speed of light m/ns

# ── Geometry loader ───────────────────────────────────────────────────────────

def load_geo(geo_file):
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
                rows.append([float(p) for p in parts[:5]])
            except ValueError:
                continue
    arr = np.array(rows)
    return {
        "x":         arr[:, 0],   # metres
        "y":         arr[:, 1],
        "z":         arr[:, 2],   # metres, negative = deeper
        "string_id": arr[:, 3].astype(int),
        "sensor_id": arr[:, 4].astype(int),
    }

# ── Ice model loader ──────────────────────────────────────────────────────────

def load_icemodel(icemodel_file):
    """
    Returns z_km (BLO convention, negative) and be400 (eff. scattering coeff m^-1).
    Columns: depth_from_surface_m, be400, adust400, delta_tau
    BLO z = -depth_from_surface_m / 1000
    """
    data = []
    with open(icemodel_file) as f:
        for line in f:
            parts = line.split()
            if len(parts) >= 2:
                try:
                    depth_m = float(parts[0])
                    be400   = float(parts[1])
                    adust   = float(parts[2]) if len(parts) > 2 else 0.0
                    if adust < 10.0:   # skip 999 sentinel values
                        data.append((-depth_m / 1000.0, be400))
                except ValueError:
                    continue
    data = np.array(data)
    return data[:, 0], data[:, 1]   # z_km, be400

# ── LineFit reconstructions (ported from run_linefit.py) ─────────────────────

def _wm(vals, ws):
    W = sum(ws); return sum(v*w for v,w in zip(vals,ws))/W if W else 0.0

def _dot(a, b): return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]
def _sub(a, b): return (a[0]-b[0], a[1]-b[1], a[2]-b[2])
def _scale(s,a): return (s*a[0], s*a[1], s*a[2])

def ic_linefit(xs, ys, zs, ts, ws):
    if len(xs) < 4: return None
    cx,cy,cz = _wm(xs,ws),_wm(ys,ws),_wm(zs,ws)
    tb = _wm(ts,ws)
    dts = [t-tb for t in ts]
    den = sum(w*dt*dt for w,dt in zip(ws,dts))
    if not den: return None
    vx = sum(w*dt*(x-cx) for w,dt,x in zip(ws,dts,xs))/den
    vy = sum(w*dt*(y-cy) for w,dt,y in zip(ws,dts,ys))/den
    vz = sum(w*dt*(z-cz) for w,dt,z in zip(ws,dts,zs))/den
    spd = math.sqrt(vx*vx+vy*vy+vz*vz)
    if not spd: return None
    return dict(dx=vx/spd, dy=vy/spd, dz=vz/spd, cx=cx, cy=cy, cz=cz)

def pivot_linefit(xs, ys, zs, ts, ws, dm_pos_m, mc_dir):
    cx,cy,cz = _wm(xs,ws),_wm(ys,ws),_wm(zs,ws)
    tb = _wm(ts,ws)
    d_proj = _dot(_sub(dm_pos_m,(cx,cy,cz)), mc_dir)
    t_dm   = tb + d_proj/C_M_NS
    dts  = [t-t_dm for t in ts]
    drxs = [x-dm_pos_m[0] for x in xs]
    drys = [y-dm_pos_m[1] for y in ys]
    drzs = [z-dm_pos_m[2] for z in zs]
    den  = sum(w*dt*dt for w,dt in zip(ws,dts))
    if not den: return None
    vx = sum(w*dt*dr for w,dt,dr in zip(ws,dts,drxs))/den
    vy = sum(w*dt*dr for w,dt,dr in zip(ws,dts,drys))/den
    vz = sum(w*dt*dr for w,dt,dr in zip(ws,dts,drzs))/den
    spd = math.sqrt(vx*vx+vy*vy+vz*vz)
    if not spd: return None
    return dict(dx=vx/spd, dy=vy/spd, dz=vz/spd, cx=cx, cy=cy, cz=cz)

import math

# ── Track projection helper ───────────────────────────────────────────────────

def track_xz(pos_m, dir_xyz, z_range_km):
    """
    Return x [km], z [km] arrays for the muon track projected onto x-z plane,
    clipped to z_range_km = (z_min, z_max).
    """
    x0, y0, z0 = np.array(pos_m) / 1000.0   # km
    dx, dy, dz = dir_xyz

    if abs(dz) < 1e-9:
        return np.array([x0, x0]), np.array(z_range_km)

    t_min = (z_range_km[0] - z0) / dz
    t_max = (z_range_km[1] - z0) / dz
    if t_min > t_max:
        t_min, t_max = t_max, t_min

    t_vals = np.linspace(t_min, t_max, 200)
    x_vals = x0 + dx * t_vals
    z_vals = z0 + dz * t_vals
    return x_vals, z_vals

# ── Reconstruct injection position from stored direction ──────────────────────

# Downgoing sim geometry: inject from above, aimed through DM-Ice detector
INJECT_Z_KM = -1.3   # km — top of IC86

def reconstruct_injection(zen_blo_rad, azi_rad, target_det_id=0):
    """Reconstruct injection position [m] aimed through the target DM-Ice detector."""
    dz = np.cos(zen_blo_rad)
    sin_zen = np.sin(zen_blo_rad)
    dx = sin_zen * np.cos(azi_rad)
    dy = sin_zen * np.sin(azi_rad)

    # target DM-Ice position in km
    dm_km = np.array(list(DMICE_POS.values())[target_det_id]) / 1000.0

    if abs(dz) < 1e-9:
        return [0.0, 0.0, INJECT_Z_KM * 1e3], [dx, dy, dz]

    # back-project from DM-Ice position to injection height
    t_km = (dm_km[2] - INJECT_Z_KM) / dz
    x0_m = (dm_km[0] - dx * t_km) * 1e3
    y0_m = (dm_km[1] - dy * t_km) * 1e3
    z0_m = INJECT_Z_KM * 1e3
    return [x0_m, y0_m, z0_m], [dx, dy, dz]

# ── Plot one event ────────────────────────────────────────────────────────────

def plot_event(ev_idx, data, det, ice_z_km, ice_be, out_dir):
    ene_GeV = float(data["energy_GeV"][ev_idx])
    zen_rad = float(data["zenith_rad"][ev_idx])
    azi_rad = float(data["azimuth_rad"][ev_idx])

    dom_x   = np.array(data["dom_x"][ev_idx])     # m
    dom_y   = np.array(data["dom_y"][ev_idx])
    dom_z   = np.array(data["dom_z"][ev_idx])
    dom_t   = np.array(data["dom_t"][ev_idx])      # ns
    nhits   = np.array(data["dom_nhits"][ev_idx])

    tgt_id_early = int(data["target_det"][ev_idx]) if "target_det" in data else 0
    pos_m, dir_xyz = reconstruct_injection(zen_rad, azi_rad, target_det_id=tgt_id_early)

    # ── z range for plot ──────────────────────────────────────────────────────
    z_all_km = det["z"] / 1000.0
    z_min_km = min(z_all_km.min(), -2.55)
    z_max_km = max(z_all_km.max() + 0.05, -1.15)

    # ── Time colormap ─────────────────────────────────────────────────────────
    if len(dom_t) > 0:
        t_us    = dom_t / 1e3                      # ns → μs
        t_min   = np.percentile(t_us, 2)
        t_max   = np.percentile(t_us, 98)
        norm    = mcolors.Normalize(vmin=t_min, vmax=t_max)
        cmap    = plt.cm.jet
        colors  = cmap(norm(t_us))
    else:
        t_us = np.array([])
        norm = mcolors.Normalize(vmin=0, vmax=1)
        cmap = plt.cm.jet
        colors = np.array([]).reshape(0, 4)

    # ── Marker sizes proportional to nhits (sqrt scaling, capped to avoid
    #    DM-Ice strings dominating the display with millions of PPC photons)
    SIZE_CAP = 50.0
    if len(nhits) > 0:
        sizes = (4.0 * np.minimum(nhits, SIZE_CAP) ** (1.0 / 3.0)) ** 2
    else:
        sizes = np.array([])

    # ── Figure layout ─────────────────────────────────────────────────────────
    fig = plt.figure(figsize=(6, 7))
    gs  = fig.add_gridspec(2, 1, height_ratios=[1, 20], hspace=0.05)

    ax_title = fig.add_subplot(gs[0])
    ax_left  = fig.add_subplot(gs[1])
    ax_right = None  # removed

    # ── Title row ─────────────────────────────────────────────────────────────
    ene_label  = f"E = {ene_GeV/1e3:.1f} TeV" if ene_GeV >= 1e3 else f"E = {ene_GeV:.0f} GeV"
    zen_std    = 180.0 - np.degrees(zen_rad)
    bin_label  = f"  bin {int(data['bin_id'][ev_idx])}" if "bin_id" in data else ""
    tgt_id     = int(data["target_det"][ev_idx]) if "target_det" in data else -1
    tgt_label  = f"  → det{'12'[tgt_id]}" if tgt_id >= 0 else ""
    ax_title.set_axis_off()
    ax_title.text(0.5, 0.5,
                  rf"$\times$  Injected: $\mu^-$ with {ene_label}  "
                  rf"zen$_{{std}}$={zen_std:.1f}°{bin_label}{tgt_label}",
                  ha="center", va="center", fontsize=11,
                  color="magenta", transform=ax_title.transAxes)

    # ── Left panel: x–z event display ────────────────────────────────────────
    # All DOMs (gray background)
    ax_left.scatter(det["x"] / 1000.0, det["z"] / 1000.0,
                    s=2, c="lightgray", zorder=1, linewidths=0)

    # Hit DOMs
    if len(dom_x) > 0:
        sc = ax_left.scatter(dom_x / 1000.0, dom_z / 1000.0,
                             s=sizes, c=colors, zorder=3,
                             edgecolors="none", alpha=0.85)

    # Muon track projection (x–z plane)
    trk_x, trk_z = track_xz(pos_m, dir_xyz, (z_min_km, z_max_km))
    ax_left.plot(trk_x, trk_z, "m--", lw=1.2, zorder=2, alpha=0.8)

    # Injection point
    inj_x_km = pos_m[0] / 1000.0
    inj_z_km = pos_m[2] / 1000.0
    ax_left.scatter([inj_x_km], [inj_z_km], marker="x", s=120,
                    c="magenta", linewidths=2.0, zorder=4)

    # DM-Ice detector markers
    for dname, dpos in DMICE_POS.items():
        dx_km = dpos[0] / 1000.0
        dz_km = dpos[2] / 1000.0
        color = "deepskyblue" if dname == "det1" else "orange"
        ax_left.scatter([dx_km], [dz_km], marker="*", s=180,
                        c=color, zorder=5, edgecolors="k", linewidths=0.5,
                        label=dname)

    # ── LineFit reconstructions ───────────────────────────────────────────────
    if len(dom_x) >= 4:
        xs = dom_x.tolist(); ys = dom_y.tolist(); zs = dom_z.tolist()
        ts = dom_t.tolist(); ws = nhits.tolist()

        # IC LineFit
        ic_fit = ic_linefit(xs, ys, zs, ts, ws)
        if ic_fit:
            trk_x, trk_z = track_xz(
                [ic_fit["cx"], ic_fit["cy"], ic_fit["cz"]],
                [ic_fit["dx"], ic_fit["dy"], ic_fit["dz"]],
                (z_min_km, z_max_km)
            )
            ax_left.plot(trk_x, trk_z, "b-", lw=1.5, zorder=6,
                         alpha=0.8, label="IC LineFit")

        # Pivot LineFit (needs MC truth direction and target DM-Ice position)
        dx_mc = math.sin(zen_rad) * math.cos(azi_rad)
        dy_mc = math.sin(zen_rad) * math.sin(azi_rad)
        dz_mc = math.cos(zen_rad)
        mc_dir = (dx_mc, dy_mc, dz_mc)
        dm_key = "det1" if tgt_id != 1 else "det2"
        piv_fit = pivot_linefit(xs, ys, zs, ts, ws, DMICE_POS[dm_key].tolist(), mc_dir)
        if piv_fit:
            trk_x, trk_z = track_xz(
                [piv_fit["cx"], piv_fit["cy"], piv_fit["cz"]],
                [piv_fit["dx"], piv_fit["dy"], piv_fit["dz"]],
                (z_min_km, z_max_km)
            )
            ax_left.plot(trk_x, trk_z, "r-", lw=1.5, zorder=6,
                         alpha=0.8, label="Pivot LineFit")

    ax_left.legend(loc="upper right", fontsize=7, framealpha=0.7)
    ax_left.set_xlabel("x [km]")
    ax_left.set_ylabel("z [km]")
    ax_left.set_xlim(-0.6, 0.6)
    ax_left.set_ylim(z_min_km, z_max_km)

    # Right panel removed (ax_right is None)

    # ── Colorbar ──────────────────────────────────────────────────────────────
    if len(t_us) > 0:
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar = fig.colorbar(sm, ax=ax_left,
                            orientation="horizontal", fraction=0.05, pad=0.08)
        cbar.set_label(r"time [$\mu$s]")

    plt.suptitle("", y=0.98)

    out_path = os.path.join(out_dir, f"event_display_ev{ev_idx:03d}.png")
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out_path}")
    return out_path

# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--npz",      default=DEFAULT_NPZ)
    parser.add_argument("--geo",      default=DEFAULT_GEO)
    parser.add_argument("--icemodel", default=DEFAULT_ICE)
    parser.add_argument("--out",      default=os.path.expanduser("~/dmice/output"))
    parser.add_argument("--per-bin",  type=int, default=None,
                        help="Plot top N events by nhits per bin_id (requires bin_id in npz)")
    args = parser.parse_args()

    if args.npz is None or not os.path.exists(args.npz):
        sys.exit(f"ERROR: npz not found: {args.npz}")
    if args.geo is None or not os.path.exists(args.geo):
        sys.exit(f"ERROR: geo file not found: {args.geo}")

    os.makedirs(args.out, exist_ok=True)

    print(f"Loading events: {args.npz}")
    data = np.load(args.npz, allow_pickle=True)
    n_events = len(data["energy_GeV"])
    print(f"  {n_events} events")

    print(f"Loading geometry: {args.geo}")
    det = load_geo(args.geo)
    print(f"  {len(det['x'])} DOMs")

    ice_z_km = ice_be = None
    if args.icemodel and os.path.exists(args.icemodel):
        print(f"Loading ice model: {args.icemodel}")
        ice_z_km, ice_be = load_icemodel(args.icemodel)

    # Select which event indices to plot
    if args.per_bin is not None and "bin_id" in data:
        bin_ids  = data["bin_id"]
        n_hits   = data["n_hits"]
        indices  = []
        for b in np.unique(bin_ids):
            mask     = np.where(bin_ids == b)[0]
            top_n    = mask[np.argsort(n_hits[mask])[::-1][:args.per_bin]]
            indices.extend(sorted(top_n))
        print(f"  Plotting top {args.per_bin}/bin → {len(indices)} events total")
    else:
        indices = list(range(n_events))

    for i in indices:
        ene = float(data["energy_GeV"][i])
        print(f"\n[Event {i}]  E={ene/1e3:.3f} TeV  "
              f"nhits={int(data['n_hits'][i])}  ndoms={int(data['n_doms'][i])}")
        plot_event(i, data, det, ice_z_km, ice_be, args.out)

    print("\nDone.")

if __name__ == "__main__":
    main()
