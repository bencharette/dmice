#!/usr/bin/env python3
"""
plot_event_display_itermpe.py

Event display showing MC truth track and IterMPE reconstruction only.
Reads hit data from a BLO NPZ and reco directions from the IterMPE CSV.

Usage:
    python3 plot_event_display_itermpe.py [--npz PATH] [--recos CSV] [--out DIR]
"""

import os, sys, argparse, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

parser = argparse.ArgumentParser()
parser.add_argument("--npz", default=os.path.expanduser(
    "~/dmice_work/output/muons_binned_5bins_1000pbin_repacked.npz"))
parser.add_argument("--recos", default=os.path.expanduser(
    "~/dmice_work/output/itermpe_events.csv"))
parser.add_argument("--geo", default=os.path.expanduser(
    "~/dmice/BLO/icecube_with_dmice.geo"))
parser.add_argument("--out", default=os.path.expanduser(
    "~/dmice_work/output/event_displays_itermpe"))
args = parser.parse_args()

import csv as _csv

recos = {}
with open(args.recos) as _f:
    for _row in _csv.DictReader(_f):
        _idx = int(_row["ev_idx"])
        recos[_idx] = {k: float(v) for k, v in _row.items() if k != "ev_idx"}

class _RowProxy(dict):
    def __getattr__(self, name):
        return self[name]

recos = {k: _RowProxy(v) for k, v in recos.items()}

DMICE_POS = {
    "det1": np.array([ 31.25,  -72.93, -2459.12]),
    "det2": np.array([-334.80, -424.50, -2459.33]),
}
C_M_NS    = 0.2998
INJECT_Z_KM = -1.3
SIZE_CAP  = 50.0

# ── Geometry ──────────────────────────────────────────────────────────────────

def load_geo(path):
    rows = []
    with open(path) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 5:
                continue
            try:
                rows.append([float(p) for p in parts[:5]])
            except ValueError:
                continue
    arr = np.array(rows)
    return {"x": arr[:,0], "y": arr[:,1], "z": arr[:,2],
            "string_id": arr[:,3].astype(int), "sensor_id": arr[:,4].astype(int)}

det = load_geo(args.geo)
data = np.load(args.npz, allow_pickle=True)

def load_ragged(key, i):
    if f"{key}_flat" in data:
        flat    = data[f"{key}_flat"]
        offsets = data[f"{key}_offsets"]
        return flat[offsets[i]:offsets[i+1]]
    return data[key][i]

def track_xz(pos_m, dir_xyz, z_range_km):
    x0, y0, z0 = np.array(pos_m) / 1000.0
    dx, dy, dz = dir_xyz
    if abs(dz) < 1e-9:
        return np.array([x0, x0]), np.array(z_range_km)
    t_min = (z_range_km[0] - z0) / dz
    t_max = (z_range_km[1] - z0) / dz
    if t_min > t_max:
        t_min, t_max = t_max, t_min
    t_vals = np.linspace(t_min, t_max, 200)
    return x0 + dx * t_vals, z0 + dz * t_vals

os.makedirs(args.out, exist_ok=True)

for ev_idx, row in recos.items():
    ene = float(data["energy_GeV"][ev_idx])
    zen = float(data["zenith_rad"][ev_idx])
    azi = float(data["azimuth_rad"][ev_idx])
    tgt_id = int(data["target_det"][ev_idx]) if "target_det" in data else 0
    bin_id = int(data["bin_id"][ev_idx]) if "bin_id" in data else -1

    dom_x  = np.array(load_ragged("dom_x",     ev_idx))
    dom_y  = np.array(load_ragged("dom_y",     ev_idx))
    dom_z  = np.array(load_ragged("dom_z",     ev_idx))
    dom_t  = np.array(load_ragged("dom_t",     ev_idx))
    nhits  = np.array(load_ragged("dom_nhits", ev_idx))

    # MC truth direction and injection point
    mc_dx, mc_dy, mc_dz = float(row.mc_dx), float(row.mc_dy), float(row.mc_dz)
    dm_km = np.array(list(DMICE_POS.values())[tgt_id]) / 1000.0
    dz_mc = mc_dz
    if abs(dz_mc) > 1e-9:
        t_km = (dm_km[2] - INJECT_Z_KM) / dz_mc
        inj_pos_m = np.array([
            (dm_km[0] - mc_dx * t_km) * 1e3,
            (dm_km[1] - mc_dy * t_km) * 1e3,
            INJECT_Z_KM * 1e3,
        ])
    else:
        inj_pos_m = np.array([0., 0., INJECT_Z_KM * 1e3])

    # IterMPE direction
    im_dx, im_dy, im_dz = float(row.itermpe_dx), float(row.itermpe_dy), float(row.itermpe_dz)
    err = float(row.ang_err_deg)

    z_all_km = det["z"] / 1000.0
    z_min_km = min(z_all_km.min(), -2.55)
    z_max_km = max(z_all_km.max() + 0.05, -1.15)

    # Time colormap
    if len(dom_t) > 0:
        t_us  = dom_t / 1e3
        norm  = mcolors.Normalize(vmin=np.percentile(t_us, 2),
                                   vmax=np.percentile(t_us, 98))
        cmap  = plt.cm.jet
        colors = cmap(norm(t_us))
    else:
        t_us = np.array([])
        norm = mcolors.Normalize(vmin=0, vmax=1)
        cmap = plt.cm.jet
        colors = np.array([]).reshape(0, 4)

    sizes = (4.0 * np.minimum(nhits, SIZE_CAP) ** (1.0/3.0)) ** 2 if len(nhits) else np.array([])

    fig = plt.figure(figsize=(6, 7))
    gs  = fig.add_gridspec(2, 1, height_ratios=[1, 20], hspace=0.05)
    ax_title = fig.add_subplot(gs[0])
    ax       = fig.add_subplot(gs[1])

    ene_label = f"E = {ene/1e3:.1f} TeV" if ene >= 1e3 else f"E = {ene:.0f} GeV"
    zen_std   = 180.0 - np.degrees(zen)
    det_label = f"det{'12'[tgt_id]}"
    ax_title.set_axis_off()
    ax_title.text(0.5, 0.5,
        rf"$\times$  {ene_label}  zen$_{{std}}$={zen_std:.1f}°  "
        rf"bin {bin_id} → {det_label}  |  IterMPE err={err:.2f}°",
        ha="center", va="center", fontsize=10,
        color="magenta", transform=ax_title.transAxes)

    # All DOMs (gray)
    ax.scatter(det["x"]/1000., det["z"]/1000., s=2, c="lightgray", zorder=1, linewidths=0)

    # Hit DOMs
    if len(dom_x):
        ax.scatter(dom_x/1000., dom_z/1000., s=sizes, c=colors,
                   zorder=3, edgecolors="none", alpha=0.85)

    # MC truth track (magenta dashed)
    tx, tz = track_xz(inj_pos_m, (mc_dx, mc_dy, mc_dz), (z_min_km, z_max_km))
    ax.plot(tx, tz, "m--", lw=1.5, zorder=5, alpha=0.9, label="MC truth")
    ax.scatter([inj_pos_m[0]/1000.], [inj_pos_m[2]/1000.],
               marker="x", s=120, c="magenta", linewidths=2, zorder=6)

    # IterMPE track (green solid) — anchor at charge-weighted centroid of hit DOMs
    if len(dom_x) >= 4:
        ws  = nhits.astype(float)
        cx  = np.average(dom_x, weights=ws)
        cy  = np.average(dom_y, weights=ws)
        cz  = np.average(dom_z, weights=ws)
        tx2, tz2 = track_xz([cx, cy, cz], (im_dx, im_dy, im_dz), (z_min_km, z_max_km))
        ax.plot(tx2, tz2, color="limegreen", lw=2.0, zorder=6, alpha=0.9, label=f"IterMPE ({err:.2f}°)")

    # DM-Ice markers
    for dname, dpos in DMICE_POS.items():
        color = "deepskyblue" if dname == "det1" else "orange"
        ax.scatter([dpos[0]/1000.], [dpos[2]/1000.], marker="*", s=180,
                   c=color, zorder=7, edgecolors="k", linewidths=0.5, label=dname)

    ax.legend(loc="upper right", fontsize=8, framealpha=0.7)
    ax.set_xlabel("x [km]")
    ax.set_ylabel("z [km]")
    ax.set_xlim(-0.6, 0.6)
    ax.set_ylim(z_min_km, z_max_km)

    if len(t_us):
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar = fig.colorbar(sm, ax=ax, orientation="horizontal", fraction=0.05, pad=0.08)
        cbar.set_label(r"time [$\mu$s]")

    out_path = os.path.join(args.out, f"itermpe_ev{ev_idx:04d}.png")
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out_path}")

print("Done.")
