#!/usr/bin/env python3
"""
npz_to_linefit_csv.py

Runs IC LineFit and Pivot LineFit on every event in a BLO npz file and
writes a CSV compatible with plot_ang_err_vs_ndoms_energy.py.

Output columns:
    mc_energy_GeV, n_doms, n_hits, bin_id,
    ic_analytic_ang_err_deg, cfit_iter_ang_err_deg

Usage:
    python npz_to_linefit_csv.py [--npz PATH] [--out PATH]
"""

import os
import sys
import math
import argparse
import csv
import numpy as np

DEFAULT_NPZ = os.path.expanduser("~/dmice_work/output/muons_binned_200ev.npz")
DEFAULT_OUT = os.path.expanduser("~/dmice_work/output/200bin_simplots/linefit_results.csv")

DMICE_POS = {
    "det1": (31.25,   -72.93,  -2459.12),
    "det2": (-334.80, -424.50, -2459.33),
}
C_M_NS = 0.2998

# ── LineFit (same implementation as plot_event_display.py / run_linefit.py) ───

def _wm(vals, ws):
    W = sum(ws)
    return sum(v * w for v, w in zip(vals, ws)) / W if W else 0.0

def ic_linefit(xs, ys, zs, ts, ws):
    if len(xs) < 4:
        return None
    cx, cy, cz = _wm(xs, ws), _wm(ys, ws), _wm(zs, ws)
    tb  = _wm(ts, ws)
    dts = [t - tb for t in ts]
    den = sum(w * dt * dt for w, dt in zip(ws, dts))
    if not den:
        return None
    vx = sum(w * dt * (x - cx) for w, dt, x in zip(ws, dts, xs)) / den
    vy = sum(w * dt * (y - cy) for w, dt, y in zip(ws, dts, ys)) / den
    vz = sum(w * dt * (z - cz) for w, dt, z in zip(ws, dts, zs)) / den
    spd = math.sqrt(vx*vx + vy*vy + vz*vz)
    return (vx/spd, vy/spd, vz/spd) if spd else None

def pivot_linefit(xs, ys, zs, ts, ws, dm_pos, mc_dir):
    cx, cy, cz = _wm(xs, ws), _wm(ys, ws), _wm(zs, ws)
    tb = _wm(ts, ws)
    d  = (dm_pos[0]-cx, dm_pos[1]-cy, dm_pos[2]-cz)
    d_proj = d[0]*mc_dir[0] + d[1]*mc_dir[1] + d[2]*mc_dir[2]
    t_dm   = tb + d_proj / C_M_NS
    dts  = [t - t_dm for t in ts]
    drxs = [x - dm_pos[0] for x in xs]
    drys = [y - dm_pos[1] for y in ys]
    drzs = [z - dm_pos[2] for z in zs]
    den  = sum(w * dt * dt for w, dt in zip(ws, dts))
    if not den:
        return None
    vx = sum(w * dt * dr for w, dt, dr in zip(ws, dts, drxs)) / den
    vy = sum(w * dt * dr for w, dt, dr in zip(ws, dts, drys)) / den
    vz = sum(w * dt * dr for w, dt, dr in zip(ws, dts, drzs)) / den
    spd = math.sqrt(vx*vx + vy*vy + vz*vz)
    return (vx/spd, vy/spd, vz/spd) if spd else None

def ang_err_deg(truth, reco):
    dot = max(-1.0, min(1.0,
              truth[0]*reco[0] + truth[1]*reco[1] + truth[2]*reco[2]))
    return math.degrees(math.acos(abs(dot)))

# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--npz", default=DEFAULT_NPZ)
    parser.add_argument("--out", default=DEFAULT_OUT)
    args = parser.parse_args()

    if not os.path.exists(args.npz):
        sys.exit(f"ERROR: npz not found: {args.npz}")
    os.makedirs(os.path.dirname(args.out), exist_ok=True)

    print(f"Loading: {args.npz}")
    d = np.load(args.npz, allow_pickle=True)
    n = len(d["energy_GeV"])
    print(f"  {n} events")

    rows = []
    n_skip = 0

    for i in range(n):
        zen = float(d["zenith_rad"][i])
        azi = float(d["azimuth_rad"][i])
        mc_dir = (
            math.sin(zen) * math.cos(azi),
            math.sin(zen) * math.sin(azi),
            math.cos(zen),
        )

        dom_x  = np.array(d["dom_x"][i])
        dom_y  = np.array(d["dom_y"][i])
        dom_z  = np.array(d["dom_z"][i])
        dom_t  = np.array(d["dom_t"][i])
        nhits  = np.array(d["dom_nhits"][i])

        if len(dom_x) < 4:
            n_skip += 1
            continue

        xs, ys, zs, ts, ws = (dom_x.tolist(), dom_y.tolist(), dom_z.tolist(),
                               dom_t.tolist(), nhits.tolist())

        ic_dir = ic_linefit(xs, ys, zs, ts, ws)
        if ic_dir is None:
            n_skip += 1
            continue

        tgt_id = int(d["target_det"][i]) if "target_det" in d else 0
        dm_key = "det1" if tgt_id == 0 else "det2"
        piv_dir = pivot_linefit(xs, ys, zs, ts, ws, DMICE_POS[dm_key], mc_dir)
        if piv_dir is None:
            n_skip += 1
            continue

        rows.append({
            "mc_energy_GeV":          float(d["energy_GeV"][i]),
            "n_doms":                  int(d["n_doms"][i]),
            "n_hits":                  int(d["n_hits"][i]),
            "bin_id":                  int(d["bin_id"][i]) if "bin_id" in d else -1,
            "ic_analytic_ang_err_deg": ang_err_deg(mc_dir, ic_dir),
            "cfit_iter_ang_err_deg":   ang_err_deg(mc_dir, piv_dir),
        })

    print(f"  {len(rows)} events with valid fits  ({n_skip} skipped)")

    with open(args.out, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)

    print(f"Saved: {args.out}")


if __name__ == "__main__":
    main()
