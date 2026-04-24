#!/usr/bin/env python3
"""
compare_splinempe_seeds.py

Compare SplineMPE seeded two ways on BLO simulation:
  A) Standard seed: existing LineFit (no DM-Ice information)
  B) Pivot seed:    LineFit anchored to DM-Ice transit time

For each event:
  1. Compute Pivot LineFit direction from IC pulses + DM-Ice hit time
  2. Run SplineMPE with seed A and seed B
  3. Score both SplineMPE results against DM-Ice log L
  4. Report angular errors vs MC truth

Ground truth is available because input is BLO simulation.

Run on Cobalt (requires IceTray + SplineMPE tables):
  /cvmfs/.../icetray/v1.12.1/env-shell.sh python3 ~/dmice/compare_splinempe_seeds.py \
    --i3 ~/dmice_work/output/blo_dmice_targeted_det1det2_both_1000events.i3.zst \
    --out ~/dmice_work/output/splinempe_seed_comparison.csv

SplineMPE tables (SPICEMie, used for BLO sim):
  /cvmfs/icecube.opensciencegrid.org/data/photon-tables/splines/ems_mie_z20_a10.abs.fits
  /cvmfs/icecube.opensciencegrid.org/data/photon-tables/splines/ems_mie_z20_a10.prob.fits
"""

import os, math, argparse
import numpy as np

# ── Constants ─────────────────────────────────────────────────────────────────
C_M_NS  = 0.2998
N_ICE   = 1.3195
THETA_C = math.acos(1.0 / N_ICE)

DMICE_POS_IC = {
    0: np.array([ 31.25,  -72.93, -511.05]),
    1: np.array([-334.80, -424.50, -511.26]),
}
DMICE_OMKEYS = {0: (87, 1), 1: (88, 1)}

SPLINE_DIR = "/cvmfs/icecube.opensciencegrid.org/data/photon-tables/splines"
ABS_TABLE  = os.path.join(SPLINE_DIR, "ems_mie_z20_a10.abs.fits")
PROB_TABLE = os.path.join(SPLINE_DIR, "ems_mie_z20_a10.prob.fits")

# ── Geometry helpers ──────────────────────────────────────────────────────────

def t_geometric_ic(track_pos, track_dir, t0, dm_pos):
    r     = np.asarray(dm_pos) - np.asarray(track_pos)
    d_hat = np.asarray(track_dir, dtype=float)
    d_hat = d_hat / np.linalg.norm(d_hat)
    s     = np.dot(r, d_hat)
    d_perp = math.sqrt(max(0.0, np.dot(r, r) - s**2))
    t_pca  = t0 + s / C_M_NS
    t_geo  = t_pca if d_perp < 0.01 else t_pca + d_perp / (C_M_NS * math.sin(THETA_C))
    return t_geo, d_perp


def score_track(track_pos, track_dir, t0, dm_hit_time, dm_pos,
                mu, sigma, eps0, d_perp_max, d_perp_override=None):
    t_geo, d_perp_geom = t_geometric_ic(track_pos, track_dir, t0, dm_pos)
    d_perp = d_perp_override if d_perp_override is not None else d_perp_geom
    delta_t = dm_hit_time - t_geo
    if d_perp > d_perp_max:
        return float("nan"), delta_t, d_perp
    log_l_t = (-0.5 * ((delta_t - mu) / sigma)**2
               - math.log(sigma * math.sqrt(2 * math.pi)))
    return log_l_t + math.log(max(eps0, 1e-9)), delta_t, d_perp


def ang_err_deg(d1, d2):
    dot = max(-1.0, min(1.0, float(np.dot(np.asarray(d1), np.asarray(d2)))))
    return math.degrees(math.acos(abs(dot)))


def compute_pivot_lf(xs, ys, zs, ts, ws, seed_dir, dm_pos):
    """Pivot LineFit: anchor centroid time to DM-Ice transit, refit direction."""
    ws = np.asarray(ws, dtype=float)
    W  = ws.sum()
    if W == 0:
        return None
    cx = np.dot(ws, xs) / W
    cy = np.dot(ws, ys) / W
    cz = np.dot(ws, zs) / W
    tb = np.dot(ws, ts) / W
    # Project DM-Ice onto seed direction from centroid
    d_proj = ((dm_pos[0]-cx)*seed_dir[0] + (dm_pos[1]-cy)*seed_dir[1]
              + (dm_pos[2]-cz)*seed_dir[2])
    t_dm = tb + d_proj / C_M_NS
    # Refit relative to DM-Ice position + time
    dts  = ts - t_dm
    drxs = xs - dm_pos[0]
    drys = ys - dm_pos[1]
    drzs = zs - dm_pos[2]
    den  = np.dot(ws, dts**2)
    if den < 1e-10:
        return None
    vx = np.dot(ws * dts, drxs) / den
    vy = np.dot(ws * dts, drys) / den
    vz = np.dot(ws * dts, drzs) / den
    spd = math.sqrt(vx**2 + vy**2 + vz**2)
    if spd < 1e-6:
        return None
    return (vx/spd, vy/spd, vz/spd)


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--i3", default=os.path.expanduser(
        "~/dmice_work/output/blo_dmice_targeted_det1det2_both_1000events.i3.zst"))
    parser.add_argument("--model", default=os.path.expanduser(
        "~/dmice_work/output/dmice_timing_model.npz"))
    parser.add_argument("--out", default=os.path.expanduser(
        "~/dmice_work/output/splinempe_seed_comparison.csv"))
    parser.add_argument("--pulse-key", default="InIcePulses")
    args = parser.parse_args()

    # ── Load timing model ─────────────────────────────────────────────────
    m     = np.load(args.model, allow_pickle=True)
    mu    = float(m["mu_ns"])
    sigma = float(m["sigma_ns"])
    eps0  = float(m["efficiency"])
    d_max = float(m["d_perp_max_m"])
    print(f"Model: μ={mu:+.1f} ns  σ={sigma:.1f} ns  ε={eps0:.3f}  d⊥_max={d_max:.0f} m")

    # ── IceTray setup ─────────────────────────────────────────────────────
    from icecube import icetray, dataio, dataclasses, recclasses
    from icecube.icetray import I3Units

    # SplineMPE
    try:
        from icecube import spline_reco
        spline_abs  = spline_reco.I3SplineRecoAmplitudeTableService(ABS_TABLE)
        spline_prob = spline_reco.I3SplineRecoProbabilityTableService(PROB_TABLE)
        has_spline  = True
        print(f"SplineMPE loaded: {ABS_TABLE}")
    except Exception as e:
        has_spline = False
        print(f"SplineMPE unavailable: {e}")

    # Geometry
    try:
        from icecube import phys_services
    except ImportError:
        pass

    rows = []
    n_total = n_dm_hit = 0

    f = dataio.I3File(args.i3)
    geo = None

    while f.more():
        frame = f.pop_frame()

        if frame.Stop == icetray.I3Frame.Geometry:
            geo = frame["I3Geometry"]
            continue
        if frame.Stop != icetray.I3Frame.Physics:
            continue

        n_total += 1

        # ── MC truth ─────────────────────────────────────────────────────
        mctree  = frame.Get("I3MCTree") if "I3MCTree" in frame else None
        if mctree is None:
            continue
        primary = mctree.get_primaries()[0]
        mc_dir  = np.array([primary.dir.x, primary.dir.y, primary.dir.z])
        mc_pos  = np.array([primary.pos.x, primary.pos.y, primary.pos.z])
        mc_t0   = primary.time

        # ── Target DM-Ice detector ────────────────────────────────────────
        blo_det = str(frame["BLO_DetId"]) if "BLO_DetId" in frame else ""
        tgt_id  = 0 if "det1" in blo_det else 1
        dm_pos  = DMICE_POS_IC[tgt_id]
        dm_omkey_tup = DMICE_OMKEYS[tgt_id]
        dm_omkey = icetray.OMKey(dm_omkey_tup[0], dm_omkey_tup[1])

        # ── Extract pulses ────────────────────────────────────────────────
        if args.pulse_key not in frame:
            continue
        try:
            pulses_raw = dataclasses.I3RecoPulseSeriesMap.from_frame(frame, args.pulse_key)
        except Exception:
            pulses_raw = frame[args.pulse_key]

        xs, ys, zs, ts, ws = [], [], [], [], []
        dm_t = None
        for omkey, ps in pulses_raw:
            if not ps:
                continue
            if (omkey.string, omkey.om) == dm_omkey_tup:
                dm_t = min(p.time for p in ps)
                continue
            if geo is None:
                continue
            omgeo = geo.omgeo[omkey]
            pos   = omgeo.position
            nhits = sum(p.charge for p in ps)
            xs.append(pos.x); ys.append(pos.y); zs.append(pos.z)
            ts.append(min(p.time for p in ps)); ws.append(nhits)

        if dm_t is None or len(xs) < 4:
            continue
        n_dm_hit += 1

        xs = np.array(xs); ys = np.array(ys); zs = np.array(zs)
        ts = np.array(ts); ws = np.array(ws, dtype=float)

        # ── LineFit (standard seed) ───────────────────────────────────────
        if "LineFit" in frame:
            lf = frame["LineFit"]
            lf_pos = np.array([lf.pos.x, lf.pos.y, lf.pos.z])
            lf_dir = np.array([lf.dir.x, lf.dir.y, lf.dir.z])
            lf_t0  = lf.time
        else:
            # Compute LineFit inline
            W  = ws.sum()
            cx = np.dot(ws, xs)/W; cy = np.dot(ws, ys)/W; cz = np.dot(ws, zs)/W
            tb = np.dot(ws, ts)/W
            dts_lf = ts - tb; den = np.dot(ws, dts_lf**2)
            if den < 1e-10:
                continue
            vx = np.dot(ws*dts_lf, xs-cx)/den
            vy = np.dot(ws*dts_lf, ys-cy)/den
            vz = np.dot(ws*dts_lf, zs-cz)/den
            spd = math.sqrt(vx**2+vy**2+vz**2)
            if spd < 1e-6:
                continue
            lf_dir = np.array([vx/spd, vy/spd, vz/spd])
            lf_pos = np.array([cx, cy, cz]); lf_t0 = tb

        # ── Pivot LineFit (DM-Ice seeded) ─────────────────────────────────
        piv_dir_tup = compute_pivot_lf(xs, ys, zs, ts, ws, lf_dir, dm_pos)
        if piv_dir_tup is None:
            piv_dir = lf_dir  # fallback
        else:
            piv_dir = np.array(piv_dir_tup)

        # ── Existing MPEFit (standard full reco) ──────────────────────────
        # NOTE: SplineMPE with pivot seed is Phase D (requires IceTray tray).
        # For now we compare LineFit / Pivot LF / MPEFit against DM-Ice log L.
        mpe_dir = mpe_pos = mpe_t0 = None
        if "MPEFit" in frame:
            mpe = frame["MPEFit"]
            if mpe.fit_status == dataclasses.I3Particle.OK:
                mpe_dir = np.array([mpe.dir.x, mpe.dir.y, mpe.dir.z])
                mpe_pos = np.array([mpe.pos.x, mpe.pos.y, mpe.pos.z])
                mpe_t0  = mpe.time

        spe_std_dir = spe_std_pos = spe_std_t0 = None
        spe_piv_dir = spe_piv_pos = spe_piv_t0 = None

        # ── DM-Ice log L scoring ──────────────────────────────────────────
        def sc(pos, dirv, t0, d_override=None):
            if pos is None:
                return float("nan"), float("nan"), float("nan")
            return score_track(pos, dirv, t0, dm_t, dm_pos,
                               mu, sigma, eps0, d_max, d_override)

        ll_mc,  dt_mc,  dp_mc  = sc(mc_pos,      mc_dir,      mc_t0,  0.0)
        ll_lf,  dt_lf,  dp_lf  = sc(lf_pos,      lf_dir,      lf_t0)
        ll_piv, dt_piv, dp_piv = sc(lf_pos,       piv_dir,     lf_t0,  0.0)
        ll_spe_std, _, _        = sc(spe_std_pos, spe_std_dir, spe_std_t0)
        ll_spe_piv, _, _        = sc(spe_piv_pos, spe_piv_dir, spe_piv_t0)
        ll_mpe, _, _            = sc(mpe_pos,     mpe_dir,     mpe_t0)

        # Winner = track with highest DM-Ice log L
        candidates = {
            "lf":      (ll_lf,      lf_dir,      ang_err_deg(mc_dir, lf_dir)),
            "pivot":   (ll_piv,     piv_dir,     ang_err_deg(mc_dir, piv_dir)),
        }
        if spe_std_dir is not None:
            candidates["spe_std"] = (ll_spe_std, spe_std_dir,
                                      ang_err_deg(mc_dir, spe_std_dir))
        if spe_piv_dir is not None:
            candidates["spe_piv"] = (ll_spe_piv, spe_piv_dir,
                                      ang_err_deg(mc_dir, spe_piv_dir))
        valid = {k: v for k, v in candidates.items() if np.isfinite(v[0])}
        if valid:
            winner_key = max(valid, key=lambda k: valid[k][0])
            winner_ang = valid[winner_key][2]
        else:
            winner_key = "none"; winner_ang = float("nan")

        row = dict(
            tgt_id          = tgt_id,
            dp_lf_m         = dp_lf,
            # Angular errors
            lf_ang_err      = ang_err_deg(mc_dir, lf_dir),
            piv_ang_err     = ang_err_deg(mc_dir, piv_dir),
            mpe_ang_err     = ang_err_deg(mc_dir, mpe_dir) if mpe_dir is not None else float("nan"),
            spe_std_ang_err = ang_err_deg(mc_dir, spe_std_dir) if spe_std_dir is not None else float("nan"),
            spe_piv_ang_err = ang_err_deg(mc_dir, spe_piv_dir) if spe_piv_dir is not None else float("nan"),
            winner_key      = winner_key,
            winner_ang_err  = winner_ang,
            # Log likelihoods
            ll_mc           = ll_mc,
            ll_lf           = ll_lf,
            ll_piv          = ll_piv,
            ll_mpe          = ll_mpe,
            ll_spe_std      = ll_spe_std,
            ll_spe_piv      = ll_spe_piv,
            ll_gain_piv_lf  = ll_piv - ll_lf,
        )
        rows.append(row)

        if n_dm_hit % 50 == 0:
            print(f"  [{n_dm_hit}/{n_total}]  winner={winner_key}  "
                  f"piv_ang={ang_err_deg(mc_dir, piv_dir):.1f}°  "
                  f"lf_ang={ang_err_deg(mc_dir, lf_dir):.1f}°")

    f.close()
    print(f"\nDone. {n_total} events, {n_dm_hit} with DM-Ice hit → {len(rows)} rows.")

    import pandas as pd
    df = pd.DataFrame(rows)
    df.to_csv(args.out, index=False)
    print(f"CSV: {args.out}")

    # ── Summary ───────────────────────────────────────────────────────────
    print(f"\nAngular error medians:")
    for col, label in [("lf_ang_err", "LineFit"),
                        ("piv_ang_err", "Pivot LF"),
                        ("mpe_ang_err", "MPEFit"),
                        ("spe_std_ang_err", "SplineMPE(std seed)"),
                        ("spe_piv_ang_err", "SplineMPE(piv seed)"),
                        ("winner_ang_err", "DM-Ice selected winner")]:
        vals = df[col].dropna()
        if len(vals):
            print(f"  {label:30s}: {vals.median():.2f}°  (n={len(vals)})")

    # ── Plot ──────────────────────────────────────────────────────────────
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    bins = np.linspace(0, 15, 31)

    ax = axes[0]
    for col, label, color in [
        ("lf_ang_err",      "LineFit",              "steelblue"),
        ("piv_ang_err",     "Pivot LF",             "darkorange"),
        ("mpe_ang_err",     "MPEFit",               "green"),
        ("spe_std_ang_err", "SplineMPE (std)",      "purple"),
        ("spe_piv_ang_err", "SplineMPE (piv seed)", "red"),
    ]:
        vals = df[col].dropna()
        if len(vals) > 2:
            ax.hist(vals.clip(0, 15), bins=bins, histtype="step", lw=2,
                    label=f"{label} ({vals.median():.1f}°)", color=color)
    ax.set_xlabel("Angular error (°)"); ax.set_ylabel("Events")
    ax.set_title("All events"); ax.legend(fontsize=7); ax.grid(True, alpha=0.3)

    ax = axes[1]
    ax.scatter(df.lf_ang_err.clip(0,15), df.piv_ang_err.clip(0,15),
               c=df.ll_gain_piv_lf.fillna(0), cmap="RdYlGn", s=8, alpha=0.6)
    ax.plot([0,15],[0,15],"k--",lw=1,alpha=0.5)
    ax.set_xlabel("LineFit ang err (°)"); ax.set_ylabel("Pivot LF ang err (°)")
    ax.set_title("Pivot vs LineFit\n(colour = ΔlogL = logL_piv − logL_LF)")
    ax.grid(True, alpha=0.3)

    ax = axes[2]
    if df.spe_piv_ang_err.notna().sum() > 2:
        ax.scatter(df.spe_std_ang_err.clip(0,15), df.spe_piv_ang_err.clip(0,15),
                   c=df.ll_spe_piv.fillna(-999), cmap="viridis", s=8, alpha=0.6)
        ax.plot([0,15],[0,15],"k--",lw=1,alpha=0.5)
        ax.set_xlabel("SplineMPE(std seed) (°)"); ax.set_ylabel("SplineMPE(piv seed) (°)")
        ax.set_title("SplineMPE: pivot seed vs standard seed\n(colour = DM-Ice logL)")
    else:
        ax.text(0.5, 0.5, "SplineMPE not available", transform=ax.transAxes, ha="center")
    ax.grid(True, alpha=0.3)

    fig.suptitle(
        f"DM-Ice seed comparison on BLO sim  ({n_dm_hit} events with DM-Ice hit)\n"
        f"Model: μ={mu:+.0f}ns σ={sigma:.0f}ns  d⊥<{d_max:.0f}m",
        fontsize=10
    )
    plt.tight_layout()
    out_png = args.out.replace(".csv", ".png")
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Plot: {out_png}")


if __name__ == "__main__":
    main()
