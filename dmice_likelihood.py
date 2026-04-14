#!/usr/bin/env python3
"""
dmice_likelihood.py

DM-Ice NaI photon detection likelihood scorer.

Given a muon track hypothesis and DM-Ice hit times, computes:
    log L_DM-Ice = log p(Δt | d_perp) + log ε(d_perp)

where:
    Δt    = t_observed - t_geometric(track, DM-Ice position)
    d_perp = perpendicular distance from track to DM-Ice detector
    p(Δt | d_perp) = Gaussian(μ(d_perp), σ(d_perp))  [from BLO sim fit]
    ε(d_perp)      = detection efficiency              [exponential decay fit]

The model parameters are loaded from dmice_timing_model.npz (output of
build_dmice_timing_model.py).

Standalone usage:
    from dmice_likelihood import DMIceLikelihood
    model = DMIceLikelihood("~/dmice_work/output/dmice_timing_model.npz")
    log_l = model.score_track(track_pos, track_dir, t0, dm_hit_time, det_id=0)

IceTray usage (as a scoring module on existing fits):
    python3 dmice_likelihood.py --npz SIM_NPZ --model MODEL_NPZ --track MPEFit
"""

import os
import math
import argparse
import numpy as np

C_M_NS  = 0.2998
N_ICE   = 1.3195
THETA_C = math.acos(1.0 / N_ICE)   # ≈ 40.8°

# DM-Ice detector positions (BLO coordinates, metres)
DMICE_POS_BLO = {
    0: np.array([ 31.25,  -72.93, -2459.12]),   # det1
    1: np.array([-334.80, -424.50, -2459.33]),   # det2
}
# DM-Ice OMKeys in the geofile (string, sensor)
DMICE_OMKEYS = {0: (87, 1), 1: (88, 1)}


# ── Model class ───────────────────────────────────────────────────────────────

class DMIceLikelihood:
    """
    Loads the timing model and evaluates log L_DM-Ice for a track hypothesis.

    Only applied when d_perp < d_perp_max (through-detector events).
    Model: single Gaussian N(μ, σ) with scalar efficiency ε₀.
    """

    def __init__(self, model_path):
        model_path = os.path.expanduser(model_path)
        m = np.load(model_path, allow_pickle=True)
        self.mu         = float(m["mu_ns"])
        self.sigma      = float(m["sigma_ns"])
        self.eps0       = float(m["efficiency"])
        self.d_perp_max = float(m["d_perp_max_m"])
        print(f"[DMIceLikelihood] loaded model from {model_path}")
        print(f"  μ  = {self.mu:+.1f} ns")
        print(f"  σ  = {self.sigma:.1f} ns")
        print(f"  ε₀ = {self.eps0:.3f}")
        print(f"  d⊥ cut = {self.d_perp_max:.0f} m  (DM-Ice term ignored beyond this)")

    def t_geometric(self, track_pos, track_dir, t0, dm_pos):
        """
        Expected first-photon arrival time at dm_pos for a muon track.
        track_pos, track_dir: 3-vectors (BLO or IceCube coords — must match dm_pos)
        t0: muon time at track_pos [ns]
        """
        r      = np.asarray(dm_pos) - np.asarray(track_pos)
        d_hat  = np.asarray(track_dir, dtype=float)
        d_hat  = d_hat / np.linalg.norm(d_hat)
        s      = np.dot(r, d_hat)                          # along-track distance
        d_perp = math.sqrt(max(0.0, np.dot(r, r) - s**2))  # perpendicular distance
        t_pca  = t0 + s / C_M_NS
        if d_perp < 0.01:
            t_geo = t_pca
        else:
            t_geo = t_pca + d_perp / (C_M_NS * math.sin(THETA_C))
        return t_geo, d_perp

    def score_track(self, track_pos, track_dir, t0, dm_hit_time, det_id=0,
                    d_perp_true=None):
        """
        Compute log L_DM-Ice for a single track hypothesis and DM-Ice hit.

        Parameters
        ----------
        track_pos   : (x, y, z) reference point on track [m, BLO coords]
        track_dir   : (dx, dy, dz) unit vector (travel direction)
        t0          : muon time at track_pos [ns]
        dm_hit_time : observed first-photon time at DM-Ice [ns]
        det_id      : 0=det1, 1=det2
        d_perp_true : if provided, override the d_perp computed from geometry
                      (use for MC truth where track passes exactly through DM-Ice)

        Returns
        -------
        log_l  : float (log likelihood; nan if d_perp > cut)
        delta_t: float (timing residual = t_obs - t_geo) [ns]
        d_perp : float [m]
        """
        dm_pos  = DMICE_POS_BLO[det_id]
        t_geo, d_perp_geom = self.t_geometric(track_pos, track_dir, t0, dm_pos)
        d_perp  = d_perp_true if d_perp_true is not None else d_perp_geom
        delta_t = dm_hit_time - t_geo

        # Only apply DM-Ice term for through-detector tracks
        if d_perp > self.d_perp_max:
            return float("nan"), delta_t, d_perp

        # Gaussian log-likelihood on timing residual
        log_l_timing = (-0.5 * ((delta_t - self.mu) / self.sigma)**2
                        - math.log(self.sigma * math.sqrt(2 * math.pi)))

        # Efficiency term (Bernoulli: we observed a hit)
        log_l_eff = math.log(max(self.eps0, 1e-9))

        return log_l_timing + log_l_eff, delta_t, d_perp

    def score_no_hit(self, track_pos, track_dir, t0, det_id=0):
        """
        Log likelihood when DM-Ice did NOT fire (absence of hit).
        Only meaningful for through-detector tracks (d_perp < d_perp_max).
        """
        dm_pos = DMICE_POS_BLO[det_id]
        _, d_perp = self.t_geometric(track_pos, track_dir, t0, dm_pos)
        if d_perp > self.d_perp_max:
            return 0.0, d_perp   # no information from non-through tracks
        return math.log(max(1 - self.eps0, 1e-9)), d_perp


# ── Standalone scoring on BLO sim ────────────────────────────────────────────

def load_ragged(d, key, N):
    if f"{key}_flat" in d:
        flat    = d[f"{key}_flat"]
        offsets = d[f"{key}_offsets"]
        return [flat[offsets[i]:offsets[i+1]] for i in range(N)]
    return list(d[key])


def score_sim_npz(npz_path, model, track_key="linefit"):
    """
    Score every event in a BLO sim npz using IC LineFit (or MC truth) as the
    track hypothesis, then compare log L between:
      - MC truth track
      - IC LineFit track (Python pivot fit)
      - Pivot LineFit track

    Returns a list of dicts with per-event scores.
    """
    d = np.load(npz_path, allow_pickle=True)
    N = len(d["energy_GeV"])

    xs_all  = load_ragged(d, "dom_x",      N)
    ys_all  = load_ragged(d, "dom_y",      N)
    zs_all  = load_ragged(d, "dom_z",      N)
    ts_all  = load_ragged(d, "dom_t",      N)
    ws_all  = load_ragged(d, "dom_nhits",  N)
    str_all = load_ragged(d, "dom_string", N)
    sen_all = load_ragged(d, "dom_sensor", N)

    def linefit_dir(xs, ys, zs, ts, ws):
        """Analytic IC LineFit (weighted least squares)."""
        W  = ws.sum()
        if W == 0:
            return None, None, None
        cx = np.dot(ws, xs) / W
        cy = np.dot(ws, ys) / W
        cz = np.dot(ws, zs) / W
        tb = np.dot(ws, ts) / W
        dts = ts - tb
        den = np.dot(ws, dts**2)
        if den == 0:
            return None, None, None
        vx = np.dot(ws * dts, xs - cx) / den
        vy = np.dot(ws * dts, ys - cy) / den
        vz = np.dot(ws * dts, zs - cz) / den
        spd = math.sqrt(vx**2 + vy**2 + vz**2)
        if spd < 1e-6:
            return None, None, None
        return (vx/spd, vy/spd, vz/spd), np.array([cx, cy, cz]), tb

    def pivot_dir(xs, ys, zs, ts, ws, seed_dir, dm_pos):
        """Pivot LineFit anchored to DM-Ice transit time."""
        W  = ws.sum()
        if W == 0:
            return None
        cx = np.dot(ws, xs) / W
        cy = np.dot(ws, ys) / W
        cz = np.dot(ws, zs) / W
        tb = np.dot(ws, ts) / W
        d_proj = ((dm_pos[0]-cx)*seed_dir[0] + (dm_pos[1]-cy)*seed_dir[1]
                  + (dm_pos[2]-cz)*seed_dir[2])
        t_dm = tb + d_proj / C_M_NS
        dts  = ts - t_dm
        drxs = xs - dm_pos[0]
        drys = ys - dm_pos[1]
        drzs = zs - dm_pos[2]
        den  = np.dot(ws, dts**2)
        if den == 0:
            return None
        vx = np.dot(ws * dts, drxs) / den
        vy = np.dot(ws * dts, drys) / den
        vz = np.dot(ws * dts, drzs) / den
        spd = math.sqrt(vx**2 + vy**2 + vz**2)
        return (vx/spd, vy/spd, vz/spd) if spd > 1e-6 else None

    def ang_err(d1, d2):
        dot = max(-1.0, min(1.0, d1[0]*d2[0] + d1[1]*d2[1] + d1[2]*d2[2]))
        return math.degrees(math.acos(abs(dot)))

    rows = []
    for i in range(N):
        xs      = np.asarray(xs_all[i],  dtype=float)
        ys      = np.asarray(ys_all[i],  dtype=float)
        zs      = np.asarray(zs_all[i],  dtype=float)
        ts      = np.asarray(ts_all[i],  dtype=float)
        ws      = np.asarray(ws_all[i],  dtype=float)
        strings = np.asarray(str_all[i], dtype=int)
        sensors = np.asarray(sen_all[i], dtype=int)

        zen = float(d["zenith_rad"][i])
        azi = float(d["azimuth_rad"][i])
        mc_dir = (math.sin(zen)*math.cos(azi),
                  math.sin(zen)*math.sin(azi),
                  math.cos(zen))
        tgt_id = int(d["target_det"][i]) if "target_det" in d else 0
        dm_pos = DMICE_POS_BLO[tgt_id]
        s_dm, sen_dm = DMICE_OMKEYS[tgt_id]

        # DM-Ice hit
        dm_mask = (strings == s_dm) & (sensors == sen_dm)
        has_dm_hit = dm_mask.any()
        dm_t = float(ts[dm_mask].min()) if has_dm_hit else float("nan")

        # IC-only hits for LineFit
        ic_mask = ~dm_mask
        if ic_mask.sum() < 4:
            continue

        lf_result = linefit_dir(xs[ic_mask], ys[ic_mask], zs[ic_mask],
                                ts[ic_mask], ws[ic_mask])
        lf_dir, lf_pos, lf_t0 = lf_result

        piv_dir = None
        if lf_dir is not None:
            piv_dir = pivot_dir(xs[ic_mask], ys[ic_mask], zs[ic_mask],
                                ts[ic_mask], ws[ic_mask], lf_dir, dm_pos)

        row = dict(
            event_idx  = i,
            energy     = float(d["energy_GeV"][i]),
            zenith_deg = math.degrees(zen),
            has_dm_hit = has_dm_hit,
            tgt_id     = tgt_id,
            # Angular errors vs MC truth
            lf_ang_err_deg    = ang_err(mc_dir, lf_dir)   if lf_dir  else float("nan"),
            pivot_ang_err_deg = ang_err(mc_dir, piv_dir)  if piv_dir else float("nan"),
        )

        if has_dm_hit and lf_dir is not None:
            # Score MC truth track — d_perp=0 by construction (track aimed through DM-Ice)
            ll_mc, dt_mc, dp_mc = model.score_track(
                lf_pos, mc_dir, lf_t0, dm_t, tgt_id, d_perp_true=0.0)
            # Score LineFit track — use computed d_perp (reflects how much LF misses DM-Ice)
            ll_lf, dt_lf, dp_lf = model.score_track(
                lf_pos, lf_dir, lf_t0, dm_t, tgt_id)
            # Score Pivot LineFit — anchored to DM-Ice, so d_perp≈0
            ll_piv = float("nan")
            if piv_dir:
                ll_piv, _, _ = model.score_track(
                    lf_pos, piv_dir, lf_t0, dm_t, tgt_id, d_perp_true=0.0)

            row.update(dict(
                d_perp_m       = dp_lf,
                delta_t_mc_ns  = dt_mc,
                delta_t_lf_ns  = dt_lf,
                log_l_mc       = ll_mc,
                log_l_lf       = ll_lf,
                log_l_pivot    = ll_piv,
                log_l_gain     = ll_piv - ll_lf,   # positive → pivot preferred
            ))
        rows.append(row)

    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", default=os.path.expanduser(
        "~/dmice_work/output/dmice_timing_model.npz"))
    parser.add_argument("--npz", default=os.path.expanduser(
        "~/dmice_work/output/muons_binned_200ev.npz"))
    parser.add_argument("--out", default=os.path.expanduser(
        "~/dmice_work/output"))
    args = parser.parse_args()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    model = DMIceLikelihood(args.model)
    print(f"\nScoring events from {os.path.basename(args.npz)} ...")
    rows = score_sim_npz(args.npz, model)

    valid = [r for r in rows if "log_l_lf" in r]
    print(f"Events scored: {len(valid)} (of {len(rows)} with ≥4 IC hits)")

    ll_mc  = np.array([r["log_l_mc"]       for r in valid])
    ll_lf  = np.array([r["log_l_lf"]       for r in valid])
    ll_piv = np.array([r["log_l_pivot"]     for r in valid])
    gain   = np.array([r["log_l_gain"]      for r in valid])
    lf_err = np.array([r["lf_ang_err_deg"]  for r in valid])
    pv_err = np.array([r["pivot_ang_err_deg"] for r in valid])
    dp     = np.array([r["d_perp_m"]        for r in valid])
    en     = np.array([r["energy"]          for r in valid])

    print(f"\nMedian log L (MC truth):    {np.median(ll_mc):.2f}")
    print(f"Median log L (LineFit):     {np.median(ll_lf):.2f}")
    print(f"Median log L (Pivot LF):    {np.nanmedian(ll_piv):.2f}")
    print(f"Median log L gain (pivot-LF): {np.median(gain):.2f}")
    print(f"\nAngular errors:")
    print(f"  LineFit median:    {np.nanmedian(lf_err):.2f}°")
    print(f"  Pivot LF median:   {np.nanmedian(pv_err):.2f}°")
    piv_better = np.nansum(pv_err < lf_err)
    print(f"  Pivot better:      {piv_better}/{len(valid)} ({100*piv_better/len(valid):.0f}%)")

    # ── Plots ─────────────────────────────────────────────────────────────────
    fig, axes = plt.subplots(2, 3, figsize=(16, 9))

    # 1. Log L distribution: MC vs LF vs Pivot
    ax = axes[0, 0]
    bins_ll = np.linspace(np.percentile(ll_lf, 2), np.percentile(ll_mc, 98), 50)
    ax.hist(ll_mc,  bins=bins_ll, histtype="step", lw=2, color="seagreen",
            density=True, label=f"MC truth  (med={np.median(ll_mc):.1f})")
    ax.hist(ll_lf,  bins=bins_ll, histtype="step", lw=2, color="steelblue",
            density=True, label=f"LineFit   (med={np.median(ll_lf):.1f})")
    ax.hist(ll_piv[~np.isnan(ll_piv)], bins=bins_ll, histtype="step", lw=2,
            color="crimson", density=True,
            label=f"Pivot LF  (med={np.nanmedian(ll_piv):.1f})")
    ax.set_xlabel("log L_DM-Ice")
    ax.set_ylabel("Normalised / bin")
    ax.set_title("DM-Ice likelihood: MC truth vs reconstructions")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # 2. log L gain (pivot - LF) vs angular improvement
    ax = axes[0, 1]
    ang_gain = lf_err - pv_err   # positive = pivot improved
    sc = ax.scatter(gain, ang_gain, c=np.log10(np.maximum(en, 1)),
                    s=12, alpha=0.5, cmap="viridis")
    fig.colorbar(sc, ax=ax, label="log₁₀(E/GeV)")
    ax.axhline(0, color="k", lw=0.8, ls="--")
    ax.axvline(0, color="k", lw=0.8, ls="--")
    ax.set_xlabel("log L gain (pivot − LF)")
    ax.set_ylabel("Angular improvement: LF err − Pivot err (°)")
    ax.set_title("DM-Ice likelihood gain vs angular improvement")
    ax.grid(True, alpha=0.3)

    # 3. log L vs angular error (LF)
    ax = axes[0, 2]
    ax.scatter(ll_lf, lf_err, c="steelblue", s=10, alpha=0.5, label="LineFit")
    ax.scatter(ll_piv[~np.isnan(ll_piv)], pv_err[~np.isnan(ll_piv)],
               c="crimson", s=10, alpha=0.5, label="Pivot LF")
    ax.set_xlabel("log L_DM-Ice")
    ax.set_ylabel("Angular error vs MC truth (°)")
    ax.set_title("DM-Ice likelihood vs reconstruction accuracy")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # 4. log L gain vs d_perp
    ax = axes[1, 0]
    ax.scatter(dp, gain, c=np.log10(np.maximum(en, 1)),
               s=12, alpha=0.5, cmap="viridis")
    ax.axhline(0, color="k", lw=0.8, ls="--")
    ax.set_xlabel("d_perp (m)")
    ax.set_ylabel("log L gain (pivot − LF)")
    ax.set_title("Likelihood gain vs closest approach")
    ax.grid(True, alpha=0.3)

    # 5. Angular error: LF vs Pivot scatter
    ax = axes[1, 1]
    lim = max(np.nanpercentile(lf_err, 98), np.nanpercentile(pv_err, 98))
    ax.scatter(lf_err, pv_err, c=gain, s=12, alpha=0.5,
               cmap="RdYlGn", vmin=-2, vmax=2)
    ax.plot([0, lim], [0, lim], "k--", lw=0.8)
    ax.set_xlabel("LineFit angular error (°)")
    ax.set_ylabel("Pivot LF angular error (°)")
    ax.set_title("Pivot vs LineFit accuracy\n(green=Pivot preferred by DM-Ice L)")
    ax.set_xlim(0, lim); ax.set_ylim(0, lim)
    ax.grid(True, alpha=0.3)

    # 6. Δt(MC truth) vs Δt(LF) — model discrimination
    ax = axes[1, 2]
    dt_mc_arr = np.array([r["delta_t_mc_ns"] for r in valid])
    dt_lf_arr = np.array([r["delta_t_lf_ns"] for r in valid])
    ax.scatter(dt_mc_arr, dt_lf_arr, c=lf_err, s=12, alpha=0.5, cmap="hot_r")
    ax.plot([dt_mc_arr.min(), dt_mc_arr.max()],
            [dt_mc_arr.min(), dt_mc_arr.max()], "k--", lw=0.8)
    ax.set_xlabel("Δt (MC truth track) [ns]")
    ax.set_ylabel("Δt (LineFit track) [ns]")
    ax.set_title("Timing residual: truth vs fit\n(coloured by LineFit ang error)")
    ax.grid(True, alpha=0.3)

    fig.suptitle(
        "DM-Ice likelihood scoring — BLO sim with MC truth\n"
        f"({len(valid)} events with DM-Ice hit and valid LineFit)",
        fontsize=12
    )
    plt.tight_layout()

    out_path = os.path.join(args.out, "dmice_likelihood_scoring.png")
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"\nPlot saved: {out_path}")


if __name__ == "__main__":
    main()
