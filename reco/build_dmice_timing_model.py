#!/usr/bin/env python3
"""
build_dmice_timing_model.py

Fits a DM-Ice NaI timing model using only events where the muon passes
DIRECTLY THROUGH the detector (d_perp < D_PERP_MAX, default 15 m).

Model: single Gaussian  p(Δt) = N(μ, σ)
       detection efficiency ε₀ (scalar)

Output:
  ~/dmice_work/output/dmice_timing_model.npz
  ~/dmice_work/output/dmice_timing_residuals.png

Run on WARD:
  python3 build_dmice_timing_model.py [--dperp-max 15] [--out DIR]
"""

import os, math, argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import norm as sp_norm

OUT_DEFAULT = os.path.expanduser("~/dmice_work/output")
D_PERP_MAX_DEFAULT = 15.0   # metres — "directly through" cut

NPZ_ON_AXIS = [
    os.path.expanduser("~/dmice_work/output/muons_binned_200ev.npz"),
    os.path.expanduser("~/dmice_work/output/muons_offset_0m_1000ev.npz"),
]

DMICE_DOMS = {
    (87, 1): np.array([ 31.25,  -72.93, -2459.12]),
    (88, 1): np.array([-334.80, -424.50, -2459.33]),
}
C_M_NS  = 0.2998
N_ICE   = 1.3195
THETA_C = math.acos(1.0 / N_ICE)


# ── Helpers ───────────────────────────────────────────────────────────────────

def load_ragged(d, key, N):
    if f"{key}_flat" in d:
        flat    = d[f"{key}_flat"]
        offsets = d[f"{key}_offsets"]
        return [flat[offsets[i]:offsets[i+1]] for i in range(N)]
    return list(d[key])


def t_geometric_dmice(xs, ys, zs, ts, ws, dm_pos, mc_dir):
    """
    Cherenkov-corrected estimate of muon transit time past dm_pos,
    using weighted IC DOM hits and MC truth direction.
    Returns (t_geo, d_perp).
    """
    ws = np.asarray(ws, dtype=float)
    W  = ws.sum()
    if W == 0:
        return None, None
    d_hat    = np.asarray(mc_dir, dtype=float)
    r_dm     = np.asarray(dm_pos, dtype=float)
    dom_pos  = np.column_stack([np.asarray(xs, float),
                                np.asarray(ys, float),
                                np.asarray(zs, float)])
    delta    = dom_pos - r_dm
    s_i      = delta @ d_hat
    dp_i     = np.sqrt(np.maximum(np.sum(delta**2, axis=1) - s_i**2, 0.0))
    t_corr   = np.asarray(ts, float) - s_i / C_M_NS - dp_i / (C_M_NS * math.sin(THETA_C))
    t_geo    = np.dot(ws, t_corr) / W

    # d_perp of DM-Ice from centroid-estimated track
    cx = np.dot(ws, dom_pos[:, 0]) / W
    cy = np.dot(ws, dom_pos[:, 1]) / W
    cz = np.dot(ws, dom_pos[:, 2]) / W
    delta_dm = r_dm - np.array([cx, cy, cz])
    s_dm     = np.dot(delta_dm, d_hat)
    d_perp   = np.linalg.norm(delta_dm - s_dm * d_hat)
    return t_geo, d_perp


def extract_hits(npz_path, d_perp_max):
    d  = np.load(npz_path, allow_pickle=True)
    N  = len(d["energy_GeV"])
    xs_all  = load_ragged(d, "dom_x",      N)
    ys_all  = load_ragged(d, "dom_y",      N)
    zs_all  = load_ragged(d, "dom_z",      N)
    ts_all  = load_ragged(d, "dom_t",      N)
    ws_all  = load_ragged(d, "dom_nhits",  N)
    str_all = load_ragged(d, "dom_string", N)
    sen_all = load_ragged(d, "dom_sensor", N)

    records = []
    for i in range(N):
        xs = np.asarray(xs_all[i],  float); ys = np.asarray(ys_all[i],  float)
        zs = np.asarray(zs_all[i],  float); ts = np.asarray(ts_all[i],  float)
        ws = np.asarray(ws_all[i],  float)
        strings = np.asarray(str_all[i], int)
        sensors = np.asarray(sen_all[i], int)

        zen = float(d["zenith_rad"][i]);  azi = float(d["azimuth_rad"][i])
        mc_dir = (math.sin(zen)*math.cos(azi),
                  math.sin(zen)*math.sin(azi),
                  math.cos(zen))

        for (s_dm, sen_dm), dm_pos in DMICE_DOMS.items():
            mask = (strings == s_dm) & (sensors == sen_dm)
            if not mask.any():
                continue
            ic_mask = ~mask
            if ic_mask.sum() < 4:
                continue
            t_geo, d_perp = t_geometric_dmice(
                xs[ic_mask], ys[ic_mask], zs[ic_mask],
                ts[ic_mask], ws[ic_mask], dm_pos, mc_dir)
            if t_geo is None or d_perp > d_perp_max:
                continue
            delta_t = float(ts[mask].min()) - t_geo
            if abs(delta_t) > 5000:
                continue
            records.append(dict(
                delta_t   = delta_t,
                d_perp    = d_perp,
                n_photons = float(ws[mask].sum()),
                energy    = float(d["energy_GeV"][i]),
                zenith    = math.degrees(zen),
            ))
    return records, N


def fit_gaussian_robust(dt_vals, clip_sigma=3.0):
    """Iterative Gaussian fit with outlier clipping."""
    mu, sig = np.median(dt_vals), np.std(dt_vals)
    for _ in range(5):
        keep = np.abs(dt_vals - mu) < clip_sigma * sig
        if keep.sum() < 5:
            break
        mu, sig = np.mean(dt_vals[keep]), np.std(dt_vals[keep])
    return mu, sig


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dperp-max", type=float, default=D_PERP_MAX_DEFAULT)
    parser.add_argument("--out",       default=OUT_DEFAULT)
    args = parser.parse_args()
    os.makedirs(args.out, exist_ok=True)

    # ── Load on-axis sims ─────────────────────────────────────────────────────
    all_records = []
    all_n_total = 0
    for path in NPZ_ON_AXIS:
        if not os.path.exists(path):
            continue
        recs, n_tot = extract_hits(path, args.dperp_max)
        all_records.extend(recs)
        all_n_total += n_tot
        print(f"  {os.path.basename(path):42s}  "
              f"{n_tot} events → {len(recs)} through hits (d⊥<{args.dperp_max:.0f}m)")

    print(f"\nTotal: {all_n_total} events, {len(all_records)} through-detector hits")
    if not all_records:
        print("No hits found."); return

    dt_arr  = np.array([r["delta_t"]   for r in all_records])
    dp_arr  = np.array([r["d_perp"]    for r in all_records])
    np_arr  = np.array([r["n_photons"] for r in all_records])
    en_arr  = np.array([r["energy"]    for r in all_records])
    zen_arr = np.array([r["zenith"]    for r in all_records])

    # ── Fit single Gaussian ───────────────────────────────────────────────────
    mu_fit, sig_fit = fit_gaussian_robust(dt_arr)
    eff_0   = len(all_records) / all_n_total
    eff_err = math.sqrt(eff_0 * (1 - eff_0) / all_n_total)

    print(f"\nGaussian fit (d⊥ < {args.dperp_max:.0f} m):")
    print(f"  μ  = {mu_fit:+.1f} ns")
    print(f"  σ  = {sig_fit:.1f} ns")
    print(f"  ε₀ = {eff_0:.3f} ± {eff_err:.3f}")
    print(f"  n  = {len(all_records)}")

    # ── Save model ────────────────────────────────────────────────────────────
    model_path = os.path.join(args.out, "dmice_timing_model.npz")
    np.savez(model_path,
        delta_t        = dt_arr,
        d_perp         = dp_arr,
        n_photons      = np_arr,
        energy         = en_arr,
        zenith         = zen_arr,
        mu_ns          = np.float64(mu_fit),
        sigma_ns       = np.float64(sig_fit),
        efficiency     = np.float64(eff_0),
        efficiency_err = np.float64(eff_err),
        d_perp_max_m   = np.float64(args.dperp_max),
        n_total_events = np.int64(all_n_total),
        n_dmice_hits   = np.int64(len(all_records)),
    )
    print(f"\nModel saved: {model_path}")

    # ── Plots ─────────────────────────────────────────────────────────────────
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    # 1. Δt distribution + Gaussian fit
    ax = axes[0]
    dt_lo, dt_hi = np.percentile(dt_arr, [0.5, 99.5])
    bins = np.linspace(dt_lo - 50, dt_hi + 50, 55)
    ax.hist(dt_arr, bins=bins, density=True, histtype="stepfilled",
            color="steelblue", alpha=0.5, label=f"Observed (n={len(dt_arr)})")
    ax.hist(dt_arr, bins=bins, density=True, histtype="step",
            color="steelblue", lw=1.5)
    x_fine = np.linspace(bins[0], bins[-1], 400)
    ax.plot(x_fine, sp_norm.pdf(x_fine, mu_fit, sig_fit),
            "r-", lw=2.5, label=f"Gaussian fit\nμ={mu_fit:+.0f} ns, σ={sig_fit:.0f} ns")
    ax.axvline(0,      color="k",   lw=0.8, ls="--", alpha=0.5, label="t_geo")
    ax.axvline(mu_fit, color="red", lw=1.2, ls=":",  alpha=0.8, label=f"μ={mu_fit:+.0f} ns")
    ax.set_xlabel("Δt = t_observed − t_geometric (ns)")
    ax.set_ylabel("Normalised events / bin")
    ax.set_title(f"DM-Ice timing residual\n(d⊥ < {args.dperp_max:.0f} m, through-detector only)")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # 2. Δt vs energy
    ax = axes[1]
    sc = ax.scatter(en_arr, dt_arr, c=zen_arr, s=14, alpha=0.6, cmap="plasma")
    fig.colorbar(sc, ax=ax, label="Zenith (°)")
    ax.set_xscale("log")
    ax.axhline(mu_fit, color="red", lw=1.5, ls="--", label=f"μ={mu_fit:+.0f} ns")
    ax.axhline(mu_fit + sig_fit, color="red", lw=0.8, ls=":", alpha=0.6)
    ax.axhline(mu_fit - sig_fit, color="red", lw=0.8, ls=":", alpha=0.6)
    ax.set_xlabel("Muon energy (GeV)")
    ax.set_ylabel("Δt (ns)")
    ax.set_title("Timing residual vs energy\n(coloured by zenith)")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3, which="both")

    # 3. Photon yield vs energy
    ax = axes[2]
    ax.scatter(en_arr, np_arr, c=dp_arr, s=14, alpha=0.6, cmap="viridis")
    sm = plt.cm.ScalarMappable(cmap="viridis",
                               norm=matplotlib.colors.Normalize(dp_arr.min(), dp_arr.max()))
    fig.colorbar(sm, ax=ax, label="d⊥ (m)")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("Muon energy (GeV)")
    ax.set_ylabel("DM-Ice photon count")
    ax.set_title("Photon yield vs energy\n(coloured by d⊥)")
    ax.grid(True, alpha=0.3, which="both")

    fig.suptitle(
        f"DM-Ice NaI timing model — through-detector events (d⊥ < {args.dperp_max:.0f} m)\n"
        f"{all_n_total} simulated, {len(all_records)} hits  |  "
        f"μ={mu_fit:+.0f} ns  σ={sig_fit:.0f} ns  ε={eff_0:.3f}",
        fontsize=11
    )
    plt.tight_layout()
    plot_path = os.path.join(args.out, "dmice_timing_residuals.png")
    fig.savefig(plot_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Plot saved:  {plot_path}")


if __name__ == "__main__":
    main()
