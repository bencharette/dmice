#!/usr/bin/env python3
"""
Pure-numpy LineFit vs DM-Ice Pivot LineFit on Prometheus simulation output.

No IceTray required. Works directly on the Prometheus parquet output.

For each event:
  1. Load photon hits (sensor positions + times) from parquet
  2. Reduce to first-photon time per DOM (charge weight = 1 per DOM)
  3. Run IC-only LineFit (exclude DM-Ice string 87/88)
  4. Run DM-Ice Pivot LineFit using earliest DM-Ice hit as the pivot
     (or MC truth transit time if no DM-Ice hit)
  5. Compare angular errors vs MC truth muon direction

Usage:
    python npz_linefit_comparison.py [--parquet FILE] [--min-hits N] [--plot OUT.png]
"""

import argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd

# DM-Ice string IDs in Prometheus simulation
DMICE_STRINGS = {87, 88}

# Speed of light
C_M_NS = 0.2998  # m/ns


# ── Direction utilities ──────────────────────────────────────────────────────

def mc_travel_dir(zenith_rad, azimuth_rad):
    """
    Convert Prometheus MC truth zenith/azimuth to travel direction vector.

    Prometheus parquet stores initial_state_zenith/azimuth in the Prometheus
    momentum convention (zenith=0 = traveling up, zenith=180 = traveling down).
    The momentum direction IS the travel direction — no negation needed.
    """
    sz = np.sin(zenith_rad)
    cz = np.cos(zenith_rad)
    sa = np.sin(azimuth_rad)
    ca = np.cos(azimuth_rad)
    return np.array([sz * ca, sz * sa, cz])


def angular_diff_deg(d1, d2):
    """Great-circle angle between two unit vectors (degrees)."""
    dot = float(np.clip(np.dot(d1, d2), -1.0, 1.0))
    return np.degrees(np.arccos(dot))


# ── LineFit implementations ──────────────────────────────────────────────────

def run_linefit(x, y, z, t, w=None):
    """
    Standard analytic LineFit (charge-weighted centre-of-gravity pivot).
    Returns velocity direction (travel direction) dict, or None.
    """
    if w is None:
        w = np.ones(len(t))
    W = np.sum(w)
    if W == 0 or len(t) < 3:
        return None
    x_bar = np.dot(w, x) / W
    y_bar = np.dot(w, y) / W
    z_bar = np.dot(w, z) / W
    t_bar = np.dot(w, t) / W
    dt = t - t_bar
    denom = np.dot(w, dt ** 2)
    if denom == 0:
        return None
    vx = np.dot(w * dt, x - x_bar) / denom
    vy = np.dot(w * dt, y - y_bar) / denom
    vz = np.dot(w * dt, z - z_bar) / denom
    speed = np.sqrt(vx ** 2 + vy ** 2 + vz ** 2)
    if speed == 0:
        return None
    dx, dy, dz = vx / speed, vy / speed, vz / speed
    return dict(dx=dx, dy=dy, dz=dz, speed_m_ns=speed)


def run_pivot_linefit(x_ic, y_ic, z_ic, t_ic, w_ic,
                      x_dm, y_dm, z_dm, t_dm):
    """
    Pivot LineFit: DM-Ice detector is the fixed space+time reference.

        dt_i = t_i - t_dm
        dr_i = r_i - r_dm
        v = Σ(w_i * dt_i * dr_i) / Σ(w_i * dt_i²)

    Returns velocity direction dict, or None.
    """
    dt = t_ic - t_dm
    dr_x = x_ic - x_dm
    dr_y = y_ic - y_dm
    dr_z = z_ic - z_dm
    denom = np.dot(w_ic, dt ** 2)
    if denom == 0:
        return None
    vx = np.dot(w_ic * dt, dr_x) / denom
    vy = np.dot(w_ic * dt, dr_y) / denom
    vz = np.dot(w_ic * dt, dr_z) / denom
    speed = np.sqrt(vx ** 2 + vy ** 2 + vz ** 2)
    if speed == 0:
        return None
    dx, dy, dz = vx / speed, vy / speed, vz / speed
    return dict(dx=dx, dy=dy, dz=dz, speed_m_ns=speed)


def run_pivot_linefit_iterative(x_ic, y_ic, z_ic, t_ic, w_ic,
                                x_dm, y_dm, z_dm, t_dm,
                                huber_sigma=200.0, max_iter=10):
    """
    Iterative Pivot LineFit with Huber-style outlier downweighting.
    Converges when direction shift < 0.01 degrees.
    """
    w = w_ic.copy()
    prev_dir = None

    for _ in range(max_iter):
        dt = t_ic - t_dm
        dr_x = x_ic - x_dm
        dr_y = y_ic - y_dm
        dr_z = z_ic - z_dm
        denom = np.dot(w, dt ** 2)
        if denom == 0:
            return None
        vx = np.dot(w * dt, dr_x) / denom
        vy = np.dot(w * dt, dr_y) / denom
        vz = np.dot(w * dt, dr_z) / denom
        speed = np.sqrt(vx ** 2 + vy ** 2 + vz ** 2)
        if speed == 0:
            return None

        t_expected = (dr_x * vx + dr_y * vy + dr_z * vz) / speed ** 2
        residual = dt - t_expected
        w = w_ic * np.minimum(1.0, huber_sigma / (np.abs(residual) + 1e-6))

        v_hat = np.array([vx, vy, vz]) / speed
        if prev_dir is not None:
            shift = np.degrees(np.arccos(np.clip(np.dot(v_hat, prev_dir), -1.0, 1.0)))
            if shift < 0.01:
                break
        prev_dir = v_hat

    dx, dy, dz = v_hat
    return dict(dx=dx, dy=dy, dz=dz, speed_m_ns=speed)


# ── DM-Ice MC truth transit time ─────────────────────────────────────────────

def compute_mc_dmice_time(x_ic, y_ic, z_ic, t_ic, w_ic, dm_pos, mc_dir_travel):
    """
    Compute MC truth DM-Ice hit time in the same time frame as IC hits.

    Uses charge-weighted IC centroid as reference; projects (dm_pos - centroid)
    onto the MC travel direction to get the time offset.
    """
    W = np.sum(w_ic)
    if W == 0:
        return None
    r_bar = np.array([np.dot(w_ic, x_ic), np.dot(w_ic, y_ic), np.dot(w_ic, z_ic)]) / W
    t_bar = np.dot(w_ic, t_ic) / W
    d = np.dot(dm_pos - r_bar, mc_dir_travel)   # signed distance along track (m)
    return t_bar + d / C_M_NS


# ── Main ─────────────────────────────────────────────────────────────────────

def get_muon_state(mc):
    """Return (zenith_rad, azimuth_rad, energy, x, y, z) for the highest-energy muon."""
    types = np.array(mc['final_state_type'])
    mu_mask = np.abs(types) == 13
    if not mu_mask.any():
        # Fall back to initial state
        return (mc['initial_state_zenith'], mc['initial_state_azimuth'],
                mc['initial_state_energy'],
                mc['initial_state_x'], mc['initial_state_y'], mc['initial_state_z'])
    energies = np.array(mc['final_state_energy'])
    idx = np.where(mu_mask)[0][np.argmax(energies[mu_mask])]
    zens = np.array(mc['final_state_zenith'])
    azs  = np.array(mc['final_state_azimuth'])
    xs   = np.array(mc['final_state_x'])
    ys   = np.array(mc['final_state_y'])
    zs   = np.array(mc['final_state_z'])
    return (zens[idx], azs[idx], energies[idx], xs[idx], ys[idx], zs[idx])


def main():
    parser = argparse.ArgumentParser(
        description='LineFit vs DM-Ice Pivot LineFit on Prometheus NPZ/parquet data')
    parser.add_argument('--parquet', default='output/sim_run8_photons.parquet',
                        help='Prometheus parquet file')
    parser.add_argument('--min-hits', type=int, default=100,
                        help='Minimum photon hits per event (default: 100)')
    parser.add_argument('--min-doms', type=int, default=4,
                        help='Minimum unique DOMs after first-hit reduction (default: 4)')
    parser.add_argument('--output', default='output/npz_linefit_results.csv')
    parser.add_argument('--plot', default=None, metavar='FILE.png')
    args = parser.parse_args()

    df = pd.read_parquet(args.parquet)
    n_total = len(df)
    print(f"Loaded {n_total} events from {args.parquet}")

    results = []
    n_pass_hits = 0
    n_pass_doms = 0
    n_with_dm   = 0

    for ev_idx, row in df.iterrows():
        ph = row['photons']
        mc = row['mc_truth']

        t_all   = np.array(ph['t'])
        sx_all  = np.array(ph['sensor_pos_x'])
        sy_all  = np.array(ph['sensor_pos_y'])
        sz_all  = np.array(ph['sensor_pos_z'])
        str_all = np.array(ph['string_id'])
        sns_all = np.array(ph['sensor_id'])

        # ── Filter: minimum photon hits ──────────────────────────────────────
        if len(t_all) < args.min_hits:
            continue
        n_pass_hits += 1

        # ── Reduce to first-hit per DOM ──────────────────────────────────────
        dom_first = {}
        for i in range(len(t_all)):
            key = (int(str_all[i]), int(sns_all[i]))
            if key not in dom_first or t_all[i] < dom_first[key][0]:
                dom_first[key] = (t_all[i], sx_all[i], sy_all[i], sz_all[i])

        # Separate IC and DM-Ice DOMs
        ic_hits  = [(t, x, y, z) for (s, n), (t, x, y, z) in dom_first.items()
                    if s not in DMICE_STRINGS]
        dm_hits  = [(t, x, y, z) for (s, n), (t, x, y, z) in dom_first.items()
                    if s in DMICE_STRINGS]

        if len(ic_hits) < args.min_doms:
            continue
        n_pass_doms += 1

        ic_t = np.array([h[0] for h in ic_hits])
        ic_x = np.array([h[1] for h in ic_hits])
        ic_y = np.array([h[2] for h in ic_hits])
        ic_z = np.array([h[3] for h in ic_hits])
        ic_w = np.ones(len(ic_hits))  # unit charge (one hit per DOM)

        # ── MC truth muon direction ──────────────────────────────────────────
        zen, az, energy, mu_x, mu_y, mu_z = get_muon_state(mc)
        mc_dir = mc_travel_dir(zen, az)  # travel direction

        # ── IC-only LineFit ──────────────────────────────────────────────────
        lf_ic = run_linefit(ic_x, ic_y, ic_z, ic_t, ic_w)
        ic_ang_err = np.nan
        ic_speed   = np.nan
        if lf_ic is not None:
            ic_vel = np.array([lf_ic['dx'], lf_ic['dy'], lf_ic['dz']])
            ic_ang_err = angular_diff_deg(mc_dir, ic_vel)
            ic_speed   = lf_ic['speed_m_ns']

        # ── DM-Ice Pivot LineFit ─────────────────────────────────────────────
        pivot_ang_err      = np.nan
        pivot_iter_ang_err = np.nan
        pivot_speed        = np.nan
        dm_hit_source      = 'none'
        dm_pos             = None
        t_dm               = None

        if dm_hits:
            # Use earliest real simulated DM-Ice hit
            dm_hits.sort(key=lambda h: h[0])
            t_dm, x_dm, y_dm, z_dm = dm_hits[0]
            dm_pos = np.array([x_dm, y_dm, z_dm])
            dm_hit_source = 'simulated'
            n_with_dm += 1
        elif lf_ic is not None:
            # Fall back to MC truth transit time
            # Find closest DM-Ice position from the geo file
            # String 87 (det1): (31.25, -72.93, -2459.12) in Prometheus coords
            dm_candidates = {
                'det1': np.array([31.25,    -72.93,   -2459.12]),
                'det2': np.array([-334.80, -424.50,   -2459.26]),
            }
            # Closest approach to muon track
            mu_pos = np.array([mu_x, mu_y, mu_z])
            best_dist = np.inf
            for det, pos in dm_candidates.items():
                dp = pos - mu_pos
                proj = np.dot(dp, mc_dir)
                closest = mu_pos + proj * mc_dir
                dist = np.linalg.norm(pos - closest)
                if dist < best_dist:
                    best_dist = dist
                    dm_pos = pos

            t_dm = compute_mc_dmice_time(ic_x, ic_y, ic_z, ic_t, ic_w, dm_pos, mc_dir)
            dm_hit_source = 'mc_truth'

        if dm_pos is not None and t_dm is not None and lf_ic is not None:
            x_dm, y_dm, z_dm = dm_pos

            lf_piv = run_pivot_linefit(
                ic_x, ic_y, ic_z, ic_t, ic_w,
                x_dm, y_dm, z_dm, t_dm)
            if lf_piv is not None:
                piv_vel = np.array([lf_piv['dx'], lf_piv['dy'], lf_piv['dz']])
                pivot_ang_err = angular_diff_deg(mc_dir, piv_vel)
                pivot_speed   = lf_piv['speed_m_ns']

            lf_piv_iter = run_pivot_linefit_iterative(
                ic_x, ic_y, ic_z, ic_t, ic_w,
                x_dm, y_dm, z_dm, t_dm)
            if lf_piv_iter is not None:
                piv_iter_vel = np.array([lf_piv_iter['dx'], lf_piv_iter['dy'], lf_piv_iter['dz']])
                pivot_iter_ang_err = angular_diff_deg(mc_dir, piv_iter_vel)

        results.append(dict(
            ev_idx           = ev_idx,
            n_photons        = len(t_all),
            n_ic_doms        = len(ic_hits),
            n_dm_hits        = len(dm_hits),
            dm_hit_source    = dm_hit_source,
            mc_zenith_deg    = np.degrees(zen),
            mc_azimuth_deg   = np.degrees(az),
            mc_energy_GeV    = energy,
            ic_ang_err_deg   = ic_ang_err,
            ic_speed_m_ns    = ic_speed,
            pivot_ang_err_deg      = pivot_ang_err,
            pivot_iter_ang_err_deg = pivot_iter_ang_err,
            pivot_speed_m_ns = pivot_speed,
        ))

    print(f"\nTotal events:             {n_total}")
    print(f"Pass >= {args.min_hits} photon hits:  {n_pass_hits}")
    print(f"Pass >= {args.min_doms} IC DOMs:       {n_pass_doms}")
    print(f"With real DM-Ice hit:     {n_with_dm}")

    if not results:
        print("No events passed cuts.")
        return

    out = pd.DataFrame(results)
    out.to_csv(args.output, index=False)
    print(f"\nSaved {len(out)} rows → {args.output}")

    # ── Summary ──────────────────────────────────────────────────────────────
    has_ic  = out['ic_ang_err_deg'].notna()
    has_piv = out['pivot_ang_err_deg'].notna()
    has_both = has_ic & has_piv

    print("\n── IC-only LineFit angular error vs MC truth ────────────────────────")
    print(out.loc[has_ic, 'ic_ang_err_deg'].describe().to_string())

    if has_piv.any():
        print("\n── Pivot LineFit angular error vs MC truth ──────────────────────────")
        print(out.loc[has_piv, 'pivot_ang_err_deg'].describe().to_string())

    if has_both.any():
        improved = out.loc[has_both, 'pivot_ang_err_deg'] < out.loc[has_both, 'ic_ang_err_deg']
        print(f"\n── DM-Ice improves direction: {improved.sum()}/{has_both.sum()} events ──")
        print("  IC-only LineFit median:          {:.2f} deg".format(
            out.loc[has_both, 'ic_ang_err_deg'].median()))
        print("  DM-Ice Pivot LineFit median:     {:.2f} deg".format(
            out.loc[has_both, 'pivot_ang_err_deg'].median()))
        has_iter = out['pivot_iter_ang_err_deg'].notna()
        if has_iter.any():
            print("  DM-Ice Pivot Iterative median:   {:.2f} deg".format(
                out.loc[has_iter, 'pivot_iter_ang_err_deg'].median()))

    print("\n── Per-event breakdown ──────────────────────────────────────────────")
    cols = ['ev_idx', 'n_photons', 'n_ic_doms', 'n_dm_hits', 'dm_hit_source',
            'mc_zenith_deg', 'mc_energy_GeV',
            'ic_ang_err_deg', 'pivot_ang_err_deg', 'pivot_iter_ang_err_deg']
    pd.set_option('display.float_format', '{:.2f}'.format)
    pd.set_option('display.max_columns', 20)
    pd.set_option('display.width', 140)
    print(out[cols].to_string(index=False))

    # ── Plot ─────────────────────────────────────────────────────────────────
    if args.plot:
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        fig.suptitle(
            f'LineFit vs DM-Ice Pivot LineFit — Prometheus simulation\n'
            f'(n={len(out)} events, ≥{args.min_hits} photon hits, ≥{args.min_doms} IC DOMs)',
            fontsize=12)

        # Left: angular error histograms
        ax = axes[0]
        ic_err   = out.loc[has_ic,  'ic_ang_err_deg']
        piv_err  = out.loc[has_piv, 'pivot_ang_err_deg']
        piv_iter = out.loc[out['pivot_iter_ang_err_deg'].notna(), 'pivot_iter_ang_err_deg']

        all_vals = pd.concat([ic_err, piv_err, piv_iter]).dropna()
        if len(all_vals) > 0:
            bins = np.linspace(0, min(all_vals.max() * 1.05, 90), 30)
            if len(ic_err) > 0:
                ax.hist(ic_err, bins=bins, histtype='stepfilled', alpha=0.55,
                        color='steelblue', edgecolor='steelblue',
                        label=f'IC-only LineFit  median={ic_err.median():.1f}°')
                ax.axvline(ic_err.median(), color='navy', lw=2, ls='--')
            if len(piv_err) > 0:
                ax.hist(piv_err, bins=bins, histtype='stepfilled', alpha=0.45,
                        color='tomato', edgecolor='tomato',
                        label=f'DM-Ice Pivot LineFit  median={piv_err.median():.1f}°')
                ax.axvline(piv_err.median(), color='darkred', lw=2, ls='--')
            if len(piv_iter) > 0:
                ax.hist(piv_iter, bins=bins, histtype='step', lw=2,
                        color='darkorange',
                        label=f'Pivot Iterative  median={piv_iter.median():.1f}°')
                ax.axvline(piv_iter.median(), color='darkorange', lw=2, ls='--')

        ax.set_xlabel('Angular error vs MC truth (deg)', fontsize=11)
        ax.set_ylabel('Events', fontsize=11)
        ax.legend(fontsize=9)
        ax.set_title('Angular resolution', fontsize=11)

        # Right: scatter IC err vs Pivot err per event
        ax2 = axes[1]
        if has_both.any():
            x_vals = out.loc[has_both, 'ic_ang_err_deg']
            y_vals = out.loc[has_both, 'pivot_ang_err_deg']
            # Color by DM-Ice hit source
            colors = ['green' if src == 'simulated' else 'steelblue'
                      for src in out.loc[has_both, 'dm_hit_source']]
            ax2.scatter(x_vals, y_vals, c=colors, s=60, zorder=5)
            lim = max(x_vals.max(), y_vals.max()) * 1.05
            ax2.plot([0, lim], [0, lim], 'k--', lw=1, label='No improvement')
            ax2.set_xlabel('IC-only angular error (deg)', fontsize=11)
            ax2.set_ylabel('Pivot angular error (deg)', fontsize=11)
            ax2.set_xlim(0, lim); ax2.set_ylim(0, lim)
            # legend for colors
            from matplotlib.patches import Patch
            legend_els = [Patch(fc='green', label='Real DM-Ice hit'),
                          Patch(fc='steelblue', label='MC truth time')]
            ax2.legend(handles=legend_els, fontsize=9)
            ax2.set_title('IC-only vs Pivot (below diagonal = DM-Ice improves)', fontsize=11)

        plt.tight_layout()
        plt.savefig(args.plot, dpi=150)
        print(f"Saved plot → {args.plot}")


if __name__ == '__main__':
    main()
