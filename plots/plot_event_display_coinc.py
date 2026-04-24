#!/usr/bin/env python3
"""
plot_event_display_coinc.py

Create event-display plots for 10 events that PASS and 10 that FAIL
the MPEFit-based Gaussian timing cut, to visually confirm which are
genuine muon-DM-Ice coincidences vs radioactivity accidentals.

For each event, plots 3 projections (x-z, y-z, x-y):
  - IC DOM hits: circle size ∝ charge, colour ∝ first-hit time (early=blue, late=red)
  - LineFit track extended through the full detector
  - DM-Ice crystal marked as a large gold star
  - Title includes delta_t_mpe, pass/fail, dm_t_ns, zenith

Cut used:
  |delta_t_mpe − 280| < 243 ns  →  PASS (genuine muon-DM-Ice)
  else                           →  FAIL (accidental radioactivity)

where delta_t_mpe = dm_t_ns − (mpe.time + dot(dm_pos − mpe.pos, d_hat) / C)

Run on cobalt:
  /cvmfs/icecube.opensciencegrid.org/py3-v4.3.0/RHEL_9_x86_64/metaprojects/icetray/v1.12.1/env-shell.sh \\
    python3 -u ~/dmice/plot_event_display_coinc.py [--year 2012] [--n-events 10]
"""

import os, sys, math, argparse, collections
import numpy as np

# ── Parameters ────────────────────────────────────────────────────────────────

MU_NS      = 280.0
SIGMA_NS   = 81.0
N_SIGMA    = 3.0
C_M_NS     = 0.2998   # m/ns

DMICE_POS_IC = {
    "det1": np.array([ 31.25,  -72.93, -511.05]),
    "det2": np.array([-334.80, -424.50, -511.26]),
}

GCD_FILE = ("/cvmfs/icecube.opensciencegrid.org/data/GCD/"
            "GeoCalibDetectorStatus_2013.56429_V1.i3.gz")
I3_FILE  = ("/data/user/bcharett/dmice_coincidences_2011_2022/"
            "all_dmice_coincidences_2011_2022_fixed.i3")
OUT_DIR  = os.path.expanduser("~/dmice_work/output/event_displays")

PULSE_PRIORITY = [
    "SplitInIcePulses", "OnlineL2_CleanedMuonPulses",
    "OfflinePulses", "SRTInIcePulses",
    "ReextractedInIcePulses", "InIcePulses",
]
MPE_KEYS     = ["MPEFit", "PoleMuonLlhFit"]
LF_KEYS      = ["LineFit", "PoleMuonLinefit"]
IC_STRINGS   = set(range(1, 87))
MUON_STREAMS = {'', 'in_ice', 'InIceSplit'}

# ── Args ──────────────────────────────────────────────────────────────────────

parser = argparse.ArgumentParser()
parser.add_argument("--year",     type=int, default=None,
                    help="Restrict to year (default: all, picks from 2012 first)")
parser.add_argument("--n-events", type=int, default=10,
                    help="How many pass + fail events to plot (default: 10 each)")
parser.add_argument("--gcd", default=GCD_FILE)
args = parser.parse_args()
N_EVENTS = args.n_events

os.makedirs(OUT_DIR, exist_ok=True)

# ── IceTray ───────────────────────────────────────────────────────────────────

from icecube import icetray, dataio, dataclasses

# ── Helpers ───────────────────────────────────────────────────────────────────

def get_track(frame, keys):
    for k in keys:
        if k in frame:
            p = frame[k]
            if hasattr(p, "fit_status") and p.fit_status == dataclasses.I3Particle.FitStatus.OK:
                return p, k
    return None, None

def compute_delta_t_mpe(mpe, dm_pos, dm_t_ns):
    """Δt_mpe = dm_t_ns − t_PCA(mpe → dm_pos)."""
    d_hat = np.array([mpe.dir.x, mpe.dir.y, mpe.dir.z])
    r     = dm_pos - np.array([mpe.pos.x, mpe.pos.y, mpe.pos.z])
    s     = float(np.dot(r, d_hat))
    t_pca = mpe.time + s / C_M_NS
    return dm_t_ns - t_pca

def get_pulses(frame):
    """Return best available pulse series map (IC strings only)."""
    for key in PULSE_PRIORITY:
        if key not in frame:
            continue
        raw = frame[key]
        # May be an I3RecoPulseSeriesMapMask — apply it
        try:
            pmap = dataclasses.I3RecoPulseSeriesMap.from_frame(frame, key)
        except Exception:
            try:
                pmap = raw.apply(frame)
            except Exception:
                pmap = raw
        return pmap, key
    return None, None

def get_om_positions(frame):
    """Return dict (string, om) → (x, y, z) for IC strings."""
    if "I3Geometry" not in frame:
        return {}
    geo = frame["I3Geometry"].omgeo
    pos = {}
    for omk, omg in geo.items():
        if omk.string in IC_STRINGS:
            pos[(omk.string, omk.om)] = np.array([omg.position.x,
                                                   omg.position.y,
                                                   omg.position.z])
    return pos

# ── Scan the i3 file ──────────────────────────────────────────────────────────

print(f"Opening: {I3_FILE}")
print(f"GCD:     {args.gcd}")

pass_events = []   # each entry: dict with all we need to plot
fail_events = []
om_pos      = {}   # loaded once from first geometry frame
seen        = set()
n_scanned   = 0
n_target    = N_EVENTS  # how many we want in each group

f = dataio.I3File([args.gcd, I3_FILE])
while f.more() and (len(pass_events) < n_target or len(fail_events) < n_target):
    frame = f.pop_frame()

    # Load geometry from G frame
    if frame.Stop == icetray.I3Frame.Geometry and not om_pos:
        om_pos = get_om_positions(frame)
        print(f"  Geometry loaded: {len(om_pos)} IC OMs")
        continue

    if frame.Stop != icetray.I3Frame.Physics:
        continue

    hdr    = frame["I3EventHeader"]
    stream = getattr(hdr, "sub_event_stream", "")
    if stream not in MUON_STREAMS:
        continue

    year = hdr.start_time.utc_year
    if args.year and year != args.year:
        continue

    uid = (hdr.run_id, hdr.event_id, stream)
    if uid in seen:
        continue
    seen.add(uid)
    n_scanned += 1

    # Need MPEFit for Gaussian cut
    mpe, mpe_key = get_track(frame, MPE_KEYS)
    if mpe is None:
        continue  # skip events without MPEFit

    lf, lf_key = get_track(frame, LF_KEYS)
    if lf is None:
        continue  # need LineFit for the track display

    # DM-Ice hit time
    if "DMIce_detection_time" not in frame:
        continue

    det_str = str(frame["DMIce_detector"]) if "DMIce_detector" in frame else "det1"
    det_key = "det1" if "det1" in det_str else "det2"
    dm_pos  = DMICE_POS_IC[det_key]

    event_start_daq = hdr.start_time.utc_daq_time
    dm_t_ns = (frame["DMIce_detection_time"].value - event_start_daq) * 0.1

    delta_t = compute_delta_t_mpe(mpe, dm_pos, dm_t_ns)
    passes  = abs(delta_t - MU_NS) < N_SIGMA * SIGMA_NS

    # Only collect what we need
    if passes and len(pass_events) >= n_target:
        continue
    if not passes and len(fail_events) >= n_target:
        continue

    # Extract pulse hit pattern
    pmap, pulse_key = get_pulses(frame)
    hits = []  # list of (x, y, z, charge, first_time)
    if pmap is not None:
        for omk, pulses in pmap.items():
            if omk.string not in IC_STRINGS:
                continue
            key_t = (omk.string, omk.om)
            if key_t not in om_pos:
                continue
            x, y, z = om_pos[key_t]
            charge    = sum(p.charge for p in pulses)
            first_t   = min(p.time   for p in pulses)
            hits.append((x, y, z, charge, first_t))

    if not hits:
        continue  # skip events with no IC hits

    ev = dict(
        run_id    = hdr.run_id,
        event_id  = hdr.event_id,
        year      = year,
        passes    = passes,
        delta_t   = delta_t,
        dm_t_ns   = dm_t_ns,
        det_key   = det_key,
        dm_pos    = dm_pos.copy(),
        # LineFit
        lf_pos    = np.array([lf.pos.x,  lf.pos.y,  lf.pos.z ]),
        lf_dir    = np.array([lf.dir.x,  lf.dir.y,  lf.dir.z ]),
        lf_zen    = math.degrees(lf.dir.zenith),
        # MPEFit (for reference)
        mpe_zen   = math.degrees(mpe.dir.zenith),
        # Hits
        hits      = hits,
        pulse_key = pulse_key,
    )

    if passes:
        pass_events.append(ev)
        print(f"  PASS [{len(pass_events):2d}/{n_target}] "
              f"run={hdr.run_id} evt={hdr.event_id} "
              f"Δt={delta_t:.0f} ns  dm_t={dm_t_ns/1000:.1f} μs  "
              f"zen_lf={ev['lf_zen']:.1f}°  nhits={len(hits)}")
    else:
        fail_events.append(ev)
        print(f"  FAIL [{len(fail_events):2d}/{n_target}] "
              f"run={hdr.run_id} evt={hdr.event_id} "
              f"Δt={delta_t:.0f} ns  dm_t={dm_t_ns/1000:.1f} μs  "
              f"zen_lf={ev['lf_zen']:.1f}°  nhits={len(hits)}")

f.close()
print(f"\nScanned {n_scanned} events  "
      f"→  PASS: {len(pass_events)}  FAIL: {len(fail_events)}")

# ── Plotting ──────────────────────────────────────────────────────────────────

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

def track_endpoints(pos, d_hat, z_min=-600, z_max=600):
    """Extend LineFit track to z limits; return two endpoints."""
    # t parameter so that pos + t*d_hat has z = z_min or z_max
    if abs(d_hat[2]) < 1e-6:
        t_lo, t_hi = -1000, 1000
    else:
        t_lo = (z_min - pos[2]) / d_hat[2]
        t_hi = (z_max - pos[2]) / d_hat[2]
        if t_lo > t_hi:
            t_lo, t_hi = t_hi, t_lo
    p0 = pos + t_lo * d_hat
    p1 = pos + t_hi * d_hat
    return p0, p1

def plot_event(ev, label, idx, group):
    """Plot one event: 3 projection panels (x-z, y-z, x-y)."""
    hits = ev["hits"]
    if not hits:
        return

    xs   = np.array([h[0] for h in hits])
    ys   = np.array([h[1] for h in hits])
    zs   = np.array([h[2] for h in hits])
    qs   = np.array([h[3] for h in hits])
    ts   = np.array([h[4] for h in hits])

    # Normalise charge for marker size; normalise time for colour
    q_size = 20 + 200 * (qs / qs.max()) ** 0.5
    t_norm = Normalize(vmin=ts.min(), vmax=ts.max())
    cmap   = plt.cm.plasma_r  # early=bright/yellow, late=dark/purple

    lf_pos = ev["lf_pos"]
    lf_dir = ev["lf_dir"]
    dm_pos = ev["dm_pos"]
    p0, p1 = track_endpoints(lf_pos, lf_dir)

    status  = "PASS" if ev["passes"] else "FAIL"
    color   = "green" if ev["passes"] else "red"
    title   = (f"{status}  run={ev['run_id']} evt={ev['event_id']} ({ev['year']})\n"
               f"Δt_mpe={ev['delta_t']:.0f} ns  dm_t={ev['dm_t_ns']/1000:.1f} μs  "
               f"LF zen={ev['lf_zen']:.1f}°  MPE zen={ev['mpe_zen']:.1f}°  "
               f"N_hits={len(hits)}")

    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    fig.patch.set_facecolor("#111111")
    for ax in axes:
        ax.set_facecolor("#111111")

    # ── Panel 0: x–z view ────────────────────────────────────────────────────
    ax = axes[0]
    sc = ax.scatter(xs, zs, c=ts, cmap=cmap, norm=t_norm,
                    s=q_size, alpha=0.85, zorder=3, linewidths=0)
    ax.plot([p0[0], p1[0]], [p0[2], p1[2]],
            color="cyan", lw=1.5, alpha=0.8, label="LineFit", zorder=4)
    ax.scatter([dm_pos[0]], [dm_pos[2]], marker="*", s=600,
               color="gold", zorder=5, label=f"DM-Ice {ev['det_key']}")
    ax.set_xlabel("x [m]", color="white"); ax.set_ylabel("z [m]", color="white")
    ax.set_title("x–z", color="white")
    ax.tick_params(colors="white")
    for spine in ax.spines.values(): spine.set_edgecolor("#555555")
    ax.legend(fontsize=8, labelcolor="white",
              facecolor="#222222", edgecolor="#555555")

    # ── Panel 1: y–z view ────────────────────────────────────────────────────
    ax = axes[1]
    ax.scatter(ys, zs, c=ts, cmap=cmap, norm=t_norm,
               s=q_size, alpha=0.85, zorder=3, linewidths=0)
    ax.plot([p0[1], p1[1]], [p0[2], p1[2]],
            color="cyan", lw=1.5, alpha=0.8, label="LineFit", zorder=4)
    ax.scatter([dm_pos[1]], [dm_pos[2]], marker="*", s=600,
               color="gold", zorder=5, label=f"DM-Ice {ev['det_key']}")
    ax.set_xlabel("y [m]", color="white"); ax.set_ylabel("z [m]", color="white")
    ax.set_title("y–z", color="white")
    ax.tick_params(colors="white")
    for spine in ax.spines.values(): spine.set_edgecolor("#555555")

    # ── Panel 2: x–y view (top-down) ─────────────────────────────────────────
    ax = axes[2]
    ax.scatter(xs, ys, c=ts, cmap=cmap, norm=t_norm,
               s=q_size, alpha=0.85, zorder=3, linewidths=0)
    ax.plot([p0[0], p1[0]], [p0[1], p1[1]],
            color="cyan", lw=1.5, alpha=0.8, label="LineFit", zorder=4)
    ax.scatter([dm_pos[0]], [dm_pos[1]], marker="*", s=600,
               color="gold", zorder=5, label=f"DM-Ice {ev['det_key']}")
    ax.set_xlabel("x [m]", color="white"); ax.set_ylabel("y [m]", color="white")
    ax.set_title("x–y (top-down)", color="white")
    ax.tick_params(colors="white")
    for spine in ax.spines.values(): spine.set_edgecolor("#555555")

    # Colourbar
    sm = ScalarMappable(cmap=cmap, norm=t_norm)
    sm.set_array([])
    cb = fig.colorbar(sm, ax=axes[2], pad=0.01, fraction=0.046)
    cb.set_label("First hit time [ns]", color="white")
    cb.ax.yaxis.set_tick_params(color="white")
    plt.setp(cb.ax.yaxis.get_ticklabels(), color="white")

    # Title
    fig.suptitle(title, fontsize=10, color=color, y=1.01)
    plt.tight_layout()

    fname = os.path.join(OUT_DIR, f"{group}_{idx:02d}_{status.lower()}_"
                                   f"r{ev['run_id']}_e{ev['event_id']}.png")
    fig.savefig(fname, dpi=120, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)
    print(f"  Saved: {os.path.basename(fname)}")

print(f"\nPlotting {len(pass_events)} PASS events …")
for i, ev in enumerate(pass_events):
    plot_event(ev, label="PASS", idx=i+1, group="pass")

print(f"\nPlotting {len(fail_events)} FAIL events …")
for i, ev in enumerate(fail_events):
    plot_event(ev, label="FAIL", idx=i+1, group="fail")

# ── Summary montage ───────────────────────────────────────────────────────────

print("\nBuilding comparison montage …")
all_ev  = pass_events[:5] + fail_events[:5]
n_shown = len(all_ev)

fig, axes = plt.subplots(n_shown, 3, figsize=(18, 3.5 * n_shown))
fig.patch.set_facecolor("#111111")

for row, ev in enumerate(all_ev):
    hits  = ev["hits"]
    if not hits:
        continue
    xs   = np.array([h[0] for h in hits])
    ys   = np.array([h[1] for h in hits])
    zs   = np.array([h[2] for h in hits])
    qs   = np.array([h[3] for h in hits])
    ts   = np.array([h[4] for h in hits])
    q_size = 10 + 80 * (qs / qs.max()) ** 0.5
    t_norm = Normalize(vmin=ts.min(), vmax=ts.max())
    cmap   = plt.cm.plasma_r

    lf_pos = ev["lf_pos"]
    lf_dir = ev["lf_dir"]
    dm_pos = ev["dm_pos"]
    p0, p1 = track_endpoints(lf_pos, lf_dir)

    status = "PASS" if ev["passes"] else "FAIL"
    color  = "lime"  if ev["passes"] else "red"

    for col, (xd, yd, xl, yl, xdm, ydm, xlabel, ylabel, title) in enumerate([
        (xs, zs, [p0[0], p1[0]], [p0[2], p1[2]], dm_pos[0], dm_pos[2], "x", "z", "x–z"),
        (ys, zs, [p0[1], p1[1]], [p0[2], p1[2]], dm_pos[1], dm_pos[2], "y", "z", "y–z"),
        (xs, ys, [p0[0], p1[0]], [p0[1], p1[1]], dm_pos[0], dm_pos[1], "x", "y", "x–y"),
    ]):
        ax = axes[row, col]
        ax.set_facecolor("#111111")
        ax.scatter(xd, yd, c=ts, cmap=cmap, norm=t_norm,
                   s=q_size, alpha=0.85, zorder=3, linewidths=0)
        ax.plot(xl, yl, color="cyan", lw=1.2, alpha=0.8, zorder=4)
        ax.scatter([xdm], [ydm], marker="*", s=350, color="gold", zorder=5)
        ax.tick_params(colors="white", labelsize=7)
        for spine in ax.spines.values():
            spine.set_edgecolor("#444444")
        if col == 0:
            ax.set_ylabel(f"{status}\nΔt={ev['delta_t']:.0f}ns\n"
                          f"dm_t={ev['dm_t_ns']/1000:.1f}μs\n"
                          f"zen={ev['lf_zen']:.0f}°",
                          color=color, fontsize=8, rotation=0, ha="right",
                          va="center", labelpad=60)
        if row == 0:
            ax.set_title(title, color="white", fontsize=9)

fig.suptitle(
    f"Event display comparison — top 5 PASS (green) vs top 5 FAIL (red)\n"
    f"Gaussian cut: |Δt_mpe − {MU_NS:.0f}| < {N_SIGMA*SIGMA_NS:.0f} ns  "
    f"◆ cyan = LineFit track  ★ = DM-Ice crystal  colour = hit time",
    fontsize=10, color="white", y=1.005
)
plt.tight_layout()
montage_path = os.path.join(OUT_DIR, "montage_pass_vs_fail.png")
fig.savefig(montage_path, dpi=120, bbox_inches="tight", facecolor=fig.get_facecolor())
plt.close(fig)
print(f"Montage: {montage_path}")
print("\nDone.")
