"""
plot_dmice_timing_diagram.py

Two-slide diagram explaining DM-Ice NaI timing for muon reconstruction.

Slide 1 — Geometry:
  Schematic of IceCube + DM-Ice showing a muon track, DOM hit sequence,
  and the DM-Ice crystal at the bottom. Illustrates d_perp ≈ 0 (direct transit).

Slide 2 — Timing chain:
  Timeline showing t_geo → scintillation delay → t_observed, with the
  Gaussian(μ=+280ns, σ=81ns) residual distribution from real data.

Run locally (no IceTray):
    python3 plot_dmice_timing_diagram.py

Saves:
    ~/dmice_work/output/dmice_timing_slide1_geometry.png
    ~/dmice_work/output/dmice_timing_slide2_timing.png
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
from matplotlib.lines import Line2D

OUT = os.path.expanduser("~/dmice_work/output")
os.makedirs(OUT, exist_ok=True)

MU_NS    = 280.0
SIGMA_NS =  81.0

# ══════════════════════════════════════════════════════════════════════════════
# SLIDE 1 — Geometry diagram
# ══════════════════════════════════════════════════════════════════════════════

fig, ax = plt.subplots(figsize=(10, 8))
ax.set_xlim(-1, 11)
ax.set_ylim(-1, 12)
ax.set_aspect("equal")
ax.axis("off")
fig.patch.set_facecolor("#0d1117")
ax.set_facecolor("#0d1117")

WHITE  = "#e8eaed"
BLUE   = "#4fc3f7"
ORANGE = "#ffb74d"
GREEN  = "#81c784"
RED    = "#ef5350"
CYAN   = "#80deea"
GRAY   = "#78909c"
YELLOW = "#fff176"

# ── IceCube strings (3 schematic strings) ────────────────────────────────────
string_xs = [2.0, 5.0, 8.0]
dom_ys    = np.linspace(2.5, 9.5, 14)   # DOM positions along string

for sx in string_xs:
    ax.plot([sx, sx], [2.3, 10.0], color=GRAY, lw=1.5, zorder=1, alpha=0.6)
    for dy in dom_ys:
        circle = plt.Circle((sx, dy), 0.12, color=GRAY, zorder=2, alpha=0.5)
        ax.add_patch(circle)

# ── DM-Ice crystal at the bottom ─────────────────────────────────────────────
crystal_x, crystal_y = 5.0, 1.0
crystal = mpatches.FancyBboxPatch(
    (crystal_x - 0.45, crystal_y - 0.35), 0.9, 0.7,
    boxstyle="round,pad=0.05",
    facecolor=CYAN, edgecolor=WHITE, lw=2, zorder=5
)
ax.add_patch(crystal)
ax.text(crystal_x, crystal_y, "NaI", ha="center", va="center",
        fontsize=10, fontweight="bold", color="#0d1117", zorder=6)
ax.text(crystal_x, crystal_y - 0.72, "DM-Ice det1", ha="center", va="top",
        fontsize=8, color=CYAN, zorder=6)

# ── Muon track ────────────────────────────────────────────────────────────────
# Track enters top-right, exits bottom-left, passes through crystal
track_dx, track_dz = -0.38, -1.0   # direction (normalized below)
norm = np.sqrt(track_dx**2 + track_dz**2)
track_dx /= norm; track_dz /= norm

# Parameterise: track passes through (crystal_x, crystal_y) at t=0
t_vals = np.linspace(-6.5, 6.5, 200)
track_xs = crystal_x + track_dx * t_vals
track_ys = crystal_y + track_dz * t_vals

# Clip to plot area
mask = (track_xs > -0.5) & (track_xs < 11) & (track_ys > -0.5) & (track_ys < 12)
ax.plot(track_xs[mask], track_ys[mask], color=RED, lw=2.5, zorder=4, alpha=0.9)

# Arrow showing muon direction
mid = len(t_vals) // 2 + 30
ax.annotate("", xy=(track_xs[mid+8], track_ys[mid+8]),
            xytext=(track_xs[mid], track_ys[mid]),
            arrowprops=dict(arrowstyle="-|>", color=RED, lw=2.0), zorder=5)

# Muon label
ax.text(track_xs[mid-10] + 0.3, track_ys[mid-10] + 0.3, "μ",
        fontsize=18, color=RED, fontstyle="italic", fontweight="bold", zorder=6)

# ── Light up DOMs that the muon passed near ──────────────────────────────────
for sx in string_xs:
    for dy in dom_ys:
        # Distance from DOM to track
        r = np.array([sx - crystal_x, dy - crystal_y])
        t_along = np.dot(r, [track_dx, track_dz])
        d_perp = np.sqrt(max(0, np.dot(r,r) - t_along**2))
        if d_perp < 1.2 and 0 < t_along < 8:
            # Hit DOM — glow
            glow = plt.Circle((sx, dy), 0.22, color=BLUE, alpha=0.4, zorder=3)
            ax.add_patch(glow)
            dot = plt.Circle((sx, dy), 0.12, color=BLUE, zorder=4)
            ax.add_patch(dot)
            # Timing label (earlier hits higher up)
            t_label = f"t={t_along:.0f}ns" if d_perp < 0.4 else ""

# ── DM-Ice crystal glowing ────────────────────────────────────────────────────
glow_dm = plt.Circle((crystal_x, crystal_y), 0.7, color=CYAN, alpha=0.15, zorder=3)
ax.add_patch(glow_dm)

# ── d_perp = 0 annotation ────────────────────────────────────────────────────
ax.annotate("d⊥ ≈ 0\n(direct transit)",
            xy=(crystal_x, crystal_y),
            xytext=(crystal_x + 2.2, crystal_y + 1.2),
            fontsize=9, color=CYAN,
            arrowprops=dict(arrowstyle="->", color=CYAN, lw=1.2),
            ha="left", va="center")

# ── IceCube label ─────────────────────────────────────────────────────────────
ax.text(9.5, 9.8, "IceCube\nDOMs", ha="center", va="top",
        fontsize=9, color=BLUE, alpha=0.8)

# ── Timing sequence annotation ────────────────────────────────────────────────
ax.text(0.2, 10.8,
        "Timing sequence:",
        fontsize=10, color=WHITE, fontweight="bold")
ax.text(0.2, 10.2,
        "① Muon transits crystal  →  t_geo",
        fontsize=9, color=WHITE)
ax.text(0.2, 9.65,
        "② NaI scintillation  →  +250 ns (decay constant)",
        fontsize=9, color=ORANGE)
ax.text(0.2, 9.1,
        "③ PMT fires  →  t_observed = t_geo + 280 ns ± 81 ns",
        fontsize=9, color=CYAN)

# ── Legend ────────────────────────────────────────────────────────────────────
legend_elements = [
    Line2D([0], [0], color=RED,  lw=2.5, label="Muon track"),
    Line2D([0], [0], color=BLUE, lw=0, marker="o", ms=8, label="IceCube DOM (hit)"),
    mpatches.Patch(facecolor=CYAN, edgecolor=WHITE, label="DM-Ice NaI crystal"),
]
ax.legend(handles=legend_elements, loc="lower right",
          fontsize=8, framealpha=0.2, labelcolor=WHITE,
          facecolor="#1a1a2e", edgecolor=GRAY)

ax.set_title("DM-Ice Timing Geometry — NaI Direct Transit",
             fontsize=13, color=WHITE, pad=12)

plt.tight_layout()
path1 = os.path.join(OUT, "dmice_timing_slide1_geometry.png")
fig.savefig(path1, dpi=150, bbox_inches="tight", facecolor=fig.get_facecolor())
plt.close()
print(f"Saved: {path1}")


# ══════════════════════════════════════════════════════════════════════════════
# SLIDE 2 — Timing chain and residual distribution
# ══════════════════════════════════════════════════════════════════════════════

fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
fig.patch.set_facecolor("white")

# ── Panel A: Timeline diagram ─────────────────────────────────────────────────
ax = axes[0]
ax.set_xlim(-50, 600)
ax.set_ylim(-1, 5)
ax.axis("off")

# Timeline arrow
ax.annotate("", xy=(580, 2), xytext=(-30, 2),
            arrowprops=dict(arrowstyle="-|>", color="black", lw=2))
ax.text(590, 2, "time", va="center", fontsize=10)

# t_geo marker
ax.plot([0, 0], [1.5, 2.5], color="navy", lw=2.5)
ax.text(0, 1.1, "$t_{geo}$\n(muon at crystal)", ha="center", va="top",
        fontsize=10, color="navy", fontweight="bold")
ax.text(0, 2.7, "muon\ntransit", ha="center", va="bottom", fontsize=8, color="navy")

# Scintillation region
ax.axvspan(0, 280, alpha=0.08, color="orange")
ax.annotate("", xy=(280, 3.5), xytext=(0, 3.5),
            arrowprops=dict(arrowstyle="<->", color="orange", lw=1.5))
ax.text(140, 3.7, "NaI scintillation\n+ electronics\n≈ 280 ns",
        ha="center", va="bottom", fontsize=9, color="darkorange")

# μ = 280 ns marker
ax.plot([280, 280], [1.5, 2.5], color="red", lw=2.5, ls="--")
ax.text(280, 1.1, "$\\mu = +280$ ns\n(mean delay)", ha="center", va="top",
        fontsize=10, color="red", fontweight="bold")

# σ band
ax.axvspan(280 - 81, 280 + 81, alpha=0.12, color="red")
ax.annotate("", xy=(280+81, 2.9), xytext=(280-81, 2.9),
            arrowprops=dict(arrowstyle="<->", color="red", lw=1.5))
ax.text(280, 3.05, "$\\sigma = 81$ ns", ha="center", va="bottom",
        fontsize=9, color="red")

# t_observed cloud
ax.plot([280], [2], "o", ms=12, color="cyan", zorder=5)
ax.text(280, 2.55, "$t_{obs}$\n(PMT fires)", ha="center", va="bottom",
        fontsize=9, color="teal")

# Residual arrow
ax.annotate("", xy=(280, 2), xytext=(0, 2),
            arrowprops=dict(arrowstyle="<->", color="purple", lw=1.5,
                            connectionstyle="arc3,rad=0.3"))
ax.text(140, 0.3, r"$\Delta t = t_{obs} - t_{geo}$", ha="center", va="center",
        fontsize=10, color="purple", style="italic")

ax.set_title("A — Timing Chain", fontsize=12, fontweight="bold")

# ── Panel B: Residual distribution ───────────────────────────────────────────
ax = axes[1]

dt = np.linspace(-200, 700, 500)
pdf = np.exp(-0.5 * ((dt - MU_NS) / SIGMA_NS)**2) / (SIGMA_NS * np.sqrt(2*np.pi))

# Fill under curve
ax.fill_between(dt, pdf, alpha=0.25, color="red", label=r"$p(\Delta t)$")
ax.plot(dt, pdf, color="red", lw=2.5)

# Mark μ and σ
ax.axvline(MU_NS, color="darkred", ls="--", lw=1.8, label=f"μ = +{MU_NS:.0f} ns")
ax.axvspan(MU_NS - SIGMA_NS, MU_NS + SIGMA_NS, alpha=0.15, color="red",
           label=f"±σ = {SIGMA_NS:.0f} ns")

# Mark t_geo at 0
ax.axvline(0, color="navy", ls="-", lw=2.0, label="$t_{geo}$ (geometric prediction)")

# Annotate peak
ax.annotate(f"Peak at +{MU_NS:.0f} ns\n(NaI scintillation\n+ electronics)",
            xy=(MU_NS, pdf.max()),
            xytext=(MU_NS + 120, pdf.max() * 0.85),
            fontsize=9,
            arrowprops=dict(arrowstyle="->", color="gray", lw=1.2))

ax.set_xlabel("$\\Delta t = t_{obs} - t_{geo}$ (ns)", fontsize=11)
ax.set_ylabel("Probability density", fontsize=11)
ax.set_title("B — Timing Residual Distribution\n(real DM-Ice coincidence data, 2012–2019)",
             fontsize=11, fontweight="bold")
ax.legend(fontsize=9, loc="upper left")
ax.grid(True, alpha=0.25)
ax.set_xlim(-200, 700)

# ── Pivot equation box ────────────────────────────────────────────────────────
eq_text = (
    "Pivot correction:\n"
    r"$t_{corrected} = t_{obs} - \mu$" "\n"
    r"$= t_{obs} - 280\,\mathrm{ns}$" "\n\n"
    "Likelihood term:\n"
    r"$\log L_{DM} = -\frac{1}{2}\left(\frac{t_0 - t_{corrected}}{\sigma}\right)^2$"
)
ax.text(0.97, 0.97, eq_text, transform=ax.transAxes,
        ha="right", va="top", fontsize=8.5,
        bbox=dict(boxstyle="round,pad=0.5", facecolor="lightyellow",
                  edgecolor="gray", alpha=0.9),
        family="monospace")

plt.suptitle("DM-Ice NaI Scintillator Timing Model", fontsize=14, fontweight="bold", y=1.01)
plt.tight_layout()
path2 = os.path.join(OUT, "dmice_timing_slide2_timing.png")
fig.savefig(path2, dpi=150, bbox_inches="tight")
plt.close()
print(f"Saved: {path2}")

print("\nDone. Two diagram slides saved.")
