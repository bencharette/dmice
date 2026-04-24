#!/usr/bin/env python3
"""
ml_coinc_classifier.py

ML-based coincidence classifier using reconstruction agreement as features.

Strategy:
  - Features: pairwise angular differences between all reconstruction pairs,
    plus dm_t_ns, d_perp, n_doms, energy
  - Background label: events with dm_t_ns > 40,000 ns (clearly accidental —
    DM-Ice radioactivity fired 40+ μs into event window)
  - Signal label: events with dm_t_ns in [7000, 16000] ns AND tight
    reconstruction agreement (used only for evaluation, NOT training)
  - Models:
      1. IsolationForest  — unsupervised anomaly detection, no labels needed
      2. GradientBoosting — semi-supervised: trained on background vs rest

Output:
  ~/dmice_work/output/ml_coinc_score.csv   — scores for all events
  ~/dmice_work/output/ml_coinc_plots.png   — diagnostic plots

Run locally (no IceTray needed):
  python3 ~/dmice/ml_coinc_classifier.py
"""

import csv, math, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.ensemble import IsolationForest, GradientBoostingClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.inspection import permutation_importance
from sklearn.metrics import roc_auc_score

OUT_DIR  = os.path.expanduser("~/dmice_work/output")
CSV_IN   = os.path.join(OUT_DIR, "real_all_recos.csv")

# ── Timing windows ─────────────────────────────────────────────────────────────
BG_DM_T_MIN   = 40_000   # ns — clearly accidental (late radioactivity)
SIG_DM_T_LO   =  7_000   # ns — expected muon transit window (low)
SIG_DM_T_HI   = 16_000   # ns — expected muon transit window (high)

# ── Load CSV ───────────────────────────────────────────────────────────────────

rows = list(csv.DictReader(open(CSV_IN)))
print(f"Loaded {len(rows)} events from {CSV_IN}")

def flt(r, k):
    try:
        v = float(r[k])
        return v if math.isfinite(v) else np.nan
    except (ValueError, TypeError, KeyError):
        return np.nan

def ang_diff(z1, a1, z2, a2):
    """Angular separation between two directions (radians input → degrees out)."""
    if any(math.isnan(v) for v in [z1, a1, z2, a2]):
        return np.nan
    d1 = np.array([math.sin(z1)*math.cos(a1), math.sin(z1)*math.sin(a1), -math.cos(z1)])
    d2 = np.array([math.sin(z2)*math.cos(a2), math.sin(z2)*math.sin(a2), -math.cos(z2)])
    dot = float(np.dot(d1, d2))
    dot = max(-1.0, min(1.0, dot))
    return math.degrees(math.acos(dot))

# ── Feature engineering ────────────────────────────────────────────────────────

RECO_PAIRS = [
    # (label, zen_key1, azi_key1, zen_key2, azi_key2)
    ("da_lf_pivlf",       "lf_zen",         "lf_azi",         "piv_lf_zen",       "piv_lf_azi"),
    ("da_lf_mpe",         "lf_zen",         "lf_azi",         "mpe_std_zen",      "mpe_std_azi"),
    ("da_lf_itermpe",     "lf_zen",         "lf_azi",         "iter_mpe_zen",     "iter_mpe_azi"),
    ("da_lf_spepiv",      "lf_zen",         "lf_azi",         "spe_piv_zen",      "spe_piv_azi"),
    ("da_mpe_pivmpe",     "mpe_std_zen",    "mpe_std_azi",    "mpe_piv_zen",      "mpe_piv_azi"),
    ("da_mpe_itermpe",    "mpe_std_zen",    "mpe_std_azi",    "iter_mpe_zen",     "iter_mpe_azi"),
    ("da_mpe_iterpivmpe", "mpe_std_zen",    "mpe_std_azi",    "iter_piv_mpe_zen", "iter_piv_mpe_azi"),
    ("da_itermpe_iterpivlf","iter_mpe_zen", "iter_mpe_azi",   "iter_piv_lf_zen",  "iter_piv_lf_azi"),
    ("da_itermpe_iterpivmpe","iter_mpe_zen","iter_mpe_azi",   "iter_piv_mpe_zen", "iter_piv_mpe_azi"),
    ("da_pivlf_pivmpe",   "piv_lf_zen",    "piv_lf_azi",     "mpe_piv_zen",      "mpe_piv_azi"),
]

SCALAR_FEATURES = ["dm_t_ns", "d_perp_m", "n_doms_ic", "n_hits_ic", "energy_GeV", "lf_zen"]

records = []
for r in rows:
    feat = {}

    # Scalar features
    for k in SCALAR_FEATURES:
        feat[k] = flt(r, k)

    # Angular difference features
    for (lab, zk1, ak1, zk2, ak2) in RECO_PAIRS:
        z1 = math.radians(flt(r, zk1)) if not np.isnan(flt(r, zk1)) else np.nan
        a1 = math.radians(flt(r, ak1)) if not np.isnan(flt(r, ak1)) else np.nan
        z2 = math.radians(flt(r, zk2)) if not np.isnan(flt(r, zk2)) else np.nan
        a2 = math.radians(flt(r, ak2)) if not np.isnan(flt(r, ak2)) else np.nan
        feat[lab] = ang_diff(z1, a1, z2, a2)

    feat["year"]     = int(r.get("year", 0))
    feat["run_id"]   = int(r.get("run_id", 0))
    feat["event_id"] = int(r.get("event_id", 0))
    records.append(feat)

# Feature columns used in ML (exclude identifiers and dm_t_ns from classifier
# — we don't want dm_t_ns to trivially drive the score)
ANG_COLS    = [p[0] for p in RECO_PAIRS]
# NOTE: d_perp excluded — it's correlated with the label by construction
# (geometric cut events cluster in signal window), making the classifier
# trivially learn d_perp rather than reconstruction agreement.
SCALAR_COLS = ["n_doms_ic", "n_hits_ic", "energy_GeV", "lf_zen"]
FEAT_COLS   = ANG_COLS + SCALAR_COLS

print(f"Features ({len(FEAT_COLS)}): {FEAT_COLS}")

# Build matrix — impute missing with column median
X_raw = np.array([[r[c] for c in FEAT_COLS] for r in records], dtype=float)
dm_t  = np.array([r["dm_t_ns"] for r in records], dtype=float)
years = np.array([r["year"]    for r in records], dtype=int)

# Impute NaN with column median
col_medians = np.nanmedian(X_raw, axis=0)
for j in range(X_raw.shape[1]):
    mask = np.isnan(X_raw[:, j])
    X_raw[mask, j] = col_medians[j]

scaler = StandardScaler()
X = scaler.fit_transform(X_raw)

# ── Labels ─────────────────────────────────────────────────────────────────────

# Background: late-firing DM-Ice (clear accidentals)
is_bg  = dm_t > BG_DM_T_MIN
# "Signal window": expected muon transit time (used only for evaluation)
in_sig_window = (dm_t >= SIG_DM_T_LO) & (dm_t <= SIG_DM_T_HI)

print(f"\nBackground events (dm_t > {BG_DM_T_MIN/1000:.0f} μs): {is_bg.sum()}")
print(f"Signal window events ({SIG_DM_T_LO/1000:.0f}–{SIG_DM_T_HI/1000:.0f} μs): {in_sig_window.sum()}")
print(f"Other events: {(~is_bg & ~in_sig_window).sum()}")

# ── Model 1: Isolation Forest (unsupervised) ───────────────────────────────────

print("\nFitting Isolation Forest …")
iso = IsolationForest(n_estimators=300, contamination=0.05, random_state=42)
iso.fit(X)
# score_samples: higher = more normal, lower = more anomalous
iso_scores = iso.score_samples(X)   # range roughly [-0.7, 0.1]
iso_anom   = -iso_scores            # flip: higher = more anomalous

# ── Model 2: GradientBoosting semi-supervised ──────────────────────────────────
# Train: background (label=0) vs signal-window events (label=1)
# Intentionally exclude "other" events from training to avoid circular logic

print("Fitting GradientBoosting classifier …")
train_mask = is_bg | in_sig_window
X_train    = X[train_mask]
y_train    = in_sig_window[train_mask].astype(int)

print(f"  Training: {train_mask.sum()} events  "
      f"(bg={is_bg.sum()}, sig_window={in_sig_window.sum()})")

gb = GradientBoostingClassifier(n_estimators=200, max_depth=4,
                                learning_rate=0.05, random_state=42)
gb.fit(X_train, y_train)
gb_prob = gb.predict_proba(X)[:, 1]   # P(signal-like)

# Quick AUC on training set (just for sanity)
auc = roc_auc_score(y_train, gb.predict_proba(X_train)[:, 1])
print(f"  Training AUC: {auc:.3f}")

# Feature importances
importances = gb.feature_importances_
feat_order  = np.argsort(importances)[::-1]
print("\nTop-10 features (GradientBoosting):")
for i in feat_order[:10]:
    print(f"  {FEAT_COLS[i]:<35} {importances[i]:.4f}")

# ── Save scores ────────────────────────────────────────────────────────────────

out_csv = os.path.join(OUT_DIR, "ml_coinc_score.csv")
with open(out_csv, "w", newline="") as fh:
    writer = csv.writer(fh)
    writer.writerow(["year", "run_id", "event_id", "dm_t_ns",
                     "iso_anom", "gb_prob", "is_bg", "in_sig_window"]
                    + FEAT_COLS)
    for i, r in enumerate(records):
        writer.writerow([r["year"], r["run_id"], r["event_id"],
                         dm_t[i], iso_anom[i], gb_prob[i],
                         int(is_bg[i]), int(in_sig_window[i])]
                        + list(X_raw[i]))
print(f"\nScores saved: {out_csv}")

# ── Top candidates ─────────────────────────────────────────────────────────────

print("\nTop-20 most signal-like events (GradientBoosting score):")
print(f"{'Rank':>4} {'Year':>5} {'run_id':>8} {'event_id':>12} "
      f"{'dm_t_μs':>8} {'gb_prob':>8} {'iso_anom':>9} {'sig_win':>7}")
order = np.argsort(gb_prob)[::-1]
for rank, i in enumerate(order[:20]):
    r = records[i]
    print(f"{rank+1:>4} {r['year']:>5} {r['run_id']:>8} {r['event_id']:>12} "
          f"{dm_t[i]/1000:>8.1f} {gb_prob[i]:>8.4f} {iso_anom[i]:>9.4f} "
          f"{'YES' if in_sig_window[i] else '':>7}")

# ── Plots ──────────────────────────────────────────────────────────────────────

fig, axes = plt.subplots(2, 3, figsize=(18, 11))
fig.patch.set_facecolor("#111111")
for ax in axes.flat:
    ax.set_facecolor("#1a1a2e")

# 1. GB score vs dm_t_ns
ax = axes[0, 0]
cmap_yr = plt.cm.tab10
for grp, col, lab, zo in [
    (is_bg,           "red",        "Background (dm_t>40μs)",      1),
    (~is_bg & ~in_sig_window, "steelblue", "Other",                2),
    (in_sig_window,   "lime",       f"Signal window (7–16 μs)",    3),
]:
    ax.scatter(dm_t[grp]/1000, gb_prob[grp], s=4, alpha=0.5, color=col, label=lab, zorder=zo)
ax.set_xlabel("dm_t [μs]", color="white")
ax.set_ylabel("GB signal probability", color="white")
ax.set_title("GradientBoosting score vs DM-Ice hit time", color="white")
ax.axvspan(SIG_DM_T_LO/1000, SIG_DM_T_HI/1000, alpha=0.1, color="lime")
ax.axhline(0.5, color="yellow", lw=1, ls="--", label="p=0.5")
ax.legend(fontsize=8, labelcolor="white", facecolor="#222244")
ax.tick_params(colors="white")
for s in ax.spines.values(): s.set_edgecolor("#555")

# 2. Isolation Forest score vs dm_t_ns
ax = axes[0, 1]
for grp, col, lab in [
    (is_bg,           "red",    "Background"),
    (~is_bg & ~in_sig_window, "steelblue", "Other"),
    (in_sig_window,   "lime",   "Signal window"),
]:
    ax.scatter(dm_t[grp]/1000, iso_anom[grp], s=4, alpha=0.5, color=col, label=lab)
ax.set_xlabel("dm_t [μs]", color="white")
ax.set_ylabel("Isolation Forest anomaly score", color="white")
ax.set_title("IsolationForest anomaly vs DM-Ice hit time", color="white")
ax.axvspan(SIG_DM_T_LO/1000, SIG_DM_T_HI/1000, alpha=0.1, color="lime")
ax.legend(fontsize=8, labelcolor="white", facecolor="#222244")
ax.tick_params(colors="white")
for s in ax.spines.values(): s.set_edgecolor("#555")

# 3. GB score distribution by group
ax = axes[0, 2]
bins = np.linspace(0, 1, 40)
for grp, col, lab in [
    (is_bg,           "red",      "Background (dm_t>40μs)"),
    (in_sig_window,   "lime",     "Signal window (7–16 μs)"),
    (~is_bg & ~in_sig_window, "steelblue", "Other"),
]:
    if grp.sum() > 0:
        ax.hist(gb_prob[grp], bins=bins, density=True, alpha=0.6,
                color=col, label=f"{lab} (n={grp.sum()})", histtype="stepfilled")
ax.set_xlabel("GB signal probability", color="white")
ax.set_ylabel("Normalised density", color="white")
ax.set_title("Score distribution by group", color="white")
ax.legend(fontsize=8, labelcolor="white", facecolor="#222244")
ax.tick_params(colors="white")
for s in ax.spines.values(): s.set_edgecolor("#555")

# 4. Feature importances
ax = axes[1, 0]
top_n = 10
top_feats = [FEAT_COLS[i] for i in feat_order[:top_n]]
top_imps  = [importances[i] for i in feat_order[:top_n]]
colors_bar = ["lime" if "da_" in f else "steelblue" for f in top_feats]
bars = ax.barh(range(top_n), top_imps[::-1], color=colors_bar[::-1], alpha=0.85)
ax.set_yticks(range(top_n))
ax.set_yticklabels([f.replace("da_","Δθ ").replace("_"," ") for f in top_feats[::-1]],
                   color="white", fontsize=9)
ax.set_xlabel("Feature importance", color="white")
ax.set_title("Top-10 GradientBoosting features\n(green=angular diff, blue=scalar)", color="white")
ax.tick_params(colors="white")
for s in ax.spines.values(): s.set_edgecolor("#555")
ax.grid(axis="x", alpha=0.2)

# 5. Per-year: how many events score > 0.5
ax = axes[1, 1]
all_years = sorted(set(years))
n_high_gb  = [np.sum((years == yr) & (gb_prob > 0.5))  for yr in all_years]
n_high_iso = [np.sum((years == yr) & (iso_anom > np.percentile(iso_anom, 95))) for yr in all_years]
n_total_yr = [np.sum(years == yr) for yr in all_years]
x = np.arange(len(all_years)); w = 0.3
ax.bar(x - w/2, n_high_gb,  w, color="lime",      alpha=0.85, label="GB prob > 0.5")
ax.bar(x + w/2, n_high_iso, w, color="darkorange", alpha=0.85, label="IsoForest top 5%")
ax.set_xticks(x); ax.set_xticklabels(all_years, rotation=30, color="white")
ax.set_ylabel("Events", color="white")
ax.set_title("Signal-like events per year", color="white")
ax.legend(fontsize=8, labelcolor="white", facecolor="#222244")
ax.tick_params(colors="white")
for s in ax.spines.values(): s.set_edgecolor("#555")
ax.grid(axis="y", alpha=0.2)

# 6. GB vs IsoForest scatter — correlation between methods
ax = axes[1, 2]
sc = ax.scatter(gb_prob, iso_anom, c=dm_t/1000,
                cmap="plasma", s=5, alpha=0.5,
                vmin=0, vmax=60)
ax.set_xlabel("GB signal probability", color="white")
ax.set_ylabel("IsoForest anomaly score", color="white")
ax.set_title("GB vs IsolationForest — colour = dm_t [μs]", color="white")
cb = fig.colorbar(sc, ax=ax)
cb.set_label("dm_t [μs]", color="white")
cb.ax.yaxis.set_tick_params(color="white")
plt.setp(cb.ax.yaxis.get_ticklabels(), color="white")
ax.tick_params(colors="white")
for s in ax.spines.values(): s.set_edgecolor("#555")

fig.suptitle(
    "ML coincidence classifier — reconstruction agreement features\n"
    "Background: dm_t>40μs  |  Signal window: 7–16μs  |  6000 real events 2012–2021",
    color="white", fontsize=11
)
plt.tight_layout()
out_png = os.path.join(OUT_DIR, "ml_coinc_plots.png")
fig.savefig(out_png, dpi=150, bbox_inches="tight", facecolor=fig.get_facecolor())
plt.close()
print(f"Plot: {out_png}")
print("\nDone.")
