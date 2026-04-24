#!/usr/bin/env python3
"""
merge_phase1_results.py
Merges per-run sim_linefit_comparison CSVs into one and makes a combined plot.

Usage: python3 merge_phase1_results.py <phase1_output_dir> <output_plot.png>
"""

import sys
import os
import glob
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    if len(sys.argv) < 3:
        print("Usage: python3 merge_phase1_results.py <phase1_output_dir> <output_plot.png>")
        sys.exit(1)

    outdir = sys.argv[1]
    plot_path = sys.argv[2]

    csv_files = sorted(glob.glob(os.path.join(outdir, "run_*_results.csv")))
    if not csv_files:
        print("ERROR: no run_*_results.csv files found in {}".format(outdir))
        sys.exit(1)

    print("Found {} run CSVs".format(len(csv_files)))
    frames = []
    for f in csv_files:
        try:
            df = pd.read_csv(f)
            run_id = os.path.basename(f).replace("run_", "").replace("_results.csv", "")
            df['run'] = run_id
            frames.append(df)
        except Exception as e:
            print("  WARNING: could not read {}: {}".format(f, e))

    if not frames:
        print("ERROR: no data loaded")
        sys.exit(1)

    df = pd.concat(frames, ignore_index=True)
    merged_csv = os.path.join(outdir, "phase1_all_runs.csv")
    df.to_csv(merged_csv, index=False)
    print("Merged {} events from {} runs -> {}".format(len(df), len(frames), merged_csv))

    # Summary stats
    has_cfit = df['cfit_ang_err_deg'].notna()
    has_both = has_cfit & df['ic_analytic_ang_err_deg'].notna()

    print("\n════════════════════════════════════════════════════════")
    print("Total events:                  {}".format(len(df)))
    print("With pivot LineFit:            {}".format(has_cfit.sum()))
    if has_cfit.any():
        improved = df.loc[has_both, 'cfit_ang_err_deg'] < df.loc[has_both, 'ic_analytic_ang_err_deg']
        print("IC-only median error:          {:.2f} deg".format(
            df.loc[has_both, 'ic_analytic_ang_err_deg'].median()))
        print("DM-Ice Pivot median error:     {:.2f} deg".format(
            df.loc[has_cfit, 'cfit_ang_err_deg'].median()))
        print("DM-Ice improves:               {} / {} ({:.1f}%)".format(
            improved.sum(), has_both.sum(), 100 * improved.mean()))
    print("════════════════════════════════════════════════════════")

    # Plot
    ic_err      = df['ic_analytic_ang_err_deg'].dropna()
    cfit_err    = df['cfit_ang_err_deg'].dropna()
    cfit_iter   = df['cfit_iter_ang_err_deg'].dropna() if 'cfit_iter_ang_err_deg' in df else pd.Series([], dtype=float)

    max_err = max(
        ic_err.max() if len(ic_err) else 0,
        cfit_err.max() if len(cfit_err) else 0,
        cfit_iter.max() if len(cfit_iter) else 0,
    )
    bins = np.linspace(0, min(max_err * 1.05, 90), 46)

    fig, ax = plt.subplots(figsize=(9, 6))
    fig.suptitle(
        'Phase 1 validation: angular error vs MC truth\n'
        '{} events, {} runs'.format(len(df), len(frames)),
        fontsize=13)

    if len(ic_err) > 0:
        ax.hist(ic_err, bins=bins, histtype='stepfilled', alpha=0.5,
                color='steelblue', edgecolor='steelblue',
                label='IC-only LineFit  median={:.1f}°'.format(ic_err.median()))
        ax.axvline(ic_err.median(), color='navy', linewidth=2, linestyle='--')

    if len(cfit_err) > 0:
        ax.hist(cfit_err, bins=bins, histtype='stepfilled', alpha=0.4,
                color='tomato', edgecolor='tomato',
                label='DM-Ice Pivot LineFit  median={:.1f}°'.format(cfit_err.median()))
        ax.axvline(cfit_err.median(), color='darkred', linewidth=2, linestyle='--')

    if len(cfit_iter) > 0:
        ax.hist(cfit_iter, bins=bins, histtype='step', linewidth=2,
                color='darkorange',
                label='DM-Ice Pivot Iterative  median={:.1f}°'.format(cfit_iter.median()))
        ax.axvline(cfit_iter.median(), color='darkorange', linewidth=2, linestyle='--')

    ax.set_xlabel('Angular error from MC truth (deg)', fontsize=12)
    ax.set_ylabel('Events', fontsize=12)
    ax.legend(fontsize=11)
    plt.tight_layout()
    plt.savefig(plot_path, dpi=150)
    print("Saved validation plot to {}".format(plot_path))


if __name__ == '__main__':
    main()
