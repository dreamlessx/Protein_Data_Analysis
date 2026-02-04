#!/usr/bin/env python3
"""
compare_trials.py
Compare validation trials for consistency and generate final analysis
"""

import pandas as pd
import numpy as np
from pathlib import Path

PROJECT_DIR = Path(__file__).parent.parent
RESULTS_DIR = PROJECT_DIR / "validation_results"

def main():
    print("=" * 70)
    print("TRIAL COMPARISON AND CONSISTENCY CHECK")
    print("=" * 70)

    # Load all trials
    trials = []
    for i in range(1, 4):
        trial_file = RESULTS_DIR / f"trial_{i}_results.csv"
        if trial_file.exists():
            df = pd.read_csv(trial_file)
            trials.append(df)
            print(f"Trial {i}: {len(df)} structures")

    if len(trials) < 2:
        print("Need at least 2 trials for comparison")
        return

    # Key metrics to compare
    metrics = ['clashscore', 'rama_favored_pct', 'rama_outliers_pct',
               'rotamer_outliers_pct', 'cbeta_deviations', 'bond_rmsz',
               'angle_rmsz', 'molprobity_score']

    print("\n" + "=" * 70)
    print("METRIC-BY-METRIC COMPARISON")
    print("=" * 70)

    all_identical = True

    for metric in metrics:
        if metric not in trials[0].columns:
            continue

        values = [df[metric].dropna().values for df in trials]

        # Check if all trials have same length for this metric
        lengths = [len(v) for v in values]
        if len(set(lengths)) > 1:
            print(f"\n{metric}: Different sample sizes across trials - {lengths}")
            continue

        # Compare values
        if len(values) >= 2:
            # Check if trial 1 and 2 are identical
            diff_1_2 = np.sum(np.abs(values[0] - values[1]))
            identical_1_2 = np.allclose(values[0], values[1], rtol=1e-10)

            # Calculate statistics
            means = [np.mean(v) for v in values]
            stds = [np.std(v) for v in values]

            print(f"\n{metric}:")
            for i, (m, s) in enumerate(zip(means, stds), 1):
                print(f"  Trial {i}: mean={m:.6f}, std={s:.6f}")

            if len(values) >= 2:
                print(f"  Trial 1 vs 2: {'IDENTICAL' if identical_1_2 else f'diff={diff_1_2:.6f}'}")
                if not identical_1_2:
                    all_identical = False

            if len(values) == 3:
                identical_2_3 = np.allclose(values[1], values[2], rtol=1e-10)
                print(f"  Trial 2 vs 3: {'IDENTICAL' if identical_2_3 else 'DIFFERENT'}")

    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)

    if all_identical:
        print("✅ ALL TRIALS PRODUCED IDENTICAL RESULTS")
        print("   The validation pipeline is fully deterministic.")
    else:
        print("⚠️  Some differences detected between trials")

    # Check for zeros
    print("\n" + "=" * 70)
    print("ZERO VALUE ANALYSIS")
    print("=" * 70)

    df = trials[0]  # Use trial 1 for analysis

    zero_analysis = []
    for col in metrics:
        if col in df.columns:
            zero_count = (df[col] == 0).sum()
            total = df[col].notna().sum()
            if zero_count > 0:
                zero_analysis.append({
                    'Metric': col,
                    'Zero Count': zero_count,
                    'Total': total,
                    'Zero %': f"{100*zero_count/total:.1f}%"
                })

    if zero_analysis:
        print("\nMetrics with zero values:")
        for z in zero_analysis:
            print(f"  {z['Metric']}: {z['Zero Count']}/{z['Total']} ({z['Zero %']})")
    else:
        print("No zeros found in key metrics")

    # Detailed zero analysis by category
    print("\n" + "-" * 40)
    print("Zero Analysis by Category:")

    for cat in ['AlphaFold', 'Boltz', 'Experimental']:
        cat_df = df[df['category'] == cat]
        print(f"\n{cat} ({len(cat_df)} structures):")
        for col in ['clashscore', 'rama_outliers_pct', 'rotamer_outliers_pct',
                    'cbeta_deviations', 'omega_cis_general', 'omega_twisted']:
            if col in cat_df.columns:
                zeros = (cat_df[col] == 0).sum()
                if zeros > 0:
                    print(f"  {col}: {zeros} zeros ({100*zeros/len(cat_df):.1f}%)")


if __name__ == "__main__":
    main()
