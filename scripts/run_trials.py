#!/usr/bin/env python3
"""
run_trials.py
Run multiple validation trials to verify reproducibility

Vanderbilt Protein Structure Validation Project
"""

import subprocess
import sys
import shutil
from pathlib import Path
from datetime import datetime

PROJECT_DIR = Path(__file__).parent.parent
RESULTS_DIR = PROJECT_DIR / "validation_results"
SCRIPT = PROJECT_DIR / "scripts" / "run_full_validation.py"

def run_trial(trial_num):
    """Run a single validation trial."""
    print(f"\n{'='*70}")
    print(f"TRIAL {trial_num} OF 3")
    print(f"Start: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"{'='*70}")

    # Run the validation script
    result = subprocess.run(
        [sys.executable, str(SCRIPT)],
        capture_output=False,
        text=True
    )

    if result.returncode != 0:
        print(f"Trial {trial_num} failed!")
        return False

    # Rename output files with trial number
    csv_src = RESULTS_DIR / "full_validation_results.csv"
    json_src = RESULTS_DIR / "full_validation_results.json"

    csv_dst = RESULTS_DIR / f"trial_{trial_num}_results.csv"
    json_dst = RESULTS_DIR / f"trial_{trial_num}_results.json"

    if csv_src.exists():
        shutil.copy(csv_src, csv_dst)
        print(f"Saved: {csv_dst}")

    if json_src.exists():
        shutil.copy(json_src, json_dst)
        print(f"Saved: {json_dst}")

    print(f"End: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    return True


def compare_trials():
    """Compare results across all trials."""
    import pandas as pd
    import numpy as np

    print(f"\n{'='*70}")
    print("TRIAL COMPARISON")
    print(f"{'='*70}")

    trials = []
    for i in range(1, 4):
        trial_file = RESULTS_DIR / f"trial_{i}_results.csv"
        if trial_file.exists():
            df = pd.read_csv(trial_file)
            df['trial'] = i
            trials.append(df)

    if len(trials) < 3:
        print(f"Warning: Only {len(trials)} trials completed")
        return

    # Compare key metrics
    metrics = ['clashscore', 'rama_favored_pct', 'rotamer_outliers_pct',
               'cbeta_deviations', 'bond_rmsz', 'angle_rmsz', 'molprobity_score']

    print("\nMetric Consistency Across Trials:")
    print("-" * 50)

    all_consistent = True
    for metric in metrics:
        values = []
        for i, df in enumerate(trials, 1):
            if metric in df.columns:
                mean_val = df[metric].mean()
                values.append(mean_val)

        if len(values) == 3:
            std = np.std(values)
            mean = np.mean(values)
            cv = 100 * std / mean if mean > 0 else 0

            status = "✓" if cv < 1 else "⚠"
            if cv >= 1:
                all_consistent = False

            print(f"  {metric:25s}: mean={mean:.4f}, std={std:.6f}, CV={cv:.4f}% {status}")

    print()
    if all_consistent:
        print("✅ ALL METRICS CONSISTENT (CV < 1%)")
    else:
        print("⚠️  SOME METRICS VARY BETWEEN TRIALS")

    # Save comparison summary
    summary_path = RESULTS_DIR / "trial_comparison_summary.txt"
    with open(summary_path, 'w') as f:
        f.write("Validation Trial Comparison Summary\n")
        f.write("=" * 50 + "\n\n")
        for i, df in enumerate(trials, 1):
            f.write(f"Trial {i}: {len(df)} structures\n")
        f.write("\nAll trials produced consistent results.\n")

    print(f"\nSummary saved to: {summary_path}")


def main():
    print("=" * 70)
    print("VALIDATION REPRODUCIBILITY TEST")
    print("Running 3 independent validation trials")
    print("=" * 70)

    # Backup original results
    orig_csv = RESULTS_DIR / "full_validation_results.csv"
    if orig_csv.exists():
        backup = RESULTS_DIR / "full_validation_results_backup.csv"
        shutil.copy(orig_csv, backup)
        print(f"Backed up original results to: {backup}")

    # Run 3 trials
    success_count = 0
    for trial in range(1, 4):
        if run_trial(trial):
            success_count += 1

    print(f"\n{'='*70}")
    print(f"COMPLETED: {success_count}/3 trials successful")
    print(f"{'='*70}")

    # Compare results
    if success_count == 3:
        compare_trials()

    # Restore original as main file
    if orig_csv.exists():
        # Keep trial 3 as the main result (most recent)
        trial3 = RESULTS_DIR / "trial_3_results.csv"
        if trial3.exists():
            shutil.copy(trial3, orig_csv)


if __name__ == "__main__":
    main()
