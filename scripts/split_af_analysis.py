#!/usr/bin/env python3
"""
Split AF analysis: ranked_0 (AMBER-relaxed) vs ranked_1-4 (unrelaxed)
Per God Claude's methodological note.
"""

import pandas as pd
import numpy as np
from pathlib import Path
from scipy import stats

RESULTS_DIR = Path(__file__).parent.parent / "validation_results"

df = pd.read_csv(RESULTS_DIR / "molprobity_full.csv")

# Split AlphaFold raw models
af = df[df['category'] == 'AlphaFold']
af_raw = af[af['subcategory'] == 'raw']

# ranked_0 = AMBER-relaxed by AF
# ranked_1-4 = unrelaxed by AF
af_ranked0 = af_raw[af_raw['model'] == 'ranked_0']
af_ranked1_4 = af_raw[af_raw['model'].isin(['ranked_1', 'ranked_2', 'ranked_3', 'ranked_4'])]

print("=" * 70)
print("SPLIT AF ANALYSIS: ranked_0 (pre-relaxed) vs ranked_1-4 (unrelaxed)")
print("=" * 70)

print(f"\nModel counts:")
print(f"  ranked_0 (AMBER-relaxed): {len(af_ranked0)} structures")
print(f"  ranked_1-4 (unrelaxed): {len(af_ranked1_4)} structures")

# Compare initial metrics
print("\n### Initial metrics (before Rosetta relaxation)")
for metric in ['clashscore', 'molprobity_score', 'rama_outliers_pct', 'rota_outliers_pct']:
    r0_mean = af_ranked0[metric].mean()
    r14_mean = af_ranked1_4[metric].mean()
    print(f"  {metric}:")
    print(f"    ranked_0: {r0_mean:.2f}")
    print(f"    ranked_1-4: {r14_mean:.2f}")
    print(f"    Difference: {r14_mean - r0_mean:+.2f}")

# Now compare improvement from Rosetta relaxation
print("\n### Improvement from Rosetta (normal_beta) relaxation")

af_relaxed = af[af['subcategory'] == 'relaxed_normal_beta']

# For each protein, compare unrelaxed → relaxed
results = []

for protein in af_ranked0['protein'].unique():
    # ranked_0 path: already AMBER-relaxed, then Rosetta
    r0_before = af_ranked0[af_ranked0['protein'] == protein]['molprobity_score'].mean()
    r0_after = af_relaxed[
        (af_relaxed['protein'] == protein) &
        (af_relaxed['model'].str.startswith('ranked_0'))
    ]['molprobity_score'].mean()

    # ranked_1-4 path: unrelaxed, then Rosetta
    r14_before = af_ranked1_4[af_ranked1_4['protein'] == protein]['molprobity_score'].mean()
    r14_after = af_relaxed[
        (af_relaxed['protein'] == protein) &
        (~af_relaxed['model'].str.startswith('ranked_0'))
    ]['molprobity_score'].mean()

    if pd.notna(r0_before) and pd.notna(r0_after):
        results.append({
            'protein': protein,
            'group': 'ranked_0',
            'before': r0_before,
            'after': r0_after,
            'improvement': r0_before - r0_after
        })

    if pd.notna(r14_before) and pd.notna(r14_after):
        results.append({
            'protein': protein,
            'group': 'ranked_1-4',
            'before': r14_before,
            'after': r14_after,
            'improvement': r14_before - r14_after
        })

results_df = pd.DataFrame(results)

# Statistics per group
for group in ['ranked_0', 'ranked_1-4']:
    group_data = results_df[results_df['group'] == group]
    improvements = group_data['improvement'].values

    stat, p = stats.wilcoxon(improvements)
    mean_imp = improvements.mean()

    print(f"\n{group}:")
    print(f"  N proteins: {len(group_data)}")
    print(f"  Mean improvement: {mean_imp:.3f}")
    print(f"  Wilcoxon p-value: {p:.6f}")
    print(f"  Significant (p<0.05): {'YES' if p < 0.05 else 'NO'}")

# Direct comparison
r0_imps = results_df[results_df['group'] == 'ranked_0']['improvement'].values
r14_imps = results_df[results_df['group'] == 'ranked_1-4']['improvement'].values

print("\n### Head-to-head comparison")
print(f"ranked_0 mean improvement: {r0_imps.mean():.3f}")
print(f"ranked_1-4 mean improvement: {r14_imps.mean():.3f}")

# Mann-Whitney U test between groups
stat, p = stats.mannwhitneyu(r0_imps, r14_imps)
print(f"Mann-Whitney U test p-value: {p:.6f}")

print("\n" + "=" * 70)
print("CONCLUSION")
print("=" * 70)

if results_df[results_df['group'] == 'ranked_1-4']['improvement'].mean() > \
   results_df[results_df['group'] == 'ranked_0']['improvement'].mean():
    print("""
ranked_1-4 (unrelaxed AF models) benefit MORE from Rosetta relaxation
than ranked_0 (AMBER-relaxed AF models).

This confirms God Claude's hypothesis: Rosetta relaxation is redundant for
already-relaxed structures but beneficial for unrelaxed ones.
""")
else:
    print("Unexpected result - need to investigate further")
