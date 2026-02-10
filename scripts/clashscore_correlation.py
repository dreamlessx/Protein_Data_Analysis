#!/usr/bin/env python3
"""Check if initial clashscore predicts relaxation outcome."""

import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
from scipy import stats

RESULTS_DIR = Path(__file__).parent.parent / "validation_results"
FIGURES_DIR = Path(__file__).parent.parent / "figures"

df = pd.read_csv(RESULTS_DIR / "molprobity_full.csv")

# Collect paired data: initial clashscore vs change in clashscore
data = []

for category in ['Experimental', 'AlphaFold', 'Boltz']:
    cat_df = df[df['category'] == category]

    if category == 'Experimental':
        baseline = cat_df[cat_df['subcategory'] == 'original']
    else:
        baseline = cat_df[cat_df['subcategory'] == 'raw']

    relaxed = cat_df[cat_df['subcategory'] == 'relaxed_normal_beta']

    for protein in baseline['protein'].unique():
        base_clash = baseline[baseline['protein'] == protein]['clashscore'].mean()
        rel_clash = relaxed[relaxed['protein'] == protein]['clashscore'].mean()

        if pd.notna(base_clash) and pd.notna(rel_clash):
            data.append({
                'protein': protein,
                'category': category,
                'initial_clashscore': base_clash,
                'final_clashscore': rel_clash,
                'change': rel_clash - base_clash,  # negative = improved
                'improved': rel_clash < base_clash
            })

data_df = pd.DataFrame(data)

print("=" * 60)
print("CLASHSCORE CORRELATION ANALYSIS")
print("=" * 60)

# Correlation between initial clashscore and change
r, p = stats.pearsonr(data_df['initial_clashscore'], data_df['change'])
print(f"\nCorrelation (initial vs change): r = {r:.3f}, p = {p:.6f}")

# Split by category
for cat in ['Experimental', 'AlphaFold', 'Boltz']:
    cat_data = data_df[data_df['category'] == cat]
    r, p = stats.pearsonr(cat_data['initial_clashscore'], cat_data['change'])
    improved = cat_data['improved'].sum()
    total = len(cat_data)
    print(f"\n{cat}:")
    print(f"  Correlation: r = {r:.3f}, p = {p:.4f}")
    print(f"  Improved: {improved}/{total} ({100*improved/total:.0f}%)")
    print(f"  Mean initial clashscore: {cat_data['initial_clashscore'].mean():.1f}")
    print(f"  Mean change: {cat_data['change'].mean():+.1f}")

# Plot
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

colors = {'Experimental': 'green', 'AlphaFold': 'blue', 'Boltz': 'red'}

for i, cat in enumerate(['Experimental', 'AlphaFold', 'Boltz']):
    ax = axes[i]
    cat_data = data_df[data_df['category'] == cat]

    ax.scatter(cat_data['initial_clashscore'], cat_data['change'],
               c=colors[cat], alpha=0.7, s=80)

    # Label outliers
    for _, row in cat_data.iterrows():
        if row['change'] > 5:  # Got much worse
            ax.annotate(row['protein'], (row['initial_clashscore'], row['change']),
                        fontsize=8, alpha=0.7)

    ax.axhline(0, color='black', linestyle='--', alpha=0.5)
    ax.set_xlabel('Initial Clashscore')
    ax.set_ylabel('Change in Clashscore (- = better)')
    ax.set_title(f'{cat}')

    # Regression line
    z = np.polyfit(cat_data['initial_clashscore'], cat_data['change'], 1)
    p = np.poly1d(z)
    x_line = np.linspace(cat_data['initial_clashscore'].min(),
                         cat_data['initial_clashscore'].max(), 100)
    ax.plot(x_line, p(x_line), 'k-', alpha=0.3)

plt.suptitle('Initial Clashscore vs Change After Relaxation', fontsize=14)
plt.tight_layout()
plt.savefig(FIGURES_DIR / 'clashscore_correlation.png', dpi=150)
plt.close()
print(f"\nSaved clashscore_correlation.png")

# Summary: structures that got worse
print("\n" + "=" * 60)
print("STRUCTURES WHERE CLASHSCORE INCREASED (GOT WORSE)")
print("=" * 60)
worse = data_df[data_df['change'] > 0].sort_values('change', ascending=False)
print(worse[['protein', 'category', 'initial_clashscore', 'final_clashscore', 'change']].to_string())
