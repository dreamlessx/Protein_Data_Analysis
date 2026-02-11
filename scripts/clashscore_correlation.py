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
        relaxed = cat_df[cat_df['subcategory'] == 'relaxed_normal_beta']
        for protein in baseline['protein'].unique():
            base_clash = baseline[baseline['protein'] == protein]['clashscore'].mean()
            rel_data = relaxed[relaxed['protein'] == protein]['clashscore']
            rel_clash = rel_data.mean()
            change = rel_clash - base_clash
            # Range: min and max change across replicates
            change_min = rel_data.min() - base_clash
            change_max = rel_data.max() - base_clash
            if pd.notna(base_clash) and pd.notna(rel_clash):
                data.append({
                    'protein': protein,
                    'category': category,
                    'af_type': None,
                    'initial_clashscore': base_clash,
                    'final_clashscore': rel_clash,
                    'change': change,
                    'change_lo': change - change_min,  # error bar down
                    'change_hi': change_max - change,  # error bar up
                    'improved': rel_clash < base_clash
                })
    elif category == 'AlphaFold':
        # Split by ranked_0 (AMBER-relaxed) vs ranked_1-4 (unrelaxed)
        baseline = cat_df[cat_df['subcategory'] == 'raw']
        relaxed = cat_df[cat_df['subcategory'] == 'relaxed_normal_beta']
        for protein in baseline['protein'].unique():
            # ranked_0 = AMBER-relaxed
            base_r0 = baseline[(baseline['protein'] == protein) & (baseline['model'] == 'ranked_0')]
            rel_r0 = relaxed[(relaxed['protein'] == protein) & (relaxed['model'].str.startswith('ranked_0'))]
            if len(base_r0) and len(rel_r0):
                base_val = base_r0['clashscore'].mean()
                rel_mean = rel_r0['clashscore'].mean()
                change = rel_mean - base_val
                change_min = rel_r0['clashscore'].min() - base_val
                change_max = rel_r0['clashscore'].max() - base_val
                data.append({
                    'protein': protein,
                    'category': category,
                    'af_type': 'ranked_0 (AMBER)',
                    'initial_clashscore': base_val,
                    'final_clashscore': rel_mean,
                    'change': change,
                    'change_lo': change - change_min,
                    'change_hi': change_max - change,
                    'improved': rel_mean < base_val
                })
            # ranked_1-4 = unrelaxed
            base_rest = baseline[(baseline['protein'] == protein) & (~baseline['model'].isin(['ranked_0']))]
            rel_rest = relaxed[(relaxed['protein'] == protein) & (~relaxed['model'].str.startswith('ranked_0'))]
            if len(base_rest) and len(rel_rest):
                base_val = base_rest['clashscore'].mean()
                rel_mean = rel_rest['clashscore'].mean()
                change = rel_mean - base_val
                change_min = rel_rest['clashscore'].min() - base_val
                change_max = rel_rest['clashscore'].max() - base_val
                data.append({
                    'protein': protein,
                    'category': category,
                    'af_type': 'ranked_1-4 (unrelaxed)',
                    'initial_clashscore': base_val,
                    'final_clashscore': rel_mean,
                    'change': change,
                    'change_lo': change - change_min,
                    'change_hi': change_max - change,
                    'improved': rel_mean < base_val
                })
    else:  # Boltz
        baseline = cat_df[cat_df['subcategory'] == 'raw']
        relaxed = cat_df[cat_df['subcategory'] == 'relaxed_normal_beta']
        for protein in baseline['protein'].unique():
            base_clash = baseline[baseline['protein'] == protein]['clashscore'].mean()
            rel_data = relaxed[relaxed['protein'] == protein]['clashscore']
            rel_clash = rel_data.mean()
            change = rel_clash - base_clash
            change_min = rel_data.min() - base_clash
            change_max = rel_data.max() - base_clash
            if pd.notna(base_clash) and pd.notna(rel_clash):
                data.append({
                    'protein': protein,
                    'category': category,
                    'af_type': None,
                    'initial_clashscore': base_clash,
                    'final_clashscore': rel_clash,
                    'change': change,
                    'change_lo': change - change_min,
                    'change_hi': change_max - change,
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

    # For AlphaFold, show split
    if cat == 'AlphaFold':
        for af_type in ['ranked_0 (AMBER)', 'ranked_1-4 (unrelaxed)']:
            sub = cat_data[cat_data['af_type'] == af_type]
            if len(sub) > 1:
                r2, p2 = stats.pearsonr(sub['initial_clashscore'], sub['change'])
                print(f"    {af_type}:")
                print(f"      Improved: {sub['improved'].sum()}/{len(sub)}")
                print(f"      Mean initial: {sub['initial_clashscore'].mean():.1f}")
                print(f"      Mean change: {sub['change'].mean():+.1f}")

# Plot
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

colors = {'Experimental': 'green', 'AlphaFold': 'blue', 'Boltz': 'red'}

for i, cat in enumerate(['Experimental', 'AlphaFold', 'Boltz']):
    ax = axes[i]
    cat_data = data_df[data_df['category'] == cat]

    if cat == 'AlphaFold':
        # Split into ranked_0 (AMBER) vs ranked_1-4 (unrelaxed)
        r0_data = cat_data[cat_data['af_type'] == 'ranked_0 (AMBER)']
        rest_data = cat_data[cat_data['af_type'] == 'ranked_1-4 (unrelaxed)']

        # ranked_0 = black circles, ranked_1-4 = dark blue circles
        ax.errorbar(r0_data['initial_clashscore'], r0_data['change'],
                    yerr=[r0_data['change_lo'], r0_data['change_hi']],
                    fmt='o', c='black', alpha=0.7,
                    markersize=8, label='ranked_0 (AMBER)', capsize=3)
        ax.errorbar(rest_data['initial_clashscore'], rest_data['change'],
                    yerr=[rest_data['change_lo'], rest_data['change_hi']],
                    fmt='o', c='darkblue', alpha=0.7,
                    markersize=8, label='ranked_1-4 (unrelaxed)', capsize=3)

        ax.legend(loc='upper right', fontsize=8)
    else:
        ax.errorbar(cat_data['initial_clashscore'], cat_data['change'],
                    yerr=[cat_data['change_lo'], cat_data['change_hi']],
                    fmt='o', c=colors[cat], alpha=0.7,
                    markersize=8, capsize=3)

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
