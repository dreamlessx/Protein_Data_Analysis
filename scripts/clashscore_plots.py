#!/usr/bin/env python3
"""Clash score bar plots and correlation scatter plot."""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

df = pd.read_csv('validation_results/molprobity_full.csv')

# Remove outliers
outliers = ['1GHQ', '1F51']
df = df[~df['protein'].isin(outliers)]

# Define before states for each category
before_map = {
    'Experimental': 'original',
    'AlphaFold': 'raw',
    'Boltz': 'raw'
}

relax_protocols = [
    'relaxed_normal_beta',
    'relaxed_normal_ref15',
    'relaxed_dualspace_beta',
    'relaxed_dualspace_ref15',
    'relaxed_cartesian_beta',
    'relaxed_cartesian_ref15'
]

short_names = {
    'relaxed_normal_beta': 'normal_β',
    'relaxed_normal_ref15': 'normal_ref15',
    'relaxed_dualspace_beta': 'dual_β',
    'relaxed_dualspace_ref15': 'dual_ref15',
    'relaxed_cartesian_beta': 'cart_β',
    'relaxed_cartesian_ref15': 'cart_ref15'
}

colors = {
    'before': '#2c3e50',
    'relaxed_normal_beta': '#27ae60',
    'relaxed_normal_ref15': '#2ecc71',
    'relaxed_dualspace_beta': '#3498db',
    'relaxed_dualspace_ref15': '#5dade2',
    'relaxed_cartesian_beta': '#e74c3c',
    'relaxed_cartesian_ref15': '#ec7063'
}

Path('figures').mkdir(exist_ok=True)

# ============ 1. Per-protein bar plots ============
proteins = df['protein'].unique()
categories_expanded = ['Experimental', 'AlphaFold (ranked_0)', 'AlphaFold (ranked_1-4)', 'Boltz']

fig, axes = plt.subplots(4, 1, figsize=(16, 10), sharex=True)
fig.subplots_adjust(hspace=0)

for idx, cat_label in enumerate(categories_expanded):
    ax = axes[idx]

    if cat_label == 'AlphaFold (ranked_0)':
        cat_df = df[df['category'] == 'AlphaFold']
        before_sub = 'raw'
        model_filter = ['ranked_0']
        relax_prefix = ['ranked_0']
    elif cat_label == 'AlphaFold (ranked_1-4)':
        cat_df = df[df['category'] == 'AlphaFold']
        before_sub = 'raw'
        model_filter = ['ranked_1', 'ranked_2', 'ranked_3', 'ranked_4']
        relax_prefix = ['ranked_1', 'ranked_2', 'ranked_3', 'ranked_4']
    else:
        cat = 'Experimental' if cat_label == 'Experimental' else 'Boltz'
        cat_df = df[df['category'] == cat]
        before_sub = before_map[cat]
        model_filter = None
        relax_prefix = None

    x = np.arange(len(proteins))
    width = 0.14
    bar_group_width = width * 6

    # Before values
    before_vals = []
    for p in proteins:
        if model_filter:
            vals = cat_df[(cat_df['protein'] == p) &
                          (cat_df['subcategory'] == before_sub) &
                          (cat_df['model'].isin(model_filter))]['clashscore']
        else:
            vals = cat_df[(cat_df['protein'] == p) & (cat_df['subcategory'] == before_sub)]['clashscore']
        before_vals.append(vals.mean() if len(vals) > 0 else np.nan)

    # After values for each protocol
    for i, proto in enumerate(relax_protocols):
        after_vals = []
        for p in proteins:
            if relax_prefix:
                # Match relaxed models by prefix
                relaxed = cat_df[(cat_df['protein'] == p) & (cat_df['subcategory'] == proto)]
                matched_vals = []
                for prefix in relax_prefix:
                    matched = relaxed[relaxed['model'].str.startswith(prefix + '_')]
                    if len(matched) > 0:
                        matched_vals.extend(matched['clashscore'].tolist())
                after_vals.append(np.mean(matched_vals) if matched_vals else np.nan)
            else:
                vals = cat_df[(cat_df['protein'] == p) & (cat_df['subcategory'] == proto)]['clashscore']
                after_vals.append(vals.mean() if len(vals) > 0 else np.nan)
        ax.bar(x + (i-2.5)*width, after_vals, width, label=short_names[proto], color=colors[proto])

    # Overlay initial values as wide semi-transparent bars (in front)
    ax.bar(x, before_vals, bar_group_width, alpha=0.35, color='black', label='Initial', zorder=10)

    ax.set_ylabel(cat_label, fontsize=9)
    ax.set_ylim(0, 52)

    # Only show legend on first subplot
    if idx == 0:
        ax.legend(loc='upper left', bbox_to_anchor=(1.01, 1), fontsize=7)

# Only label x-axis on bottom subplot
axes[-1].set_xticks(x)
axes[-1].set_xticklabels(proteins, rotation=45, ha='right', fontsize=8)

# Common y-axis label
fig.text(0.01, 0.5, 'Clash Score', va='center', rotation='vertical', fontsize=11)

plt.savefig('figures/clashscore_per_protein.png', dpi=150, bbox_inches='tight')
plt.close()
print('Saved figures/clashscore_per_protein.png')

# ============ 2. Averaged bar plot (horizontal with labels inside) ============
categories = ['Experimental', 'AlphaFold', 'Boltz']
fig, ax = plt.subplots(figsize=(10, 8))

y_labels = ['Experimental', 'AF (ranked_0)', 'AF (ranked_1-4)', 'Boltz']
y = np.arange(len(y_labels))
height = 0.12


def get_af_clashscores(subset_models, subcategory, relax_prefix=None):
    """Get AlphaFold clashscores for specific model subset."""
    af_df = df[df['category'] == 'AlphaFold']
    if subcategory == 'raw':
        return af_df[(af_df['subcategory'] == 'raw') &
                     (af_df['model'].isin(subset_models))]['clashscore'].mean()
    else:
        relaxed = af_df[af_df['subcategory'] == subcategory]
        vals = []
        for prefix in relax_prefix:
            matched = relaxed[relaxed['model'].str.startswith(prefix + '_')]
            vals.extend(matched['clashscore'].tolist())
        return np.mean(vals) if vals else np.nan


# Before averages
before_avgs = [
    df[(df['category'] == 'Experimental') & (df['subcategory'] == 'original')]['clashscore'].mean(),
    get_af_clashscores(['ranked_0'], 'raw'),
    get_af_clashscores(['ranked_1', 'ranked_2', 'ranked_3', 'ranked_4'], 'raw'),
    df[(df['category'] == 'Boltz') & (df['subcategory'] == 'raw')]['clashscore'].mean()
]

# After averages - horizontal bars with labels inside
for i, proto in enumerate(relax_protocols):
    after_avgs = [
        df[(df['category'] == 'Experimental') & (df['subcategory'] == proto)]['clashscore'].mean(),
        get_af_clashscores(['ranked_0'], proto, ['ranked_0']),
        get_af_clashscores(['ranked_1', 'ranked_2', 'ranked_3', 'ranked_4'], proto,
                           ['ranked_1', 'ranked_2', 'ranked_3', 'ranked_4']),
        df[(df['category'] == 'Boltz') & (df['subcategory'] == proto)]['clashscore'].mean()
    ]
    bars = ax.barh(y + (i-2.5)*height, after_avgs, height, color=colors[proto])
    # Add protocol labels inside bars
    for bar, val in zip(bars, after_avgs):
        if not np.isnan(val) and val > 3:
            ax.text(bar.get_width() - 0.5, bar.get_y() + bar.get_height()/2,
                    short_names[proto], va='center', ha='right', fontsize=6, color='white', fontweight='bold')

# Overlay initial values as wide semi-transparent bars
bar_group_height = height * 6
ax.barh(y, before_avgs, bar_group_height, alpha=0.3, color='black', zorder=10)
# Add "Initial" labels
for yi, val in zip(y, before_avgs):
    if not np.isnan(val):
        ax.text(val + 0.5, yi, f'Initial: {val:.1f}', va='center', ha='left', fontsize=8, alpha=0.7)

ax.set_xlabel('Clash Score', labelpad=2)
ax.set_title('Average Clash Score: Before vs After Relaxation')
ax.set_yticks(y)
ax.set_yticklabels(y_labels)
ax.set_xlim(0, 45)

plt.tight_layout()
plt.savefig('figures/clashscore_averaged.png', dpi=150, bbox_inches='tight')
plt.close()
print('Saved figures/clashscore_averaged.png')

# ============ 3. Initial vs Final scatter (4 panels: Exp, AF ranked_0, AF ranked_1-4, Boltz) ============
fig, axes = plt.subplots(1, 4, figsize=(14, 3.5), sharex=True, sharey=True)
fig.subplots_adjust(wspace=0.05)

panel_configs = [
    ('Experimental', None, '#e74c3c'),
    ('AF ranked_0', ['ranked_0'], '#1a1a2e'),
    ('AF ranked_1-4', ['ranked_1', 'ranked_2', 'ranked_3', 'ranked_4'], '#3498db'),
    ('Boltz', None, '#2ecc71')
]
fixed_max = 55

for idx, (panel_name, model_subset, color) in enumerate(panel_configs):
    ax = axes[idx]

    if panel_name.startswith('AF'):
        cat_df = df[df['category'] == 'AlphaFold']
        before_sub = 'raw'
    elif panel_name == 'Experimental':
        cat_df = df[df['category'] == 'Experimental']
        before_sub = 'original'
    else:
        cat_df = df[df['category'] == 'Boltz']
        before_sub = 'raw'

    initial_means = []
    final_means = []
    final_stds = []

    for p in proteins:
        if model_subset:
            before_vals = cat_df[(cat_df['protein'] == p) &
                                 (cat_df['subcategory'] == before_sub) &
                                 (cat_df['model'].isin(model_subset))]['clashscore']
        else:
            before_vals = cat_df[(cat_df['protein'] == p) & (cat_df['subcategory'] == before_sub)]['clashscore']

        if len(before_vals) == 0:
            continue
        before_val = before_vals.mean()

        after_vals = []
        for proto in relax_protocols:
            if model_subset:
                relaxed = cat_df[(cat_df['protein'] == p) & (cat_df['subcategory'] == proto)]
                for base_model in model_subset:
                    matched = relaxed[relaxed['model'].str.startswith(base_model + '_')]
                    if len(matched) > 0:
                        after_vals.append(matched['clashscore'].mean())
            else:
                v = cat_df[(cat_df['protein'] == p) & (cat_df['subcategory'] == proto)]['clashscore'].mean()
                if not np.isnan(v):
                    after_vals.append(v)

        if not np.isnan(before_val) and len(after_vals) > 0:
            initial_means.append(before_val)
            final_means.append(np.mean(after_vals))
            final_stds.append(np.std(after_vals))

    ax.errorbar(initial_means, final_means, yerr=final_stds,
                fmt='o', color=color, alpha=0.7, capsize=3,
                markersize=7, elinewidth=1.5)

    ax.plot([0, fixed_max], [0, fixed_max], 'k--', alpha=0.5, lw=1)

    if idx == 1 or idx == 2:
        ax.set_xlabel('Initial Clash Score', labelpad=2)
    if idx == 0:
        ax.set_ylabel('Final Clash Score', labelpad=2)
    ax.set_title(panel_name, fontsize=10)
    ax.set_xlim(0, fixed_max)
    ax.set_ylim(0, fixed_max)
    ax.set_aspect('equal', adjustable='box')

plt.savefig('figures/clashscore_initial_vs_final.png', dpi=150, bbox_inches='tight')
plt.close()
print('Saved figures/clashscore_initial_vs_final.png')

print('Done!')
