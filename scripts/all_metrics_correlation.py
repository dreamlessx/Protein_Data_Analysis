#!/usr/bin/env python3
"""Generate correlation plots for all MolProbity and PoseBusters continuous metrics."""

import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
from scipy import stats

RESULTS_DIR = Path(__file__).parent.parent / "validation_results"
FIGURES_DIR = Path(__file__).parent.parent / "figures"

# Load data
mp_df = pd.read_csv(RESULTS_DIR / "molprobity_full.csv")

# Load PoseBusters from per-protein directories (has protein column)
PROTEINS_DIR = Path(__file__).parent.parent / "proteins"
pb_dfs = []
for protein_dir in sorted(PROTEINS_DIR.iterdir()):
    if not protein_dir.is_dir():
        continue
    pb_file = protein_dir / "analysis" / "posebusters_raw.csv"
    if pb_file.exists():
        df = pd.read_csv(pb_file)
        df['protein'] = protein_dir.name
        pb_dfs.append(df)
pb_df = pd.concat(pb_dfs, ignore_index=True) if pb_dfs else pd.DataFrame()

# MolProbity metrics (lower = better for all)
MP_METRICS = [
    ('clashscore', 'Clashscore', True),  # (column, display_name, lower_is_better)
    ('molprobity_score', 'MolProbity Score', True),
    ('rama_outliers', 'Ramachandran Outliers (%)', True),
    ('rota_outliers', 'Rotamer Outliers (%)', True),
    ('bond_rmsz', 'Bond Length RMSZ', True),
    ('angle_rmsz', 'Bond Angle RMSZ', True),
]

# PoseBusters continuous metrics
PB_METRICS = [
    ('raw_n_clashes', 'Steric Clashes (count)', True),
    ('raw_worst_clash', 'Worst Clash Overlap (A)', True),
    ('raw_n_bond_outliers', 'Bond Length Outliers', True),
    ('raw_mean_bond_length', 'Mean Peptide Bond (A)', False),  # ~1.33 is ideal
    ('raw_n_angle_outliers', 'Angle Outliers', True),
    ('raw_mean_backbone_angle', 'Mean N-CA-C Angle', False),  # ~111 is ideal
    ('raw_n_twisted', 'Twisted Omega Count', True),
    ('raw_rosetta_score', 'Rosetta Energy (REU)', True),
]


def collect_metric_data(df, metric_col, category_col='category', subcat_col='subcategory',
                        model_col='model', protein_col='protein', is_molprobity=True):
    """Collect paired data for a metric: initial value vs change after relaxation."""
    data = []

    for category in ['Experimental', 'AlphaFold', 'Boltz']:
        cat_df = df[df[category_col] == category]

        if category == 'Experimental':
            baseline = cat_df[cat_df[subcat_col] == 'original']
            relaxed = cat_df[cat_df[subcat_col] == 'relaxed_normal_beta']

            for protein in baseline[protein_col].unique():
                base_val = baseline[baseline[protein_col] == protein][metric_col].mean()
                rel_data = relaxed[relaxed[protein_col] == protein][metric_col]
                if len(rel_data) == 0:
                    continue
                rel_mean = rel_data.mean()
                change = rel_mean - base_val
                change_min = rel_data.min() - base_val
                change_max = rel_data.max() - base_val

                if pd.notna(base_val) and pd.notna(rel_mean):
                    data.append({
                        'protein': protein,
                        'category': category,
                        'af_type': None,
                        'initial': base_val,
                        'final': rel_mean,
                        'change': change,
                        'change_lo': change - change_min,
                        'change_hi': change_max - change,
                    })

        elif category == 'AlphaFold':
            baseline = cat_df[cat_df[subcat_col] == 'raw']
            relaxed = cat_df[cat_df[subcat_col] == 'relaxed_normal_beta']

            for protein in baseline[protein_col].unique():
                # ranked_0 (AMBER-relaxed)
                if is_molprobity:
                    base_r0 = baseline[(baseline[protein_col] == protein) & (baseline[model_col] == 'ranked_0')]
                    rel_r0 = relaxed[(relaxed[protein_col] == protein) & (relaxed[model_col].str.startswith('ranked_0'))]
                else:
                    base_r0 = baseline[(baseline[protein_col] == protein) & (baseline[model_col] == 'model0')]
                    rel_r0 = relaxed[(relaxed[protein_col] == protein) & (relaxed[model_col].str.startswith('ranked_0'))]

                if len(base_r0) and len(rel_r0):
                    base_val = base_r0[metric_col].mean()
                    rel_mean = rel_r0[metric_col].mean()
                    change = rel_mean - base_val
                    change_min = rel_r0[metric_col].min() - base_val
                    change_max = rel_r0[metric_col].max() - base_val

                    if pd.notna(base_val) and pd.notna(rel_mean):
                        data.append({
                            'protein': protein,
                            'category': category,
                            'af_type': 'ranked_0 (AMBER)',
                            'initial': base_val,
                            'final': rel_mean,
                            'change': change,
                            'change_lo': change - change_min,
                            'change_hi': change_max - change,
                        })

                # ranked_1-4 (unrelaxed)
                if is_molprobity:
                    base_rest = baseline[(baseline[protein_col] == protein) & (~baseline[model_col].isin(['ranked_0']))]
                    rel_rest = relaxed[(relaxed[protein_col] == protein) & (~relaxed[model_col].str.startswith('ranked_0'))]
                else:
                    base_rest = baseline[(baseline[protein_col] == protein) & (~baseline[model_col].isin(['model0']))]
                    rel_rest = relaxed[(relaxed[protein_col] == protein) & (~relaxed[model_col].str.startswith('ranked_0'))]

                if len(base_rest) and len(rel_rest):
                    base_val = base_rest[metric_col].mean()
                    rel_mean = rel_rest[metric_col].mean()
                    change = rel_mean - base_val
                    change_min = rel_rest[metric_col].min() - base_val
                    change_max = rel_rest[metric_col].max() - base_val

                    if pd.notna(base_val) and pd.notna(rel_mean):
                        data.append({
                            'protein': protein,
                            'category': category,
                            'af_type': 'ranked_1-4 (unrelaxed)',
                            'initial': base_val,
                            'final': rel_mean,
                            'change': change,
                            'change_lo': change - change_min,
                            'change_hi': change_max - change,
                        })
        else:  # Boltz
            baseline = cat_df[cat_df[subcat_col] == 'raw']
            relaxed = cat_df[cat_df[subcat_col] == 'relaxed_normal_beta']

            for protein in baseline[protein_col].unique():
                base_val = baseline[baseline[protein_col] == protein][metric_col].mean()
                rel_data = relaxed[relaxed[protein_col] == protein][metric_col]
                if len(rel_data) == 0:
                    continue
                rel_mean = rel_data.mean()
                change = rel_mean - base_val
                change_min = rel_data.min() - base_val
                change_max = rel_data.max() - base_val

                if pd.notna(base_val) and pd.notna(rel_mean):
                    data.append({
                        'protein': protein,
                        'category': category,
                        'af_type': None,
                        'initial': base_val,
                        'final': rel_mean,
                        'change': change,
                        'change_lo': change - change_min,
                        'change_hi': change_max - change,
                    })

    return pd.DataFrame(data)


def plot_correlation(data_df, metric_name, lower_is_better, filename):
    """Create 3-panel correlation plot."""
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    colors = {'Experimental': 'green', 'AlphaFold': 'blue', 'Boltz': 'red'}

    # Calculate overall correlation
    if len(data_df) > 2 and data_df['initial'].std() > 0 and data_df['change'].std() > 0:
        r_all, p_all = stats.pearsonr(data_df['initial'], data_df['change'])
    else:
        r_all, p_all = 0, 1

    for i, cat in enumerate(['Experimental', 'AlphaFold', 'Boltz']):
        ax = axes[i]
        cat_data = data_df[data_df['category'] == cat]

        if len(cat_data) == 0:
            ax.set_title(f'{cat} (no data)')
            continue

        if cat == 'AlphaFold':
            r0_data = cat_data[cat_data['af_type'] == 'ranked_0 (AMBER)']
            rest_data = cat_data[cat_data['af_type'] == 'ranked_1-4 (unrelaxed)']

            if len(r0_data):
                ax.errorbar(r0_data['initial'], r0_data['change'],
                            yerr=[r0_data['change_lo'].clip(lower=0), r0_data['change_hi'].clip(lower=0)],
                            fmt='o', c='black', alpha=0.7, markersize=8,
                            label='ranked_0 (AMBER)', capsize=3)
            if len(rest_data):
                ax.errorbar(rest_data['initial'], rest_data['change'],
                            yerr=[rest_data['change_lo'].clip(lower=0), rest_data['change_hi'].clip(lower=0)],
                            fmt='o', c='darkblue', alpha=0.7, markersize=8,
                            label='ranked_1-4 (unrelaxed)', capsize=3)
            ax.legend(loc='best', fontsize=7)
        else:
            ax.errorbar(cat_data['initial'], cat_data['change'],
                        yerr=[cat_data['change_lo'].clip(lower=0), cat_data['change_hi'].clip(lower=0)],
                        fmt='o', c=colors[cat], alpha=0.7, markersize=8, capsize=3)

        ax.axhline(0, color='black', linestyle='--', alpha=0.5)
        ax.set_xlabel(f'Initial {metric_name}')

        direction = '- = better' if lower_is_better else '+ = better'
        ax.set_ylabel(f'Change ({direction})')

        # Correlation for this category
        if len(cat_data) > 2 and cat_data['initial'].std() > 0 and cat_data['change'].std() > 0:
            r, p = stats.pearsonr(cat_data['initial'], cat_data['change'])
            ax.set_title(f'{cat} (r={r:.2f})')
        else:
            ax.set_title(f'{cat}')

        # Regression line
        if len(cat_data) > 1 and cat_data['initial'].std() > 0:
            try:
                z = np.polyfit(cat_data['initial'], cat_data['change'], 1)
                poly = np.poly1d(z)
                x_line = np.linspace(cat_data['initial'].min(), cat_data['initial'].max(), 100)
                ax.plot(x_line, poly(x_line), 'k-', alpha=0.3)
            except (np.linalg.LinAlgError, ValueError):
                pass  # Skip regression line if data is problematic

    plt.suptitle(f'{metric_name}: Initial vs Change (overall r={r_all:.2f})', fontsize=14)
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / filename, dpi=150)
    plt.close()

    return r_all, p_all


def main():
    print("=" * 70)
    print("GENERATING CORRELATION PLOTS FOR ALL METRICS")
    print("=" * 70)

    results = []

    # MolProbity metrics
    print("\n--- MolProbity Metrics ---")
    for col, name, lower_better in MP_METRICS:
        if col not in mp_df.columns:
            print(f"  {name}: column not found, skipping")
            continue

        data = collect_metric_data(mp_df, col, is_molprobity=True)
        if len(data) == 0:
            print(f"  {name}: no data")
            continue

        filename = f'corr_mp_{col}.png'
        r, p = plot_correlation(data, name, lower_better, filename)
        print(f"  {name}: r={r:.3f}, p={p:.6f} -> {filename}")
        results.append({'metric': name, 'source': 'MolProbity', 'r': r, 'p': p, 'lower_better': lower_better})

    # PoseBusters metrics
    print("\n--- PoseBusters Continuous Metrics ---")
    for col, name, lower_better in PB_METRICS:
        if col not in pb_df.columns:
            print(f"  {name}: column not found, skipping")
            continue

        data = collect_metric_data(pb_df, col, is_molprobity=False)
        if len(data) == 0:
            print(f"  {name}: no data")
            continue

        filename = f'corr_pb_{col.replace("raw_", "")}.png'
        r, p = plot_correlation(data, name, lower_better, filename)
        print(f"  {name}: r={r:.3f}, p={p:.6f} -> {filename}")
        results.append({'metric': name, 'source': 'PoseBusters', 'r': r, 'p': p, 'lower_better': lower_better})

    # Summary table
    print("\n" + "=" * 70)
    print("SUMMARY: Correlation between initial value and change")
    print("=" * 70)
    print(f"\n{'Metric':<35} {'r':>8} {'p-value':>12} {'Interpretation':<30}")
    print("-" * 90)

    for res in results:
        if res['r'] < -0.7:
            interp = "Strong convergence (regression to mean)"
        elif res['r'] < -0.3:
            interp = "Moderate convergence"
        elif res['r'] > 0.3:
            interp = "Divergence (gets more extreme)"
        else:
            interp = "No clear pattern"

        print(f"{res['metric']:<35} {res['r']:>8.3f} {res['p']:>12.6f} {interp:<30}")

    print(f"\nGenerated {len(results)} correlation plots in figures/")


if __name__ == "__main__":
    main()
