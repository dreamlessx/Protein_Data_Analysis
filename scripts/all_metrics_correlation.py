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

# Load PoseBusters pass/fail from per-protein directories
pb_pf_dfs = []
for protein_dir in sorted(PROTEINS_DIR.iterdir()):
    if not protein_dir.is_dir():
        continue
    pb_file = protein_dir / "analysis" / "posebusters_results.csv"
    if pb_file.exists():
        df = pd.read_csv(pb_file)
        df['protein'] = protein_dir.name
        pb_pf_dfs.append(df)
pb_pf_df = pd.concat(pb_pf_dfs, ignore_index=True) if pb_pf_dfs else pd.DataFrame()

# MolProbity metrics - ALL 36 columns
MP_METRICS = [
    # Core metrics
    ('clashscore', 'Clashscore', True),
    ('molprobity_score', 'MolProbity Score', True),
    ('clash_count', 'Clash Count (raw)', True),
    ('atom_count', 'Atom Count', False),
    # Ramachandran
    ('rama_favored', 'Rama Favored (count)', False),
    ('rama_allowed', 'Rama Allowed (count)', False),
    ('rama_outliers', 'Rama Outliers (count)', True),
    ('rama_total', 'Rama Total', False),
    ('rama_favored_pct', 'Rama Favored %', False),
    ('rama_outliers_pct', 'Rama Outliers %', True),
    # Rotamers
    ('rota_favored', 'Rotamer Favored (count)', False),
    ('rota_outliers', 'Rotamer Outliers (count)', True),
    ('rota_total', 'Rotamer Total', False),
    ('rota_favored_pct', 'Rotamer Favored %', False),
    ('rota_outliers_pct', 'Rotamer Outliers %', True),
    # C-beta deviations
    ('cbeta_deviations', 'C-beta Deviations', True),
    ('cbeta_mean', 'C-beta Mean Dev', True),
    ('cbeta_max', 'C-beta Max Dev', True),
    ('cbeta_std', 'C-beta Std Dev', True),
    ('cbeta_rmsz', 'C-beta RMSZ', True),
    ('cbeta_n', 'C-beta N', False),
    ('cbeta_outliers', 'C-beta Outliers', True),
    # Omega angles
    ('omega_cis_proline', 'Omega Cis-Pro', False),
    ('omega_cis_general', 'Omega Cis-General', True),
    ('omega_twisted', 'Omega Twisted', True),
    ('omega_mean', 'Omega Mean', False),
    ('omega_std', 'Omega Std', True),
    ('omega_min', 'Omega Min', False),
    ('omega_max', 'Omega Max', False),
    ('omega_trans', 'Omega Trans', False),
    ('omega_cis', 'Omega Cis', False),
    # Bonds and angles
    ('bond_rmsz', 'Bond Length RMSZ', True),
    ('bond_outliers', 'Bond Outliers', True),
    ('bond_n', 'Bond N', False),
    ('angle_rmsz', 'Bond Angle RMSZ', True),
    ('angle_outliers', 'Angle Outliers', True),
    ('angle_n', 'Angle N', False),
]

# PoseBusters RAW continuous metrics - ALL columns
PB_RAW_METRICS = [
    ('raw_n_atoms', 'Atom Count', False),
    ('raw_n_residue_types', 'Residue Types', False),
    ('raw_n_valid_residue_types', 'Valid Residue Types', False),
    ('raw_n_backbone_breaks', 'Backbone Breaks', True),
    ('raw_max_break_distance', 'Max Break Distance (A)', True),
    ('raw_n_peptide_bonds', 'Peptide Bonds', False),
    ('raw_n_bond_outliers', 'Bond Outliers', True),
    ('raw_mean_bond_length', 'Mean Bond Length (A)', False),
    ('raw_std_bond_length', 'Std Bond Length', True),
    ('raw_n_backbone_angles', 'Backbone Angles', False),
    ('raw_n_angle_outliers', 'Angle Outliers', True),
    ('raw_mean_backbone_angle', 'Mean N-CA-C Angle', False),
    ('raw_std_backbone_angle', 'Std N-CA-C Angle', True),
    ('raw_n_clashes', 'Steric Clashes', True),
    ('raw_worst_clash', 'Worst Clash (A)', True),
    ('raw_n_aromatic_rings', 'Aromatic Rings', False),
    ('raw_n_nonplanar_rings', 'Non-planar Rings', True),
    ('raw_max_ring_deviation', 'Max Ring Deviation', True),
    ('raw_n_omega_angles', 'Omega Angles', False),
    ('raw_n_cis', 'Cis Peptides', False),
    ('raw_n_trans', 'Trans Peptides', False),
    ('raw_n_twisted', 'Twisted Peptides', True),
    ('raw_mean_abs_omega', 'Mean |Omega|', False),
    ('raw_n_chiral_centers', 'Chiral Centers', False),
    ('raw_n_d_amino_acids', 'D-Amino Acids', True),
    ('raw_n_residues', 'Residues', False),
    ('raw_n_incomplete_residues', 'Incomplete Residues', True),
    ('raw_n_missing_backbone_atoms', 'Missing BB Atoms', True),
    ('raw_rosetta_score', 'Rosetta Energy', True),
]

# PoseBusters PASS/FAIL metrics (as 0/1)
PB_PASSFAIL_METRICS = [
    ('structure_loaded', 'Structure Loaded', False),  # higher = better (pass)
    ('valid_residues', 'Valid Residues', False),
    ('backbone_connected', 'Backbone Connected', False),
    ('bond_lengths', 'Bond Lengths OK', False),
    ('bond_angles', 'Bond Angles OK', False),
    ('steric_clashes', 'No Steric Clashes', False),
    ('aromatic_flatness', 'Aromatic Flatness', False),
    ('peptide_planarity', 'Peptide Planarity', False),
    ('chirality', 'Chirality OK', False),
    ('complete_residues', 'Complete Residues', False),
    ('internal_energy', 'Internal Energy OK', False),
    ('all_pass', 'All Tests Pass', False),
    ('n_pass', 'Tests Passed (count)', False),
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

    # PoseBusters RAW metrics
    print("\n--- PoseBusters RAW Continuous Metrics ---")
    for col, name, lower_better in PB_RAW_METRICS:
        if col not in pb_df.columns:
            print(f"  {name}: column not found, skipping")
            continue

        data = collect_metric_data(pb_df, col, is_molprobity=False)
        if len(data) == 0:
            print(f"  {name}: no data")
            continue

        filename = f'corr_pbraw_{col.replace("raw_", "")}.png'
        r, p = plot_correlation(data, name, lower_better, filename)
        print(f"  {name}: r={r:.3f}, p={p:.6f} -> {filename}")
        results.append({'metric': name, 'source': 'PB-Raw', 'r': r, 'p': p, 'lower_better': lower_better})

    # PoseBusters PASS/FAIL metrics
    print("\n--- PoseBusters Pass/Fail Metrics ---")
    for col, name, lower_better in PB_PASSFAIL_METRICS:
        if col not in pb_pf_df.columns:
            print(f"  {name}: column not found, skipping")
            continue

        data = collect_metric_data(pb_pf_df, col, is_molprobity=False)
        if len(data) == 0:
            print(f"  {name}: no data")
            continue

        filename = f'corr_pbpf_{col}.png'
        r, p = plot_correlation(data, name, lower_better, filename)
        print(f"  {name}: r={r:.3f}, p={p:.6f} -> {filename}")
        results.append({'metric': name, 'source': 'PB-PassFail', 'r': r, 'p': p, 'lower_better': lower_better})

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
