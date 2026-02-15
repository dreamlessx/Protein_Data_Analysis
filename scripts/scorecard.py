#!/usr/bin/env python3
"""Generate comprehensive scorecard for Wednesday presentation."""

import pandas as pd
import numpy as np
from pathlib import Path

RESULTS_DIR = Path(__file__).parent.parent / "validation_results"
PROTEINS_DIR = Path(__file__).parent.parent / "proteins"

# Load MolProbity
mp_df = pd.read_csv(RESULTS_DIR / "molprobity_full.csv")

# Load PoseBusters raw from per-protein directories
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

# Load PoseBusters pass/fail
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

# All relaxation protocols
PROTOCOLS = ['relaxed_normal_beta', 'relaxed_normal_ref15', 'relaxed_cartesian_beta',
             'relaxed_cartesian_ref15', 'relaxed_dualspace_beta', 'relaxed_dualspace_ref15']

SOURCE_CONFIGS = [
    # (name, baseline_subcat, model_filter)
    ('Experimental', 'original', None),
    ('AF ranked_0 (AMBER)', 'raw', 'ranked_0'),
    ('AF ranked_1-4 (unrelaxed)', 'raw', 'not_ranked_0'),
    ('Boltz', 'raw', None),
]


def get_baseline_relaxed(df, category, baseline_subcat, model_filter, is_molprobity=True):
    """Get baseline and all relaxed data for a source."""
    cat_df = df[df['category'] == category]

    # Baseline
    baseline = cat_df[cat_df['subcategory'] == baseline_subcat]
    if model_filter == 'ranked_0':
        if is_molprobity:
            baseline = baseline[baseline['model'] == 'ranked_0']
        else:
            baseline = baseline[baseline['model'] == 'model0']
    elif model_filter == 'not_ranked_0':
        if is_molprobity:
            baseline = baseline[~baseline['model'].isin(['ranked_0'])]
        else:
            baseline = baseline[~baseline['model'].isin(['model0'])]

    # Relaxed (all protocols)
    relaxed_all = {}
    for protocol in PROTOCOLS:
        rel = cat_df[cat_df['subcategory'] == protocol]
        if model_filter == 'ranked_0':
            if is_molprobity:
                rel = rel[rel['model'].str.startswith('ranked_0')]
            else:
                rel = rel[rel['model'].str.startswith('ranked_0')]
        elif model_filter == 'not_ranked_0':
            if is_molprobity:
                rel = rel[~rel['model'].str.startswith('ranked_0')]
            else:
                rel = rel[~rel['model'].str.startswith('ranked_0')]
        relaxed_all[protocol] = rel

    return baseline, relaxed_all


def compute_metric_stats(baseline, relaxed_all, metric_col):
    """Compute before/after stats with best/worst protocols."""
    if metric_col not in baseline.columns:
        return None

    before_mean = baseline[metric_col].mean()
    before_min = baseline[metric_col].min()
    before_max = baseline[metric_col].max()

    # Aggregate across all protocols
    all_relaxed = pd.concat(relaxed_all.values())
    if len(all_relaxed) == 0:
        return None

    after_mean = all_relaxed[metric_col].mean()
    after_min = all_relaxed[metric_col].min()
    after_max = all_relaxed[metric_col].max()

    # Best and worst protocol (by mean)
    protocol_means = {}
    for prot, rel_df in relaxed_all.items():
        if len(rel_df) > 0:
            protocol_means[prot] = rel_df[metric_col].mean()

    if protocol_means:
        best_protocol = min(protocol_means, key=protocol_means.get)  # Assuming lower is better
        worst_protocol = max(protocol_means, key=protocol_means.get)
        best_val = protocol_means[best_protocol]
        worst_val = protocol_means[worst_protocol]
    else:
        best_protocol = worst_protocol = None
        best_val = worst_val = np.nan

    delta = after_mean - before_mean

    return {
        'before_mean': before_mean,
        'before_min': before_min,
        'before_max': before_max,
        'after_mean': after_mean,
        'after_min': after_min,
        'after_max': after_max,
        'delta': delta,
        'best_protocol': best_protocol.replace('relaxed_', '') if best_protocol else None,
        'best_val': best_val,
        'worst_protocol': worst_protocol.replace('relaxed_', '') if worst_protocol else None,
        'worst_val': worst_val,
    }


def compute_passfail_rate(baseline, relaxed_all, metric_col):
    """Compute pass rate before/after."""
    if metric_col not in baseline.columns:
        return None

    before_rate = baseline[metric_col].mean() * 100  # Convert to %

    all_relaxed = pd.concat(relaxed_all.values())
    if len(all_relaxed) == 0:
        return None

    after_rate = all_relaxed[metric_col].mean() * 100

    delta = after_rate - before_rate

    return {
        'before': before_rate,
        'after': after_rate,
        'delta': delta,
    }


def main():
    print("=" * 100)
    print("COMPREHENSIVE SCORECARD FOR WEDNESDAY PRESENTATION")
    print("=" * 100)

    # MolProbity metrics
    MP_METRICS = [
        ('clashscore', 'Clashscore', '↓'),
        ('molprobity_score', 'MolProbity Score', '↓'),
        ('rama_outliers_pct', 'Rama Outliers %', '↓'),
        ('rama_favored_pct', 'Rama Favored %', '↑'),
        ('rota_outliers_pct', 'Rotamer Outliers %', '↓'),
        ('cbeta_deviations', 'C-beta Deviations', '↓'),
        ('bond_rmsz', 'Bond RMSZ', '↓'),
        ('angle_rmsz', 'Angle RMSZ', '↓'),
        ('clash_count', 'Clash Count (raw)', '↓'),
        ('omega_twisted', 'Omega Twisted', '↓'),
    ]

    # PoseBusters binary tests
    PB_BINARY = [
        ('bond_lengths', 'Bond Lengths', '↑'),
        ('bond_angles', 'Bond Angles', '↑'),
        ('steric_clashes', 'Steric Clashes', '↑'),
        ('aromatic_flatness', 'Aromatic Flatness', '↑'),
        ('peptide_planarity', 'Peptide Planarity', '↑'),
        ('internal_energy', 'Internal Energy', '↑'),
        ('backbone_connected', 'Backbone Connected', '↑'),
        ('chirality', 'Chirality', '↑'),
        ('complete_residues', 'Complete Residues', '↑'),
        ('all_pass', 'ALL TESTS PASS', '↑'),
    ]

    # PoseBusters continuous
    PB_CONT = [
        ('raw_n_clashes', 'PB Steric Clashes', '↓'),
        ('raw_worst_clash', 'PB Worst Clash (Å)', '↓'),
        ('raw_n_bond_outliers', 'PB Bond Outliers', '↓'),
        ('raw_n_angle_outliers', 'PB Angle Outliers', '↓'),
        ('raw_n_twisted', 'PB Twisted Omega', '↓'),
        ('raw_rosetta_score', 'PB Rosetta Energy', '↓'),
        ('raw_mean_bond_length', 'PB Mean Bond (Å)', '—'),
        ('raw_mean_backbone_angle', 'PB Mean N-CA-C (°)', '—'),
    ]

    results = []

    # Process MolProbity
    print("\n### MolProbity Metrics")
    print("-" * 100)

    for col, name, direction in MP_METRICS:
        row = {'metric': name, 'direction': direction, 'type': 'MolProbity'}

        for source_name, baseline_subcat, model_filter in SOURCE_CONFIGS:
            category = source_name.split()[0] if source_name.startswith('AF') else source_name
            if 'AF' in source_name:
                category = 'AlphaFold'

            baseline, relaxed_all = get_baseline_relaxed(mp_df, category, baseline_subcat, model_filter, True)
            stats = compute_metric_stats(baseline, relaxed_all, col)

            if stats:
                row[f'{source_name}_before'] = stats['before_mean']
                row[f'{source_name}_after'] = stats['after_mean']
                row[f'{source_name}_delta'] = stats['delta']
                row[f'{source_name}_range'] = f"[{stats['after_min']:.1f}, {stats['after_max']:.1f}]"
                row[f'{source_name}_best'] = f"{stats['best_protocol']}: {stats['best_val']:.2f}" if stats['best_protocol'] else "N/A"
                row[f'{source_name}_worst'] = f"{stats['worst_protocol']}: {stats['worst_val']:.2f}" if stats['worst_protocol'] else "N/A"

                # Flag if worse
                if direction == '↓' and stats['delta'] > 0.1:
                    row[f'{source_name}_flag'] = '⚠️ WORSE'
                elif direction == '↑' and stats['delta'] < -0.1:
                    row[f'{source_name}_flag'] = '⚠️ WORSE'
                else:
                    row[f'{source_name}_flag'] = ''

        results.append(row)
        print(f"{name}: processed")

    # Process PoseBusters binary
    print("\n### PoseBusters Binary Tests (% passing)")
    print("-" * 100)

    for col, name, direction in PB_BINARY:
        row = {'metric': name, 'direction': direction, 'type': 'PB Binary'}

        for source_name, baseline_subcat, model_filter in SOURCE_CONFIGS:
            category = source_name.split()[0] if source_name.startswith('AF') else source_name
            if 'AF' in source_name:
                category = 'AlphaFold'

            baseline, relaxed_all = get_baseline_relaxed(pb_pf_df, category, baseline_subcat, model_filter, False)
            stats = compute_passfail_rate(baseline, relaxed_all, col)

            if stats:
                row[f'{source_name}_before'] = stats['before']
                row[f'{source_name}_after'] = stats['after']
                row[f'{source_name}_delta'] = stats['delta']

                # Flag if worse
                if stats['delta'] < -5:
                    row[f'{source_name}_flag'] = '⚠️ WORSE'
                else:
                    row[f'{source_name}_flag'] = ''

        results.append(row)
        print(f"{name}: processed")

    # Process PoseBusters continuous
    print("\n### PoseBusters Continuous Metrics")
    print("-" * 100)

    for col, name, direction in PB_CONT:
        row = {'metric': name, 'direction': direction, 'type': 'PB Continuous'}

        for source_name, baseline_subcat, model_filter in SOURCE_CONFIGS:
            category = source_name.split()[0] if source_name.startswith('AF') else source_name
            if 'AF' in source_name:
                category = 'AlphaFold'

            baseline, relaxed_all = get_baseline_relaxed(pb_df, category, baseline_subcat, model_filter, False)
            stats = compute_metric_stats(baseline, relaxed_all, col)

            if stats:
                row[f'{source_name}_before'] = stats['before_mean']
                row[f'{source_name}_after'] = stats['after_mean']
                row[f'{source_name}_delta'] = stats['delta']
                row[f'{source_name}_range'] = f"[{stats['after_min']:.1f}, {stats['after_max']:.1f}]"

                # Flag if worse
                if direction == '↓' and stats['delta'] > 0.5:
                    row[f'{source_name}_flag'] = '⚠️ WORSE'
                else:
                    row[f'{source_name}_flag'] = ''

        results.append(row)
        print(f"{name}: processed")

    # Save results
    df_results = pd.DataFrame(results)
    df_results.to_csv(RESULTS_DIR / "scorecard.csv", index=False)
    print(f"\nSaved scorecard.csv")

    # Generate markdown table
    print("\n" + "=" * 100)
    print("MARKDOWN TABLE FOR MEETING_PREP.md")
    print("=" * 100)

    sources = ['Experimental', 'AF ranked_0 (AMBER)', 'AF ranked_1-4 (unrelaxed)', 'Boltz']

    print("\n### MolProbity Scorecard\n")
    print("| Metric | Dir | " + " | ".join(sources) + " |")
    print("|--------|-----|" + "|".join(["---" for _ in sources]) + "|")

    for _, row in df_results[df_results['type'] == 'MolProbity'].iterrows():
        cells = [row['metric'], row['direction']]
        for src in sources:
            before = row.get(f'{src}_before', np.nan)
            after = row.get(f'{src}_after', np.nan)
            delta = row.get(f'{src}_delta', np.nan)
            flag = row.get(f'{src}_flag', '')

            if pd.notna(before) and pd.notna(after):
                if abs(before) > 100 or abs(after) > 100:
                    cell = f"{before:.0f}→{after:.0f} ({delta:+.0f})"
                else:
                    cell = f"{before:.1f}→{after:.1f} ({delta:+.1f})"
                if flag:
                    cell += f" {flag}"
            else:
                cell = "—"
            cells.append(cell)
        print("| " + " | ".join(cells) + " |")

    print("\n### PoseBusters Binary Tests (% passing)\n")
    print("| Test | Dir | " + " | ".join(sources) + " |")
    print("|------|-----|" + "|".join(["---" for _ in sources]) + "|")

    for _, row in df_results[df_results['type'] == 'PB Binary'].iterrows():
        cells = [row['metric'], row['direction']]
        for src in sources:
            before = row.get(f'{src}_before', np.nan)
            after = row.get(f'{src}_after', np.nan)
            delta = row.get(f'{src}_delta', np.nan)
            flag = row.get(f'{src}_flag', '')

            if pd.notna(before) and pd.notna(after):
                cell = f"{before:.0f}%→{after:.0f}% ({delta:+.0f})"
                if flag:
                    cell += f" {flag}"
            else:
                cell = "—"
            cells.append(cell)
        print("| " + " | ".join(cells) + " |")

    print("\n### PoseBusters Continuous Metrics\n")
    print("| Metric | Dir | " + " | ".join(sources) + " |")
    print("|--------|-----|" + "|".join(["---" for _ in sources]) + "|")

    for _, row in df_results[df_results['type'] == 'PB Continuous'].iterrows():
        cells = [row['metric'], row['direction']]
        for src in sources:
            before = row.get(f'{src}_before', np.nan)
            after = row.get(f'{src}_after', np.nan)
            delta = row.get(f'{src}_delta', np.nan)
            flag = row.get(f'{src}_flag', '')

            if pd.notna(before) and pd.notna(after):
                if abs(before) > 1000 or abs(after) > 1000:
                    cell = f"{before:.0f}→{after:.0f}"
                elif abs(before) > 10 or abs(after) > 10:
                    cell = f"{before:.1f}→{after:.1f} ({delta:+.1f})"
                else:
                    cell = f"{before:.2f}→{after:.2f} ({delta:+.2f})"
                if flag:
                    cell += f" {flag}"
            else:
                cell = "—"
            cells.append(cell)
        print("| " + " | ".join(cells) + " |")


if __name__ == "__main__":
    main()
