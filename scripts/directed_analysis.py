#!/usr/bin/env python3
"""
directed_analysis.py
Generate pass/fail tables for PoseBusters results.
"""

import pandas as pd
from pathlib import Path

PROJECT_DIR = Path(__file__).parent.parent
INPUT_FILE = PROJECT_DIR / "analysis_output" / "protein_busters_results.csv"
OUTPUT_DIR = PROJECT_DIR / "directed_analysis"

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# All test columns from protein_busters.py
TEST_COLS = [
    'bond_pass', 'omega_pass', 'chiral_pass', 'backbone_pass',
    'phipsi_pass', 'disulfide_pass', 'all_tests_pass'
]

PROTOCOLS = ['original', 'raw', 'relaxed_cartesian_beta', 'relaxed_cartesian_ref15',
             'relaxed_dualspace_beta', 'relaxed_dualspace_ref15',
             'relaxed_normal_beta', 'relaxed_normal_ref15']


def generate_table(df: pd.DataFrame, output_name: str):
    """Generate pass percentage table grouped by protein and protocol."""
    results = []

    for protein in sorted(df['protein'].unique()):
        prot_df = df[df['protein'] == protein]

        for protocol in PROTOCOLS:
            sub_df = prot_df[prot_df['subcategory'] == protocol]
            if len(sub_df) == 0:
                continue

            n = len(sub_df)
            row = {'protein': protein, 'protocol': protocol, 'n': n}

            for col in TEST_COLS:
                if col not in sub_df.columns:
                    continue
                pct = 100 * sub_df[col].sum() / n
                col_name = col.replace('_pass', '_pass_pct').replace('all_tests_pass_pct', 'all_pass_pct')
                row[col_name] = round(pct, 1)

            # Real AllPass including phipsi
            all_cols = ['bond_pass', 'omega_pass', 'chiral_pass', 'backbone_pass', 'phipsi_pass', 'disulfide_pass']
            real_all = sub_df[all_cols].all(axis=1).sum()
            row['real_all_pass_pct'] = round(100 * real_all / n, 1)

            results.append(row)

    result_df = pd.DataFrame(results)
    output_path = OUTPUT_DIR / output_name
    result_df.to_csv(output_path, index=False)
    print(f"Saved: {output_path}")
    print(result_df.head(20).to_string(index=False))
    print()

    return result_df


def generate_combined_summary(df: pd.DataFrame):
    """Generate combined summary table for all categories and protocols."""
    results = []

    # Experimental original first
    sub_df = df[(df['category'] == 'Experimental') & (df['subcategory'] == 'original')]
    n = len(sub_df)
    all_cols = ['bond_pass', 'omega_pass', 'chiral_pass', 'backbone_pass', 'phipsi_pass', 'disulfide_pass']
    real_all = sub_df[all_cols].all(axis=1).sum()
    results.append({
        'Category': 'Experimental',
        'Protocol': 'original',
        'N': n,
        'Bond': f"{100 * sub_df['bond_pass'].sum() / n:.1f}%",
        'Omega': f"{100 * sub_df['omega_pass'].sum() / n:.1f}%",
        'Chiral': f"{100 * sub_df['chiral_pass'].sum() / n:.1f}%",
        'Backbone': f"{100 * sub_df['backbone_pass'].sum() / n:.1f}%",
        'PhiPsi': f"{100 * sub_df['phipsi_pass'].sum() / n:.1f}%",
        'Disulfide': f"{100 * sub_df['disulfide_pass'].sum() / n:.1f}%",
        'AllPass': f"{100 * sub_df['all_tests_pass'].sum() / n:.1f}%",
        'RealAllPass': f"{100 * real_all / n:.1f}%",
    })

    # AlphaFold and Boltz protocols
    protocols = ['raw', 'relaxed_cartesian_beta', 'relaxed_cartesian_ref15',
                 'relaxed_dualspace_beta', 'relaxed_dualspace_ref15',
                 'relaxed_normal_beta', 'relaxed_normal_ref15']

    for protocol in protocols:
        for cat in ['AlphaFold', 'Boltz']:
            sub_df = df[(df['category'] == cat) & (df['subcategory'] == protocol)]
            if len(sub_df) == 0:
                continue

            n = len(sub_df)
            all_cols = ['bond_pass', 'omega_pass', 'chiral_pass', 'backbone_pass', 'phipsi_pass', 'disulfide_pass']
            real_all = sub_df[all_cols].all(axis=1).sum()

            row = {
                'Category': cat,
                'Protocol': protocol.replace('relaxed_', ''),
                'N': n,
                'Bond': f"{100 * sub_df['bond_pass'].sum() / n:.1f}%",
                'Omega': f"{100 * sub_df['omega_pass'].sum() / n:.1f}%",
                'Chiral': f"{100 * sub_df['chiral_pass'].sum() / n:.1f}%",
                'Backbone': f"{100 * sub_df['backbone_pass'].sum() / n:.1f}%",
                'PhiPsi': f"{100 * sub_df['phipsi_pass'].sum() / n:.1f}%",
                'Disulfide': f"{100 * sub_df['disulfide_pass'].sum() / n:.1f}%",
                'AllPass': f"{100 * sub_df['all_tests_pass'].sum() / n:.1f}%",
                'RealAllPass': f"{100 * real_all / n:.1f}%",
            }
            results.append(row)

    result_df = pd.DataFrame(results)
    output_path = OUTPUT_DIR / "table4_summary.csv"
    result_df.to_csv(output_path, index=False)
    print(f"Saved: {output_path}")
    print(result_df.to_string(index=False))
    print()

    return result_df


def main():
    df = pd.read_csv(INPUT_FILE)

    print("=" * 80)
    print("PROTEIN BUSTERS DIRECTED ANALYSIS")
    print("=" * 80)

    # Table 1: Experimental (original + all relaxed)
    print("\n--- Table 1: Experimental (All) ---")
    exp_df = df[df['category'] == 'Experimental']
    generate_table(exp_df, "table1_experimental.csv")

    # Table 2: Boltz (raw + all relaxed)
    print("\n--- Table 2: Boltz (All) ---")
    boltz_df = df[df['category'] == 'Boltz']
    generate_table(boltz_df, "table2_boltz.csv")

    # Table 3: AlphaFold (raw + all relaxed)
    print("\n--- Table 3: AlphaFold (All) ---")
    af_df = df[df['category'] == 'AlphaFold']
    generate_table(af_df, "table3_alphafold.csv")

    # Table 4: Combined summary (all categories and protocols)
    print("\n--- Table 4: Combined Summary ---")
    generate_combined_summary(df)

    print("=" * 80)
    print("Done!")


if __name__ == "__main__":
    main()
