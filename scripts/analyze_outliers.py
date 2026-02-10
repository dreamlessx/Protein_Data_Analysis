#!/usr/bin/env python3
"""Analyze why certain structures degrade after relaxation."""

import pandas as pd
from pathlib import Path

RESULTS_DIR = Path(__file__).parent.parent / "validation_results"

df = pd.read_csv(RESULTS_DIR / "molprobity_full.csv")

outliers = ['2I25', '1AY7', '1AVX', '1VFB', '1BVN']

print("=" * 60)
print("OUTLIER ANALYSIS: Structures that got worse after relaxation")
print("=" * 60)

for prot in outliers:
    prot_df = df[df['protein'] == prot]
    print(f'\n{"=" * 40}')
    print(f'{prot}')
    print(f'{"=" * 40}')
    print(f'Total models: {prot_df["model"].nunique()}')

    for category in ['Experimental', 'AlphaFold', 'Boltz']:
        cat_df = prot_df[prot_df['category'] == category]
        if cat_df.empty:
            continue

        if category == 'Experimental':
            raw = cat_df[cat_df['subcategory'] == 'original']
        else:
            raw = cat_df[cat_df['subcategory'] == 'raw']

        relaxed = cat_df[cat_df['subcategory'] == 'relaxed_normal_beta']

        if raw.empty or relaxed.empty:
            continue

        print(f'\n  {category}:')
        print(f'    Clashscore: {raw["clashscore"].mean():.1f} → {relaxed["clashscore"].mean():.1f}')
        print(f'    Bond RMSZ:  {raw["bond_rmsz"].mean():.2f} → {relaxed["bond_rmsz"].mean():.2f}')
        print(f'    Angle RMSZ: {raw["angle_rmsz"].mean():.2f} → {relaxed["angle_rmsz"].mean():.2f}')
        print(f'    MolProbity: {raw["molprobity_score"].mean():.2f} → {relaxed["molprobity_score"].mean():.2f}')

# Summary table
print("\n" + "=" * 60)
print("METRIC BREAKDOWN FOR WORST CASE: 2I25 Experimental")
print("=" * 60)

prot_df = df[(df['protein'] == '2I25') & (df['category'] == 'Experimental')]
orig = prot_df[prot_df['subcategory'] == 'original']
relax = prot_df[prot_df['subcategory'] == 'relaxed_normal_beta']

metrics = ['clashscore', 'rama_outliers_pct', 'rota_outliers_pct', 'cbeta_outliers',
           'bond_rmsz', 'angle_rmsz', 'molprobity_score']

print(f'{"Metric":<20} {"Original":>10} {"Relaxed":>10} {"Change":>10}')
print("-" * 52)
for m in metrics:
    o = orig[m].mean()
    r = relax[m].mean()
    change = r - o
    print(f'{m:<20} {o:>10.2f} {r:>10.2f} {change:>+10.2f}')
