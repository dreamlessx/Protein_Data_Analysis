# Protein Data Analysis: AF2 vs Boltz1 Structural Validation

## What This Is
Post-prediction validation framework. Compares AlphaFold2 and Boltz1 structures
against crystallographic references using wwPDB-standard MolProbity metrics.
6,820 structures from 20 PDB complexes.

## Key Scripts
- `scripts/run_validation_parallel.py` runs MolProbity validation (12-worker parallel).
  Metrics: Ramachandran, rotamers, clashscore, C-beta deviations, omega angles.
- `scripts/molprobity_extended.py` computes extended geometry (RMSZ, bond/angle outliers).
- `scripts/posebusters.py` runs geometry checks (connectivity, planarity, chirality, clashes).

## Dependencies
pandas>=2.0, numpy>=1.24, scipy>=1.10, matplotlib>=3.7, seaborn>=0.12,
plotly>=5.15, tabulate>=0.9, jinja2>=3.1, scikit-learn>=1.3, tqdm>=4.65
System: cctbx-base, reduce, probe (MolProbity toolchain)

## Compute
- Curium (curium.csb.vanderbilt.edu): 36-core Xeon, 125GB RAM, ideal for parallel validation.
- ACCRE batch partition for large-scale runs.
- SSH: `ssh curium` (requires VPN) or `ssh accre`.

## Reference Baselines
| Source | Rama Favored | Rota Favored | Clashscore | MP Score | Bond RMSZ | Angle RMSZ |
| --- | --- | --- | --- | --- | --- | --- |
| Crystal | 97.8% | 90.5% | 4.7 | 1.18 | 0.54 | 1.05 |
| AF2 | 95.2% | 93.8% | 52.4 | 1.79 | 0.22 | 0.59 |
| Boltz1 | 98.9% | 98.4% | 7.9 | 1.11 | 0.38 | 0.69 |

Companion repo dreamlessx/Protein_Relax_Pipeline reached 100.000% Rosetta MolProbity coverage on snapshot 2026-04-27a (416,340/416,340 rows).

## Output
- validation_results/molprobity_full.csv (41 columns)
- validation_results/molprobity_extended.csv
- validation_results/posebusters_results.csv
- Per-protein summaries in proteins/{PDB_ID}/analysis/
