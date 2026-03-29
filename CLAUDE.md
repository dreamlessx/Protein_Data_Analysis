# Protein Data Analysis — AF2 vs Boltz1 Structural Validation

## What This Is
Post-prediction validation framework. Compares AlphaFold2 and Boltz1 structures
against crystallographic references using wwPDB-standard MolProbity metrics.
6,820 structures from 20 PDB complexes.

## Key Scripts
- `scripts/run_validation_parallel.py` — MolProbity validation (12-worker parallel)
  Metrics: Ramachandran, rotamers, clashscore, C-beta deviations, omega angles
- `scripts/molprobity_extended.py` — Extended geometry (RMSZ, bond/angle outliers)
- `scripts/posebusters.py` — Geometry checks (connectivity, planarity, chirality, clashes)

## Dependencies
pandas>=2.0, numpy>=1.24, scipy>=1.10, matplotlib>=3.7, seaborn>=0.12,
plotly>=5.15, tabulate>=0.9, jinja2>=3.1, scikit-learn>=1.3, tqdm>=4.65
System: cctbx-base, reduce, probe (MolProbity toolchain)

## Compute
- Curium (curium.csb.vanderbilt.edu): 36-core Xeon, 125GB RAM — ideal for parallel validation
- ACCRE batch partition for large-scale runs
- SSH: `ssh curium` (requires VPN) or `ssh accre`

## Reference Baselines
| Method | Rama Favored | Clashscore | MolProbity Score |
|--------|-------------|-----------|------------------|
| Experimental | 97.8% | 4.7 | 1.18 |
| AlphaFold2 | 95.2% | 52.4 | 1.79 |
| Boltz1 | 98.9% | 7.9 | 1.11 |

## Output
- validation_results/molprobity_full.csv (41 columns)
- validation_results/molprobity_extended.csv
- validation_results/posebusters_results.csv
- Per-protein summaries in proteins/{PDB_ID}/analysis/
