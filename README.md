# Protein_Data_Analysis: Phase 1 Validation Pilot (20 proteins)

**Phase 1 pilot for the BM5.5 docking-relaxation benchmark.** Established the validation methodology (MolProbity + PoseBusters geometry checks) on a 20-protein subset before scaling to the full 257-target dataset. Superseded by [`Protein_Relax_Pipeline`](https://github.com/dreamlessx/Protein_Relax_Pipeline) (Blue + canonical analysis) and [`Protein_Ideal`](https://github.com/dreamlessx/Protein_Ideal) (Green pipeline). Kept here for the historical record.

---

## Phase 1 vs full benchmark

| | Phase 1 (this repo) | Full BM5.5 benchmark |
|---|---|---|
| Proteins | 20 | 257 |
| Structures evaluated | 6,820 | 416,340 (Rosetta MolProbity) |
| Pipelines | 1 | 2 (Blue + Green) |
| Pre-Rosetta TM | not separately tracked | 12,065 rows |
| Post-Rosetta TM | not separately tracked | 92,700 rows |
| Rosetta total energy | not separately tracked | 416,340 rows (100% coverage) |
| Pre-Rosetta MolProbity | included in 6,820 | 13,364 rows |
| Locked snapshot | n/a | 2026-04-27a, qc_status = pass |
| Statistical analysis | descriptive | paired t + Wilcoxon + BH FDR + Cliff's d |
| PoseBusters checks | yes | not in scope (deprecated for the paper release) |

The Phase 1 pilot's primary contribution was establishing that MolProbity and PoseBusters geometry checks could be run at scale on Rosetta-relaxed predicted structures. The full benchmark dropped PoseBusters (out of paper scope) and shifted to AMBER + Rosetta as the primary relaxation matrix.

## Snapshot 2026-04-27a (full benchmark)

The companion repo `Protein_Relax_Pipeline` is fully populated under the locked snapshot:

- **rosetta_metrics**: 416,340 (locked)
- **prerosetta_metrics**: 13,364
- **tm_scores**: 104,765 (12,065 pre + 92,700 post)
- **rosetta_energy**: 416,340 (100% coverage of rosetta_metrics)
- **targets**: 257 with full metadata + parent_pdb_id for the 4 non-standard
- **qc_quarantine**: 0
- 0 orphans, 0 gaps, qc_status = pass

DB and raw TSVs in the [`db-2026-04-27a-supp`](https://github.com/dreamlessx/Protein_Relax_Pipeline/releases/tag/db-2026-04-27a-supp) Release.

Three findings (full numbers in `Protein_Relax_Pipeline/red_analysis/PAPER_FINDINGS.md`):
1. AMBER fixes local geometry without touching global fold (clashscore Cliff's d = -0.99, TM Cliff's d = -0.01).
2. Crystal structures carry the worst pre-Rosetta MolProbity (idealization artifact, not failure).
3. dualspace_beta wins integrated MolProbity at small TM cost; beta_nov16 dominates ref2015.

---

## Phase 1 dataset

20 BM5.5 protein-protein complexes:

```
1AK4, 1AKJ, 1AVX, 1AY7, 1AZS, 1BUH, 1BVN, 1E6E, 1EFN, 1EWY,
1EXB, 1F51, 1FCC, 1GHQ, 1GLA, 1HCF, 1JPS, 1K74, 1VFB, 2I25
```

Per complex: 1 experimental crystal + 5 AlphaFold 2.3.2 predictions + 5 Boltz-1 v0.4.1 predictions + 6 Rosetta protocols × every prediction = 341 structures. Total: 20 × 341 = 6,820 structures evaluated.

## Phase 1 results

| Method | Rama Favored | Rota Favored | Clashscore | MP Score | Bond RMSZ | Angle RMSZ |
|---|---|---|---|---|---|---|
| Experimental | 97.8% | 90.5% | 4.7 | 1.18 | 0.54 | 1.05 |
| AlphaFold 2 | 95.2% | 93.8% | 52.4 | 1.79 | 0.22 | 0.59 |
| Boltz 1 | 98.9% | 98.4% | 7.9 | 1.11 | 0.38 | 0.69 |

The Phase 1 pilot established that AlphaFold and Boltz predictions show near-ideal bond and angle geometry (RMSZ < 1.0) as expected from energy-minimized predictions, and that Boltz produces superior backbone geometry (Ramachandran-favored at 98.9% vs 95.2% for AlphaFold). These observations motivated the full benchmark's question: does relaxation (AMBER, Rosetta, or both) close the gap on local geometry, and at what cost to global fold?

The full benchmark's answer is yes (Finding 1: AMBER alone closes the clashscore gap with Cliff's d = -0.99 at near-zero TM cost).

---

## Validation methodology

### MolProbity metrics

Standard wwPDB validation, computed via the MolProbity toolchain.

| Metric | Description | Ideal |
|---|---|---|
| Ramachandran favored | Backbone dihedrals in favored regions | > 98% |
| Ramachandran outliers | Backbone dihedrals in disallowed regions | < 0.2% |
| Rotamer favored | Sidechain conformations in favored rotamers | > 98% |
| Clashscore | Steric overlaps per 1000 atoms | < 10 |
| MolProbity score | Composite quality metric | < 2.0 |
| C-beta deviations | Deviations from ideal C-beta positions | 0 |
| Bond RMSZ | RMS Z-score for bond lengths | < 1.0 |
| Angle RMSZ | RMS Z-score for bond angles | < 1.0 |

### PoseBusters geometry checks (Phase 1 only)

Adapted from Buttenschoen 2024.

| Check | Description |
|---|---|
| Backbone connectivity | Peptide bond distances within tolerance |
| Bond lengths | C-N peptide bonds 1.18-1.48 Å |
| Bond angles | N-CA-C angles 100-120° |
| Aromatic planarity | Ring flatness RMSD < 0.1 Å |
| Peptide planarity | Omega angles in cis or trans |
| Chirality | L-amino acid stereochemistry |
| Steric clashes | Van der Waals overlap detection |

The full benchmark dropped PoseBusters (out of paper scope) but the methodology remains in this repo for future use.

---

## Repository layout

```
Protein_Data_Analysis/
├── proteins/                  Per-protein structure files (Git LFS)
│   └── {PDB_ID}/
│       ├── {PDB_ID}.pdb      Experimental reference
│       ├── AF/                AlphaFold predictions (5 ranked)
│       ├── Boltz/             Boltz-1 predictions (5 models)
│       ├── {protocol}/        Rosetta-relaxed experimental
│       ├── relax/AF|Boltz/    Rosetta-relaxed predictions
│       └── analysis/          Per-protein validation
├── scripts/
│   ├── posebusters.py            Geometry checks
│   ├── molprobity_extended.py    Extended MolProbity (C-beta, omega, RMSZ)
│   └── run_validation_parallel.py  Parallel MolProbity runner
├── validation_results/
│   ├── molprobity_full.csv       Complete MolProbity output (41 columns)
│   ├── molprobity_extended.csv   Extended geometry (18 columns)
│   ├── posebusters_results.csv   Boolean pass/fail per check
│   └── posebusters_raw.csv       Raw metric values
└── requirements.txt
```

## Installation

```bash
conda create -n protein_validation python=3.11
conda activate protein_validation
pip install -r requirements.txt

# MolProbity via cctbx
conda install -c conda-forge cctbx-base

# Reduce + Probe for clash analysis
conda install -c conda-forge reduce probe
```

Rosetta licenses available from RosettaCommons.

## Usage

```bash
# MolProbity validation (parallel)
python scripts/run_validation_parallel.py

# Extended geometry metrics (C-beta, omega, bond/angle RMSZ)
python scripts/molprobity_extended.py

# PoseBusters geometry checks
python scripts/posebusters.py --workers 12
```

Outputs in `validation_results/`. Per-protein summaries in `proteins/{PDB_ID}/analysis/`.

## References

1. Williams, C.J., Headd, J.J., Moriarty, N.W. et al. MolProbity: More and better reference data for improved all-atom structure validation. *Protein Sci.* 27, 293-315 (2018).
2. Buttenschoen, M., Morris, G.M., Deane, C.M. PoseBusters: AI-based docking methods fail to generate physically valid ligand poses. *Chem. Sci.* 15, 3130-3139 (2024).
3. Jumper, J., Evans, R., Pritzel, A. et al. Highly accurate protein structure prediction with AlphaFold. *Nature* 596, 583-589 (2021).
4. Wohlwend, J., Corso, G., Passaro, S. et al. Boltz-1: Democratizing Biomolecular Interaction Modeling. *bioRxiv* (2024).

## License

MIT.

---

*Phase 1 pilot. Full benchmark at [`Protein_Relax_Pipeline`](https://github.com/dreamlessx/Protein_Relax_Pipeline), snapshot 2026-04-27a, locked at 100.000% with 416,340 Rosetta MolProbity rows.*
