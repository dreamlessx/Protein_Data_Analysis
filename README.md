# Protein Structure Validation Pipeline

Systematic benchmarking of AlphaFold2 and Boltz1 predictions against crystallographic reference structures. Comprehensive validation of 6,820 protein structures derived from 20 experimentally resolved complexes using wwPDB-standard metrics.

## Overview

Machine learning methods for protein structure prediction have achieved remarkable accuracy, yet rigorous validation against established structural biology metrics remains essential. This project evaluates predicted structures using wwPDB validation criteria and geometry analysis tools standardized across the structural biology community.

### Key Findings

| Method | Rama Favored | Rota Favored | Clashscore | MP Score | Bond RMSZ | Angle RMSZ |
|--------|--------------|--------------|------------|----------|-----------|------------|
| Experimental | 97.8% | 90.5% | 4.7 | 1.18 | 0.54 | 1.05 |
| AlphaFold2 | 95.2% | 93.8% | 52.4 | 1.79 | 0.22 | 0.59 |
| Boltz1 | 98.9% | 98.4% | 7.9 | 1.11 | 0.38 | 0.69 |

AlphaFold2 and Boltz1 show near-ideal bond and angle geometry (RMSZ < 1.0) as expected from energy-minimized predictions. Boltz1 demonstrates superior backbone geometry with near-ideal Ramachandran distributions.

## Dataset

**20 Protein Complexes**

PDB identifiers: 1AK4, 1AKJ, 1AVX, 1AY7, 1AZS, 1BUH, 1BVN, 1E6E, 1EFN, 1EWY, 1EXB, 1F51, 1FCC, 1GHQ, 1GLA, 1HCF, 1JPS, 1K74, 1VFB, 2I25

**Structure Categories**

| Category | Count | Description |
|----------|-------|-------------|
| Experimental | 20 | X-ray crystallographic reference structures |
| AlphaFold2 | 100 | Five ranked predictions per complex |
| Boltz1 | 100 | Five model predictions per complex |
| Relaxed | 6,600 | Six Rosetta protocols applied to all structures |

**Rosetta Relaxation Protocols**

| Protocol | Energy Function | Coordinate Space |
|----------|-----------------|------------------|
| normal_ref15 | REF15 | Torsion |
| normal_beta | Beta Nov16 | Torsion |
| cartesian_ref15 | REF15 | Cartesian |
| cartesian_beta | Beta Nov16 | Cartesian |
| dualspace_ref15 | REF15 | Dualspace |
| dualspace_beta | Beta Nov16 | Dualspace |

## Validation Metrics

### MolProbity Analysis

Standard wwPDB validation metrics computed via the MolProbity toolchain:

| Metric | Description | Ideal Value |
|--------|-------------|-------------|
| Ramachandran Favored | Backbone dihedrals in favored regions | >98% |
| Ramachandran Outliers | Backbone dihedrals in disallowed regions | <0.2% |
| Rotamer Favored | Sidechain conformations in favored rotamers | >98% |
| Clashscore | Steric overlaps per 1000 atoms | <10 |
| MolProbity Score | Composite quality metric (lower is better) | <2.0 |
| C-beta Deviations | Deviations from ideal C-beta positions | 0 |
| Bond RMSZ | Root-mean-square Z-score for bond lengths | <1.0 |
| Angle RMSZ | Root-mean-square Z-score for bond angles | <1.0 |

### Geometry Validation

Structural integrity checks adapted from the PoseBusters methodology:

| Check | Description |
|-------|-------------|
| Backbone Connectivity | Peptide bond distances within tolerance |
| Bond Lengths | C-N peptide bonds 1.18-1.48 angstroms |
| Bond Angles | N-CA-C angles 100-120 degrees |
| Aromatic Planarity | Ring flatness RMSD below 0.1 angstroms |
| Peptide Planarity | Omega angles in cis or trans conformations |
| Chirality | L-amino acid stereochemistry verification |
| Steric Clashes | Van der Waals overlap detection |

## Repository Structure

```
Protein_Analysis/
├── proteins/                      # Structure files (Git LFS)
│   └── {PDB_ID}/
│       ├── {PDB_ID}.pdb          # Experimental reference
│       ├── AF/                    # AlphaFold predictions
│       ├── Boltz/                 # Boltz1 predictions
│       ├── {protocol}/            # Relaxed experimental structures
│       ├── relax/AF|Boltz/        # Relaxed predictions
│       └── analysis/              # Per-protein validation results
├── scripts/
│   ├── posebusters.py            # Geometry validation pipeline
│   ├── molprobity_extended.py    # Extended MolProbity metrics
│   └── run_validation_parallel.py # MolProbity validation pipeline
├── validation_results/
│   ├── posebusters_results.csv   # Boolean pass/fail per check
│   ├── posebusters_raw.csv       # Raw metric values
│   ├── molprobity_full.csv       # Complete MolProbity output (41 columns)
│   └── molprobity_extended.csv   # C-beta, omega, bond/angle RMSZ
└── requirements.txt              # Python dependencies
```

## Installation

### Python Environment

```bash
conda create -n protein_validation python=3.11
conda activate protein_validation
pip install -r requirements.txt
```

### External Dependencies

**MolProbity**

```bash
conda install -c conda-forge cctbx-base
```

**Rosetta** (for energy scoring)

Academic licenses available at: https://www.rosettacommons.org/software/license-and-download

**Reduce and Probe** (for clash analysis)

```bash
conda install -c conda-forge reduce probe
```

## Usage

```bash
# MolProbity validation (parallel)
python scripts/run_validation_parallel.py

# Extended metrics (C-beta, omega, bond/angle RMSZ)
python scripts/molprobity_extended.py

# PoseBusters geometry validation
python scripts/posebusters.py --workers 12
```

## Results

All validation results in `validation_results/` as CSV files. Per-protein summaries in `proteins/{PDB_ID}/analysis/`.

| File | Contents | Columns |
|------|----------|---------|
| molprobity_full.csv | Complete validation | 41 |
| molprobity_extended.csv | Extended geometry | 18 |
| posebusters_results.csv | Pass/fail checks | 12 |
| posebusters_raw.csv | Raw metrics | 25+ |

## References

Williams CJ, et al. (2018). MolProbity: More and better reference data for improved all-atom structure validation. *Protein Science* 27:293-315.

Buttenschoen M, Morris GM, Deane CM. (2024). PoseBusters: AI-based docking methods fail to generate physically valid ligand poses. *Chemical Science* 15:3130-3139.

Jumper J, et al. (2021). Highly accurate protein structure prediction with AlphaFold. *Nature* 596:583-589.

Wohlwend J, et al. (2024). Boltz1: Democratizing Biomolecular Interaction Modeling. *bioRxiv*.

## License

MIT License
