# Protein Structure Validation Pipeline

Systematic benchmarking of AlphaFold2 and Boltz1 predictions against crystallographic reference structures. This repository contains validation scripts, analysis pipelines, and comprehensive results for 6,820 protein structures derived from 20 experimentally resolved complexes.

## Overview

Machine learning methods for protein structure prediction have achieved remarkable accuracy, yet rigorous validation against established structural biology metrics remains essential. This project evaluates predicted structures using wwPDB validation criteria and geometry analysis tools standardized across the structural biology community.

### Key Findings

| Method | Ramachandran Favored | Rotamer Favored | MolProbity Score |
|--------|---------------------|-----------------|------------------|
| Experimental (X ray) | 97.8% | 90.5% | 1.18 |
| AlphaFold2 | 95.2% | 93.8% | 1.79 |
| Boltz1 | 98.9% | 98.4% | 1.11 |

Boltz1 demonstrates superior backbone geometry with near ideal Ramachandran distributions. AlphaFold2 shows competitive performance on composite quality metrics. Rosetta relaxation protocols improve geometric quality across all prediction methods.

## Dataset

**20 Protein Complexes**

PDB identifiers: 1AK4, 1AKJ, 1AVX, 1AY7, 1AZS, 1BUH, 1BVN, 1E6E, 1EFN, 1EWY, 1EXB, 1F51, 1FCC, 1GHQ, 1GLA, 1HCF, 1JPS, 1K74, 1VFB, 2I25

**Structure Categories**
| Category | Count | Description |
|----------|-------|-------------|
| Experimental | 20 | X ray crystallographic reference structures |
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
| C beta Deviations | Deviations from ideal C beta positions | 0 |

### Geometry Validation
Structural integrity checks adapted from the PoseBusters methodology:

| Check | Description |
|-------|-------------|
| Backbone Connectivity | Peptide bond distances within tolerance |
| Bond Lengths | C N peptide bonds 1.18 to 1.48 angstroms |
| Bond Angles | N CA C angles 100 to 120 degrees |
| Aromatic Planarity | Ring flatness RMSD below 0.1 angstroms |
| Peptide Planarity | Omega angles in cis or trans conformations |
| Chirality | L amino acid stereochemistry verification |
| Steric Clashes | Van der Waals overlap detection |

## Repository Structure

```
Protein_Analysis/
├── proteins/                      # Structure files (not tracked)
│   └── {PDB_ID}/
│       ├── {PDB_ID}.pdb          # Experimental reference
│       ├── AF/                    # AlphaFold predictions
│       ├── Boltz/                 # Boltz1 predictions
│       ├── {protocol}/            # Relaxed experimental structures
│       └── relax/AF|Boltz/        # Relaxed predictions
├── scripts/
│   ├── posebusters.py            # Geometry validation pipeline
│   ├── molprobity_extended.py    # Extended MolProbity metrics
│   └── run_validation_parallel.py # MolProbity validation pipeline
├── validation_results/
│   ├── posebusters_results.csv   # Boolean pass/fail per check
│   ├── posebusters_raw.csv       # Raw metric values
│   ├── molprobity_full.csv       # Complete MolProbity output
│   └── molprobity_extended.csv   # C beta and omega distributions
└── requirements.txt              # Python dependencies
```

## Installation

### Python Environment

```bash
# Create conda environment
conda create -n protein_validation python=3.11
conda activate protein_validation

# Install dependencies
pip install -r requirements.txt
```

### External Dependencies

**MolProbity**

MolProbity provides the validation tools used by the worldwide Protein Data Bank. Installation via conda is recommended.

```bash
conda install -c conda-forge cctbx-base
```

Documentation: https://github.com/rlabduke/MolProbity

**Rosetta**

Rosetta is required for energy scoring and structure relaxation. Academic licenses are available at no cost.

Download: https://www.rosettacommons.org/software/license-and-download

**Reduce and Probe**

Required for hydrogen placement and clash analysis.

```bash
# Via conda
conda install -c conda-forge reduce probe
```

Reduce: https://github.com/rlabduke/reduce
Probe: https://github.com/rlabduke/probe

### Python Dependencies

| Package | Version | Purpose |
|---------|---------|---------|
| pandas | ≥2.0.0 | Data manipulation |
| numpy | ≥1.24.0 | Numerical computing |
| scipy | ≥1.10.0 | Statistical analysis |
| statsmodels | ≥0.14.0 | Statistical modeling |
| scikit learn | ≥1.3.0 | Machine learning utilities |
| matplotlib | ≥3.7.0 | Visualization |
| seaborn | ≥0.12.0 | Statistical visualization |
| plotly | ≥5.15.0 | Interactive plots |
| tabulate | ≥0.9.0 | Table formatting |
| jinja2 | ≥3.1.0 | Template rendering |
| openpyxl | ≥3.1.0 | Excel export |
| tqdm | ≥4.65.0 | Progress bars |

## Usage

### Run Full Validation Pipeline

```bash
# MolProbity validation (parallel execution)
python scripts/run_validation_parallel.py

# Extended geometry metrics
python scripts/molprobity_extended.py

# PoseBusters style validation
python scripts/posebusters.py --workers 12
```

### Command Line Options

```bash
# PoseBusters validation
python scripts/posebusters.py --help

Options:
  --no-energy     Skip Rosetta energy scoring
  --limit N       Process first N structures only
  -j, --workers   Number of parallel workers (default: CPU count)
```

## Results

Validation results are stored in `validation_results/` as CSV files suitable for statistical analysis. Per protein summaries are written to `proteins/{PDB_ID}/analysis/`.

### Output Files

| File | Contents |
|------|----------|
| posebusters_results.csv | Boolean pass/fail for each validation check |
| posebusters_raw.csv | Raw numerical values for all metrics |
| molprobity_full.csv | Complete MolProbity analysis |
| molprobity_extended.csv | C beta deviation and omega angle statistics |

## References

**Validation Methods**

Williams CJ, Headd JJ, Moriarty NW, et al. (2018). MolProbity: More and better reference data for improved all atom structure validation. Protein Science 27:293-315.

Buttenschoen M, Morris GM, Deane CM. (2024). PoseBusters: AI based docking methods fail to generate physically valid ligand poses or generalise to novel sequences. Chemical Science 15:3130-3139.

**Structure Prediction**

Jumper J, Evans R, Pritzel A, et al. (2021). Highly accurate protein structure prediction with AlphaFold. Nature 596:583-589.

Wohlwend J, et al. (2024). Boltz1: Democratizing Biomolecular Interaction Modeling. bioRxiv.

**Energy Functions**

Alford RF, Leaver-Fay A, Jeliazkov JR, et al. (2017). The Rosetta All Atom Energy Function for Macromolecular Modeling and Design. Journal of Chemical Theory and Computation 13:3031-3048.

## License

MIT License. See LICENSE file for details.

## Author

Mudit Agar
