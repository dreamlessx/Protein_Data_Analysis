#!/usr/bin/env python3
"""
run_full_validation_per_protein.py
Complete validation pipeline with per-protein organization

Output structure:
  proteins/{PDB}/analysis/
    molprobity_results.csv      - All MolProbity metrics
    posebusters_results.csv     - All PoseBusters metrics
    VALIDATION_SUMMARY.md       - Human-readable summary

Usage:
    conda activate molprobity
    python scripts/run_full_validation_per_protein.py
"""

import subprocess
import os
import sys
import re
import gzip
import shutil
import numpy as np
import pandas as pd
from pathlib import Path
from datetime import datetime
from tqdm import tqdm
import warnings
import csv
import json

warnings.filterwarnings('ignore')

PROJECT_DIR = Path(__file__).parent.parent
PROTEINS_DIR = PROJECT_DIR / "proteins"
TEMP_DIR = PROJECT_DIR / "temp_validation"
REDUCE_HET_DICT = Path.home() / "miniconda3/envs/molprobity/share/reduce/reduce_wwPDB_het_dict.txt"

TEMP_DIR.mkdir(parents=True, exist_ok=True)


def find_all_structures_for_protein(protein_dir):
    """Find ALL PDB structures for a single protein, including nested relaxed."""
    structures = []
    protein_id = protein_dir.name

    # 1. Experimental original
    exp_pdb = protein_dir / f"{protein_id}.pdb"
    if exp_pdb.exists():
        structures.append({
            'path': exp_pdb,
            'protein': protein_id,
            'category': 'Experimental',
            'subcategory': 'original',
            'model': 'exp',
            'source': 'PDB',
            'relaxation': 'none'
        })

    # 2. AlphaFold raw predictions
    af_dir = protein_dir / "AF"
    if af_dir.exists():
        for pdb in sorted(af_dir.glob("ranked_*.pdb")):
            model_num = pdb.stem.split('_')[1]
            structures.append({
                'path': pdb,
                'protein': protein_id,
                'category': 'AlphaFold',
                'subcategory': 'raw',
                'model': f"ranked_{model_num}",
                'source': 'AlphaFold3',
                'relaxation': 'none'
            })

    # 3. Boltz raw predictions
    boltz_dir = protein_dir / "Boltz"
    if boltz_dir.exists():
        for pdb in sorted(boltz_dir.glob("boltz_input_model_*.pdb")):
            model_num = pdb.stem.split('_')[-1]
            structures.append({
                'path': pdb,
                'protein': protein_id,
                'category': 'Boltz',
                'subcategory': 'raw',
                'model': f"model_{model_num}",
                'source': 'Boltz1',
                'relaxation': 'none'
            })

    # 4. Experimental relaxed structures
    relax_protocols = ['cartesian_beta', 'cartesian_ref15', 'dualspace_beta',
                       'dualspace_ref15', 'normal_beta', 'normal_ref15']

    for protocol in relax_protocols:
        relax_dir = protein_dir / protocol
        if relax_dir.exists():
            for pdb_gz in sorted(relax_dir.glob(f"{protein_id}_r*.pdb.gz")):
                rep = pdb_gz.stem.replace('.pdb', '').split('_r')[-1]
                structures.append({
                    'path': pdb_gz,
                    'protein': protein_id,
                    'category': 'Experimental',
                    'subcategory': f'relaxed_{protocol}',
                    'model': f"r{rep}",
                    'source': 'PDB',
                    'relaxation': protocol,
                    'compressed': True
                })

    # 5. AF/Boltz relaxed structures (in relax/ directory)
    relax_base = protein_dir / "relax"
    if relax_base.exists():
        for cat_name in ['AF', 'Boltz']:
            cat_dir = relax_base / cat_name
            if not cat_dir.exists():
                continue

            category = 'AlphaFold' if cat_name == 'AF' else 'Boltz'

            # Walk through all subdirectories to find .pdb.gz files
            for pdb_gz in cat_dir.rglob("*.pdb.gz"):
                # Parse the path to extract model and protocol info
                rel_path = pdb_gz.relative_to(cat_dir)
                parts = list(rel_path.parts)

                # First part is usually the model (ranked_0, ranked_1, etc.)
                base_model = parts[0] if parts else "unknown"

                # Find protocol (one of the relax_protocols)
                protocol = "unknown"
                for p in parts:
                    if p in relax_protocols:
                        protocol = p
                        break

                # Create a unique model identifier
                stem = pdb_gz.stem.replace('.pdb', '')

                structures.append({
                    'path': pdb_gz,
                    'protein': protein_id,
                    'category': category,
                    'subcategory': f'relaxed_{protocol}',
                    'model': stem,
                    'source': 'AlphaFold3' if cat_name == 'AF' else 'Boltz1',
                    'relaxation': protocol,
                    'base_model': base_model,
                    'compressed': True
                })

    return structures


def decompress_if_needed(structure):
    """Decompress gzipped PDB file if needed."""
    if structure.get('compressed'):
        gz_path = structure['path']
        safe_name = f"{structure['protein']}_{structure['category']}_{structure['model'][:50]}.pdb"
        safe_name = re.sub(r'[^\w\-_.]', '_', safe_name)
        pdb_path = TEMP_DIR / safe_name

        with gzip.open(gz_path, 'rb') as f_in:
            with open(pdb_path, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        return pdb_path
    return structure['path']


# =============================================================================
# MOLPROBITY FUNCTIONS
# =============================================================================

def run_ramalyze(pdb_path):
    """Run Ramachandran analysis."""
    try:
        result = subprocess.run(
            ['molprobity.ramalyze', str(pdb_path)],
            capture_output=True, text=True, timeout=120
        )
        output = result.stdout + result.stderr

        favored = allowed = outliers = total = 0
        for line in output.split('\n'):
            if 'SUMMARY:' in line:
                match = re.search(r'(\d+) Favored, (\d+) Allowed, (\d+) Outlier.* out of (\d+)', line)
                if match:
                    favored, allowed, outliers, total = map(int, match.groups())

        if total > 0:
            return {
                'rama_favored': favored,
                'rama_allowed': allowed,
                'rama_outliers': outliers,
                'rama_total': total,
                'rama_favored_pct': round(100 * favored / total, 2),
                'rama_outliers_pct': round(100 * outliers / total, 2)
            }
    except Exception:
        pass
    return {}


def run_cbetadev(pdb_path):
    """Run C-beta deviation analysis."""
    try:
        result = subprocess.run(
            ['molprobity.cbetadev', str(pdb_path)],
            capture_output=True, text=True, timeout=120
        )
        for line in (result.stdout + result.stderr).split('\n'):
            if 'SUMMARY:' in line:
                match = re.search(r'(\d+) C-beta deviation', line)
                if match:
                    return {'cbeta_deviations': int(match.group(1))}
    except Exception:
        pass
    return {}


def run_omegalyze(pdb_path):
    """Run omega angle analysis."""
    try:
        result = subprocess.run(
            ['molprobity.omegalyze', str(pdb_path)],
            capture_output=True, text=True, timeout=120
        )
        output = result.stdout + result.stderr

        cis_pro = cis_general = twisted_pro = twisted_general = 0
        for line in output.split('\n'):
            if 'SUMMARY:' in line:
                m = re.search(r'(\d+)\s+cis\s+prolines?\s+out\s+of', line, re.I)
                if m: cis_pro = int(m.group(1))
                m = re.search(r'(\d+)\s+twisted\s+prolines?\s+out\s+of', line, re.I)
                if m: twisted_pro = int(m.group(1))
                m = re.search(r'(\d+)\s+other\s+cis\s+residues?\s+out\s+of', line, re.I)
                if m: cis_general = int(m.group(1))
                m = re.search(r'(\d+)\s+other\s+twisted\s+residues?\s+out\s+of', line, re.I)
                if m: twisted_general = int(m.group(1))

        return {
            'omega_cis_proline': cis_pro,
            'omega_cis_general': cis_general,
            'omega_twisted': twisted_pro + twisted_general
        }
    except Exception:
        pass
    return {}


def run_clashscore(pdb_path):
    """Calculate clashscore using reduce + probe."""
    try:
        env = os.environ.copy()
        env['REDUCE_HET_DICT'] = str(REDUCE_HET_DICT)

        reduce_result = subprocess.run(
            ['reduce', '-build', str(pdb_path)],
            capture_output=True, text=True, timeout=120, env=env
        )

        pdb_with_h = reduce_result.stdout
        if not pdb_with_h or 'ATOM' not in pdb_with_h:
            return {}

        temp_h_pdb = TEMP_DIR / f"temp_H_{os.getpid()}.pdb"
        with open(temp_h_pdb, 'w') as f:
            f.write(pdb_with_h)

        probe_path = '/private/tmp/probe/probe'
        if not os.path.exists(probe_path):
            probe_path = 'probe'

        probe_result = subprocess.run(
            [probe_path, '-4H', '-mc', '-self', 'ALL', '-unformated', str(temp_h_pdb)],
            capture_output=True, text=True, timeout=120
        )

        temp_h_pdb.unlink(missing_ok=True)

        clash_pairs = set()
        for line in probe_result.stdout.split('\n'):
            if line.startswith(':') and ':bo:' in line:
                parts = line.split(':')
                if len(parts) >= 5:
                    pair = tuple(sorted([parts[3].strip(), parts[4].strip()]))
                    clash_pairs.add(pair)

        clash_count = len(clash_pairs)
        atom_count = pdb_with_h.count('\nATOM ') + pdb_with_h.count('\nHETATM ')

        if atom_count > 0:
            clashscore = (clash_count * 1000) / atom_count
            return {
                'clashscore': round(clashscore, 2),
                'clash_count': clash_count,
                'atom_count': atom_count
            }
    except Exception:
        pass
    return {}


def run_molprobity(pdb_path):
    """Run all MolProbity analyses."""
    results = {}
    results.update(run_ramalyze(pdb_path))
    results.update(run_cbetadev(pdb_path))
    results.update(run_omegalyze(pdb_path))
    results.update(run_clashscore(pdb_path))
    return results


# =============================================================================
# POSEBUSTERS FUNCTIONS
# =============================================================================

def parse_pdb(pdb_path):
    """Parse PDB file and extract atom coordinates."""
    atoms = []
    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                try:
                    atoms.append({
                        'name': line[12:16].strip(),
                        'resname': line[17:20].strip(),
                        'chain': line[21].strip(),
                        'resseq': int(line[22:26]),
                        'x': float(line[30:38]),
                        'y': float(line[38:46]),
                        'z': float(line[46:54]),
                        'element': line[76:78].strip() if len(line) > 76 else ''
                    })
                except (ValueError, IndexError):
                    continue
    return atoms


def distance(a1, a2):
    return np.sqrt((a1['x']-a2['x'])**2 + (a1['y']-a2['y'])**2 + (a1['z']-a2['z'])**2)


def dihedral(a1, a2, a3, a4):
    b1 = np.array([a2['x']-a1['x'], a2['y']-a1['y'], a2['z']-a1['z']])
    b2 = np.array([a3['x']-a2['x'], a3['y']-a2['y'], a3['z']-a2['z']])
    b3 = np.array([a4['x']-a3['x'], a4['y']-a3['y'], a4['z']-a3['z']])
    n1 = np.cross(b1, b2)
    n2 = np.cross(b2, b3)
    m1 = np.cross(n1, b2/np.linalg.norm(b2))
    return np.degrees(np.arctan2(np.dot(m1, n2), np.dot(n1, n2)))


def run_posebusters(pdb_path):
    """Run PoseBusters-style validation tests."""
    try:
        atoms = parse_pdb(pdb_path)
        if len(atoms) == 0:
            return {'error': 'No atoms found'}

        results = {'n_atoms': len(atoms)}

        # Group atoms by residue
        residues = {}
        for atom in atoms:
            key = (atom['chain'], atom['resseq'])
            if key not in residues:
                residues[key] = {'resname': atom['resname']}
            residues[key][atom['name']] = atom

        results['n_residues'] = len(residues)

        # Bond length check (C-N peptide bonds)
        sorted_keys = sorted(residues.keys())
        bond_lengths = []
        bad_bonds = 0

        for i in range(len(sorted_keys) - 1):
            if sorted_keys[i][0] != sorted_keys[i + 1][0]:
                continue
            if sorted_keys[i + 1][1] - sorted_keys[i][1] != 1:
                continue

            res1 = residues[sorted_keys[i]]
            res2 = residues[sorted_keys[i + 1]]

            if 'C' in res1 and 'N' in res2:
                d = distance(res1['C'], res2['N'])
                bond_lengths.append(d)
                if d < 1.2 or d > 1.5:
                    bad_bonds += 1

        results['n_peptide_bonds'] = len(bond_lengths)
        results['bad_bonds'] = bad_bonds
        results['bond_length_mean'] = round(np.mean(bond_lengths), 4) if bond_lengths else 0
        results['bond_length_std'] = round(np.std(bond_lengths), 4) if bond_lengths else 0

        # Omega angle check
        cis_count = trans_count = twisted_count = cis_proline = 0

        for i in range(len(sorted_keys) - 1):
            if sorted_keys[i][0] != sorted_keys[i + 1][0]:
                continue

            res1 = residues[sorted_keys[i]]
            res2 = residues[sorted_keys[i + 1]]

            if all(a in res1 for a in ['CA', 'C']) and all(a in res2 for a in ['N', 'CA']):
                omega = dihedral(res1['CA'], res1['C'], res2['N'], res2['CA'])
                abs_omega = abs(omega)

                if abs_omega < 30:
                    cis_count += 1
                    if res2.get('resname') == 'PRO':
                        cis_proline += 1
                elif abs_omega > 150:
                    trans_count += 1
                else:
                    twisted_count += 1

        results['pb_cis_peptides'] = cis_count
        results['pb_trans_peptides'] = trans_count
        results['pb_twisted_peptides'] = twisted_count
        results['pb_cis_proline'] = cis_proline

        # Chirality check
        d_amino_acids = 0
        chiral_checked = 0

        for key, res in residues.items():
            if res.get('resname') == 'GLY':
                continue
            if all(a in res for a in ['N', 'CA', 'C', 'CB']):
                chiral_checked += 1
                chi = dihedral(res['N'], res['CA'], res['C'], res['CB'])
                if chi < -30:
                    d_amino_acids += 1

        results['chiral_residues_checked'] = chiral_checked
        results['d_amino_acids'] = d_amino_acids

        # Backbone continuity
        backbone_atoms = ['N', 'CA', 'C', 'O']
        missing_atoms = 0
        chain_breaks = 0

        for key, res in residues.items():
            for bb in backbone_atoms:
                if bb not in res:
                    missing_atoms += 1

        for i in range(len(sorted_keys) - 1):
            if sorted_keys[i][0] != sorted_keys[i + 1][0]:
                continue

            res1 = residues[sorted_keys[i]]
            res2 = residues[sorted_keys[i + 1]]

            if 'C' in res1 and 'N' in res2:
                if distance(res1['C'], res2['N']) > 2.0:
                    chain_breaks += 1

        results['missing_backbone_atoms'] = missing_atoms
        results['chain_breaks'] = chain_breaks

        # Disulfide bonds
        cys_sg = [a for a in atoms if a['resname'] == 'CYS' and a['name'] == 'SG']
        n_disulfides = 0
        bad_disulfides = 0

        for i in range(len(cys_sg)):
            for j in range(i + 1, len(cys_sg)):
                d = distance(cys_sg[i], cys_sg[j])
                if 1.8 < d < 2.5:
                    n_disulfides += 1
                    if d < 1.95 or d > 2.15:
                        bad_disulfides += 1

        results['n_disulfides'] = n_disulfides
        results['bad_disulfide_geometry'] = bad_disulfides

        # Overall quality flags
        results['pb_bond_pass'] = bad_bonds == 0
        results['pb_omega_pass'] = twisted_count == 0
        results['pb_chiral_pass'] = d_amino_acids == 0
        results['pb_backbone_pass'] = chain_breaks == 0 and missing_atoms == 0
        results['pb_all_pass'] = all([
            bad_bonds == 0,
            twisted_count == 0,
            d_amino_acids == 0,
            chain_breaks == 0
        ])

        return results

    except Exception as e:
        return {'error': str(e)}


# =============================================================================
# VALIDATION SUMMARY GENERATION
# =============================================================================

def generate_validation_summary(protein_id, mp_df, pb_df, output_path):
    """Generate a comprehensive validation summary markdown file."""

    lines = [
        f"# Validation Summary: {protein_id}",
        f"",
        f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        f"",
        f"---",
        f"",
        f"## Overview",
        f"",
        f"| Metric | Value |",
        f"|--------|-------|",
        f"| Total Structures | {len(mp_df)} |",
    ]

    # Count by category
    if 'category' in mp_df.columns:
        for cat in mp_df['category'].unique():
            count = len(mp_df[mp_df['category'] == cat])
            lines.append(f"| {cat} Structures | {count} |")

    lines.extend([
        f"",
        f"---",
        f"",
        f"## MolProbity Results Summary",
        f"",
    ])

    # MolProbity summary by category
    if 'rama_favored_pct' in mp_df.columns:
        lines.append("### Ramachandran Analysis")
        lines.append("")
        lines.append("| Category | Subcategory | Avg Favored % | Avg Outliers % | Count |")
        lines.append("|----------|-------------|---------------|----------------|-------|")

        for cat in ['Experimental', 'AlphaFold', 'Boltz']:
            cat_df = mp_df[mp_df['category'] == cat]
            if len(cat_df) == 0:
                continue

            for subcat in cat_df['subcategory'].unique():
                sub_df = cat_df[cat_df['subcategory'] == subcat]
                if 'rama_favored_pct' in sub_df.columns:
                    avg_fav = sub_df['rama_favored_pct'].mean()
                    avg_out = sub_df['rama_outliers_pct'].mean() if 'rama_outliers_pct' in sub_df.columns else 0
                    lines.append(f"| {cat} | {subcat} | {avg_fav:.2f} | {avg_out:.2f} | {len(sub_df)} |")

        lines.append("")

    if 'clashscore' in mp_df.columns:
        lines.append("### Clashscore Analysis")
        lines.append("")
        lines.append("| Category | Subcategory | Avg Clashscore | Min | Max | Count |")
        lines.append("|----------|-------------|----------------|-----|-----|-------|")

        for cat in ['Experimental', 'AlphaFold', 'Boltz']:
            cat_df = mp_df[mp_df['category'] == cat]
            if len(cat_df) == 0:
                continue

            for subcat in cat_df['subcategory'].unique():
                sub_df = cat_df[cat_df['subcategory'] == subcat]
                if 'clashscore' in sub_df.columns and sub_df['clashscore'].notna().any():
                    avg_cs = sub_df['clashscore'].mean()
                    min_cs = sub_df['clashscore'].min()
                    max_cs = sub_df['clashscore'].max()
                    lines.append(f"| {cat} | {subcat} | {avg_cs:.2f} | {min_cs:.2f} | {max_cs:.2f} | {len(sub_df)} |")

        lines.append("")

    lines.extend([
        f"---",
        f"",
        f"## PoseBusters Results Summary",
        f"",
    ])

    if 'pb_all_pass' in pb_df.columns:
        lines.append("### Pass Rates by Category")
        lines.append("")
        lines.append("| Category | Subcategory | All Tests Pass | Bond Pass | Omega Pass | Chiral Pass | Count |")
        lines.append("|----------|-------------|----------------|-----------|------------|-------------|-------|")

        for cat in ['Experimental', 'AlphaFold', 'Boltz']:
            cat_df = pb_df[pb_df['category'] == cat]
            if len(cat_df) == 0:
                continue

            for subcat in cat_df['subcategory'].unique():
                sub_df = cat_df[cat_df['subcategory'] == subcat]
                all_pass = 100 * sub_df['pb_all_pass'].sum() / len(sub_df) if len(sub_df) > 0 else 0
                bond_pass = 100 * sub_df['pb_bond_pass'].sum() / len(sub_df) if 'pb_bond_pass' in sub_df.columns else 0
                omega_pass = 100 * sub_df['pb_omega_pass'].sum() / len(sub_df) if 'pb_omega_pass' in sub_df.columns else 0
                chiral_pass = 100 * sub_df['pb_chiral_pass'].sum() / len(sub_df) if 'pb_chiral_pass' in sub_df.columns else 0
                lines.append(f"| {cat} | {subcat} | {all_pass:.1f}% | {bond_pass:.1f}% | {omega_pass:.1f}% | {chiral_pass:.1f}% | {len(sub_df)} |")

        lines.append("")

    # Best structures
    lines.extend([
        f"---",
        f"",
        f"## Best Structures",
        f"",
    ])

    if 'rama_favored_pct' in mp_df.columns and len(mp_df) > 0:
        best_rama = mp_df.loc[mp_df['rama_favored_pct'].idxmax()]
        lines.append(f"**Best Ramachandran:** {best_rama['model']} ({best_rama['category']}/{best_rama['subcategory']}) - {best_rama['rama_favored_pct']:.2f}% favored")
        lines.append("")

    if 'clashscore' in mp_df.columns and mp_df['clashscore'].notna().any():
        best_clash = mp_df.loc[mp_df['clashscore'].idxmin()]
        lines.append(f"**Best Clashscore:** {best_clash['model']} ({best_clash['category']}/{best_clash['subcategory']}) - {best_clash['clashscore']:.2f}")
        lines.append("")

    lines.extend([
        f"---",
        f"",
        f"## Files",
        f"",
        f"- `molprobity_results.csv` - Full MolProbity validation data",
        f"- `posebusters_results.csv` - Full PoseBusters validation data",
        f"",
    ])

    with open(output_path, 'w') as f:
        f.write('\n'.join(lines))


# =============================================================================
# MAIN
# =============================================================================

def process_protein(protein_dir):
    """Process a single protein: run all validations and save results."""
    protein_id = protein_dir.name
    analysis_dir = protein_dir / "analysis"
    analysis_dir.mkdir(parents=True, exist_ok=True)

    # Find all structures
    structures = find_all_structures_for_protein(protein_dir)

    if not structures:
        return protein_id, 0, "No structures found"

    mp_results = []
    pb_results = []

    for struct in structures:
        pdb_path = decompress_if_needed(struct)

        try:
            # Run MolProbity
            mp_data = run_molprobity(pdb_path)
            mp_data.update({
                'protein': struct['protein'],
                'category': struct['category'],
                'subcategory': struct['subcategory'],
                'model': struct['model'],
                'source': struct['source'],
                'relaxation': struct['relaxation']
            })
            mp_results.append(mp_data)

            # Run PoseBusters
            pb_data = run_posebusters(pdb_path)
            pb_data.update({
                'protein': struct['protein'],
                'category': struct['category'],
                'subcategory': struct['subcategory'],
                'model': struct['model'],
                'source': struct['source'],
                'relaxation': struct['relaxation']
            })
            pb_results.append(pb_data)

        finally:
            if struct.get('compressed') and pdb_path.exists():
                pdb_path.unlink(missing_ok=True)

    # Save results
    mp_df = pd.DataFrame(mp_results)
    pb_df = pd.DataFrame(pb_results)

    # Reorder columns
    meta_cols = ['protein', 'category', 'subcategory', 'model', 'source', 'relaxation']

    mp_cols = meta_cols + [c for c in mp_df.columns if c not in meta_cols]
    mp_df = mp_df[mp_cols]

    pb_cols = meta_cols + [c for c in pb_df.columns if c not in meta_cols]
    pb_df = pb_df[pb_cols]

    mp_df.to_csv(analysis_dir / "molprobity_results.csv", index=False)
    pb_df.to_csv(analysis_dir / "posebusters_results.csv", index=False)

    # Generate summary
    generate_validation_summary(protein_id, mp_df, pb_df, analysis_dir / "VALIDATION_SUMMARY.md")

    return protein_id, len(structures), "Success"


def main():
    print("=" * 70)
    print("FULL VALIDATION PIPELINE - Per-Protein Organization")
    print("=" * 70)
    print(f"\nStart: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

    # Find all protein directories
    protein_dirs = sorted([d for d in PROTEINS_DIR.iterdir() if d.is_dir()])
    print(f"\nFound {len(protein_dirs)} proteins")

    # Process each protein
    results = []
    for protein_dir in tqdm(protein_dirs, desc="Processing proteins"):
        protein_id, n_structures, status = process_protein(protein_dir)
        results.append({
            'protein': protein_id,
            'structures': n_structures,
            'status': status
        })
        tqdm.write(f"  {protein_id}: {n_structures} structures - {status}")

    # Cleanup temp directory
    shutil.rmtree(TEMP_DIR, ignore_errors=True)

    # Summary
    print(f"\n{'=' * 70}")
    print("SUMMARY")
    print("=" * 70)

    total_structures = sum(r['structures'] for r in results)
    successful = sum(1 for r in results if r['status'] == 'Success')

    print(f"Proteins processed: {successful}/{len(protein_dirs)}")
    print(f"Total structures validated: {total_structures}")
    print(f"\nOutput structure:")
    print(f"  proteins/{{PDB}}/analysis/")
    print(f"    molprobity_results.csv")
    print(f"    posebusters_results.csv")
    print(f"    VALIDATION_SUMMARY.md")

    print(f"\nEnd: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")


if __name__ == "__main__":
    main()
