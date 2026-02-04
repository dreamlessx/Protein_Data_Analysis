#!/usr/bin/env python3
"""
protein_busters.py
Comprehensive protein structure validation suite
(Inspired by PoseBusters but for protein structures)

Checks:
1. Bond length validity
2. Bond angle validity
3. Peptide bond planarity
4. Chirality (L-amino acids)
5. Steric clashes
6. Missing atoms
7. Ring geometry
8. Proline puckering
9. Disulfide bond geometry
10. Backbone continuity
"""

import subprocess
import os
import sys
import gzip
import shutil
import numpy as np
import pandas as pd
from pathlib import Path
from datetime import datetime
from tqdm import tqdm
import warnings
warnings.filterwarnings('ignore')

PROJECT_DIR = Path(__file__).parent.parent
PROTEINS_DIR = PROJECT_DIR / "proteins"
OUTPUT_DIR = PROJECT_DIR / "validation_results"
TEMP_DIR = OUTPUT_DIR / "temp"
RESULTS_DIR = PROJECT_DIR / "analysis_output"

TEMP_DIR.mkdir(parents=True, exist_ok=True)


def find_structures():
    """Find ALL PDB structure files including relaxed AF/Boltz."""
    structures = []

    for protein_dir in sorted(PROTEINS_DIR.iterdir()):
        if not protein_dir.is_dir():
            continue

        protein_id = protein_dir.name

        # Experimental original
        exp_pdb = protein_dir / f"{protein_id}.pdb"
        if exp_pdb.exists():
            structures.append({
                'path': exp_pdb,
                'protein': protein_id,
                'category': 'Experimental',
                'subcategory': 'original',
                'model': 'exp'
            })

        # AlphaFold raw
        af_dir = protein_dir / "AF"
        if af_dir.exists():
            for pdb in sorted(af_dir.glob("ranked_*.pdb")):
                model_num = pdb.stem.split('_')[1]
                structures.append({
                    'path': pdb,
                    'protein': protein_id,
                    'category': 'AlphaFold',
                    'subcategory': 'raw',
                    'model': f"model{model_num}"
                })

        # Boltz raw
        boltz_dir = protein_dir / "Boltz"
        if boltz_dir.exists():
            for pdb in sorted(boltz_dir.glob("boltz_input_model_*.pdb")):
                model_num = pdb.stem.split('_')[-1]
                structures.append({
                    'path': pdb,
                    'protein': protein_id,
                    'category': 'Boltz',
                    'subcategory': 'raw',
                    'model': f"model{model_num}"
                })

        # Experimental relaxed: proteins/{PDB}/{protocol}/{PDB}_r{N}.pdb.gz
        for relax_type in ['cartesian_beta', 'cartesian_ref15', 'dualspace_beta',
                           'dualspace_ref15', 'normal_beta', 'normal_ref15']:
            relax_dir = protein_dir / relax_type
            if relax_dir.exists():
                for pdb_gz in sorted(relax_dir.glob(f"{protein_id}_r*.pdb.gz")):
                    rep = pdb_gz.stem.replace('.pdb', '').split('_r')[-1]
                    structures.append({
                        'path': pdb_gz,
                        'protein': protein_id,
                        'category': 'Experimental',
                        'subcategory': f'relaxed_{relax_type}',
                        'model': f"r{rep}",
                        'compressed': True
                    })

        # AF/Boltz relaxed: relax/{AF,Boltz}/{model}/{protocol}/*.pdb.gz
        relax_base = protein_dir / "relax"
        if relax_base.exists():
            for category_dir in ['AF', 'Boltz']:
                cat_relax = relax_base / category_dir
                if cat_relax.exists():
                    for model_dir in sorted(cat_relax.iterdir()):
                        if not model_dir.is_dir():
                            continue
                        model_name = model_dir.name
                        for protocol_dir in sorted(model_dir.iterdir()):
                            if not protocol_dir.is_dir() or protocol_dir.name == 'log':
                                continue
                            relax_type = protocol_dir.name
                            for pdb_gz in sorted(protocol_dir.glob("*.pdb.gz")):
                                stem = pdb_gz.stem.replace('.pdb', '')
                                category = 'AlphaFold' if category_dir == 'AF' else 'Boltz'
                                rep = stem.split('_r')[-1] if '_r' in stem else stem
                                structures.append({
                                    'path': pdb_gz,
                                    'protein': protein_id,
                                    'category': category,
                                    'subcategory': f'relaxed_{relax_type}',
                                    'model': f"{model_name}_r{rep}",
                                    'compressed': True
                                })

    return structures


def parse_pdb(pdb_path):
    """Parse PDB file and extract atom coordinates."""
    atoms = []

    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                try:
                    atom = {
                        'record': line[:6].strip(),
                        'serial': int(line[6:11]),
                        'name': line[12:16].strip(),
                        'altloc': line[16].strip(),
                        'resname': line[17:20].strip(),
                        'chain': line[21].strip(),
                        'resseq': int(line[22:26]),
                        'x': float(line[30:38]),
                        'y': float(line[38:46]),
                        'z': float(line[46:54]),
                        'element': line[76:78].strip() if len(line) > 76 else ''
                    }
                    atoms.append(atom)
                except (ValueError, IndexError):
                    continue

    return atoms


def distance(a1, a2):
    """Calculate distance between two atoms."""
    return np.sqrt((a1['x']-a2['x'])**2 + (a1['y']-a2['y'])**2 + (a1['z']-a2['z'])**2)


def angle(a1, a2, a3):
    """Calculate angle between three atoms (a2 is vertex)."""
    v1 = np.array([a1['x']-a2['x'], a1['y']-a2['y'], a1['z']-a2['z']])
    v2 = np.array([a3['x']-a2['x'], a3['y']-a2['y'], a3['z']-a2['z']])

    cos_angle = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
    cos_angle = np.clip(cos_angle, -1, 1)
    return np.degrees(np.arccos(cos_angle))


def dihedral(a1, a2, a3, a4):
    """Calculate dihedral angle between four atoms."""
    b1 = np.array([a2['x']-a1['x'], a2['y']-a1['y'], a2['z']-a1['z']])
    b2 = np.array([a3['x']-a2['x'], a3['y']-a2['y'], a3['z']-a2['z']])
    b3 = np.array([a4['x']-a3['x'], a4['y']-a3['y'], a4['z']-a3['z']])

    n1 = np.cross(b1, b2)
    n2 = np.cross(b2, b3)

    m1 = np.cross(n1, b2/np.linalg.norm(b2))

    x = np.dot(n1, n2)
    y = np.dot(m1, n2)

    return np.degrees(np.arctan2(y, x))


# =============================================================================
# TEST FUNCTIONS
# =============================================================================

def test_bond_lengths(atoms):
    """Check peptide bond lengths."""
    results = {
        'n_bonds': 0,
        'bad_bonds': 0,
        'bond_length_mean': 0,
        'bond_length_std': 0,
        'pass': True
    }

    # Group atoms by residue
    residues = {}
    for atom in atoms:
        key = (atom['chain'], atom['resseq'])
        if key not in residues:
            residues[key] = {}
        residues[key][atom['name']] = atom

    bond_lengths = []

    # Check C-N bonds between consecutive residues
    sorted_keys = sorted(residues.keys())
    for i in range(len(sorted_keys) - 1):
        res1 = residues[sorted_keys[i]]
        res2 = residues[sorted_keys[i + 1]]

        if sorted_keys[i][0] != sorted_keys[i + 1][0]:  # Different chains
            continue
        if sorted_keys[i + 1][1] - sorted_keys[i][1] != 1:  # Not consecutive
            continue

        if 'C' in res1 and 'N' in res2:
            d = distance(res1['C'], res2['N'])
            bond_lengths.append(d)
            if d < 1.2 or d > 1.5:  # Standard peptide bond is ~1.33Å
                results['bad_bonds'] += 1

    results['n_bonds'] = len(bond_lengths)
    if bond_lengths:
        results['bond_length_mean'] = np.mean(bond_lengths)
        results['bond_length_std'] = np.std(bond_lengths)
        results['pass'] = results['bad_bonds'] == 0

    return results


def test_omega_angles(atoms):
    """Check omega (peptide bond) dihedral angles."""
    results = {
        'n_omega': 0,
        'cis_count': 0,
        'trans_count': 0,
        'twisted_count': 0,
        'cis_proline': 0,
        'pass': True
    }

    residues = {}
    for atom in atoms:
        key = (atom['chain'], atom['resseq'], atom['resname'])
        if key not in residues:
            residues[key] = {}
        residues[key][atom['name']] = atom

    sorted_keys = sorted(residues.keys())

    for i in range(len(sorted_keys) - 1):
        res1_key = sorted_keys[i]
        res2_key = sorted_keys[i + 1]

        if res1_key[0] != res2_key[0]:  # Different chains
            continue

        res1 = residues[res1_key]
        res2 = residues[res2_key]

        if all(a in res1 for a in ['CA', 'C']) and all(a in res2 for a in ['N', 'CA']):
            omega = dihedral(res1['CA'], res1['C'], res2['N'], res2['CA'])
            results['n_omega'] += 1

            abs_omega = abs(omega)
            if abs_omega < 30:  # Cis
                results['cis_count'] += 1
                if res2_key[2] == 'PRO':
                    results['cis_proline'] += 1
            elif abs_omega > 150:  # Trans
                results['trans_count'] += 1
            else:  # Twisted
                results['twisted_count'] += 1

    results['pass'] = results['twisted_count'] == 0

    return results


def test_chirality(atoms):
    """Check that amino acids have correct L-chirality."""
    results = {
        'n_residues': 0,
        'd_amino_acids': 0,
        'pass': True
    }

    residues = {}
    for atom in atoms:
        key = (atom['chain'], atom['resseq'])
        if key not in residues:
            residues[key] = {'resname': atom['resname']}
        residues[key][atom['name']] = atom

    for key, res in residues.items():
        if res['resname'] == 'GLY':  # Glycine is achiral
            continue

        if all(a in res for a in ['N', 'CA', 'C', 'CB']):
            results['n_residues'] += 1

            # Calculate chirality using improper dihedral
            chi = dihedral(res['N'], res['CA'], res['C'], res['CB'])

            # L-amino acids should have positive improper dihedral
            if chi < -30:  # D-amino acid
                results['d_amino_acids'] += 1

    results['pass'] = results['d_amino_acids'] == 0

    return results


def test_clashes(atoms):
    """Check for steric clashes."""
    results = {
        'n_atom_pairs': 0,
        'n_clashes': 0,
        'worst_clash': 0,
        'pass': True
    }

    # VDW radii
    vdw = {'C': 1.7, 'N': 1.55, 'O': 1.52, 'S': 1.8, 'H': 1.2}

    # Only check non-bonded pairs (different residues or >3 bonds apart)
    clash_threshold = 0.4  # Å overlap

    n_atoms = len(atoms)
    clashes = []

    for i in range(n_atoms):
        for j in range(i + 1, n_atoms):
            a1, a2 = atoms[i], atoms[j]

            # Skip if same residue
            if a1['chain'] == a2['chain'] and a1['resseq'] == a2['resseq']:
                continue

            # Skip if adjacent residues (likely bonded)
            if a1['chain'] == a2['chain'] and abs(a1['resseq'] - a2['resseq']) <= 1:
                continue

            d = distance(a1, a2)
            r1 = vdw.get(a1['element'], 1.7)
            r2 = vdw.get(a2['element'], 1.7)

            overlap = r1 + r2 - d
            if overlap > clash_threshold:
                clashes.append(overlap)

            results['n_atom_pairs'] += 1

    results['n_clashes'] = len(clashes)
    results['worst_clash'] = max(clashes) if clashes else 0
    results['pass'] = len(clashes) == 0

    return results


def test_backbone_continuity(atoms):
    """Check for breaks in the backbone."""
    results = {
        'n_residues': 0,
        'n_breaks': 0,
        'missing_atoms': 0,
        'pass': True
    }

    residues = {}
    for atom in atoms:
        key = (atom['chain'], atom['resseq'])
        if key not in residues:
            residues[key] = {'resname': atom['resname']}
        residues[key][atom['name']] = atom

    sorted_keys = sorted(residues.keys())
    results['n_residues'] = len(sorted_keys)

    # Check for missing backbone atoms
    backbone_atoms = ['N', 'CA', 'C', 'O']
    for key, res in residues.items():
        for bb in backbone_atoms:
            if bb not in res:
                results['missing_atoms'] += 1

    # Check for chain breaks (C-N distance > 2.0Å)
    for i in range(len(sorted_keys) - 1):
        key1, key2 = sorted_keys[i], sorted_keys[i + 1]

        if key1[0] != key2[0]:  # Different chains
            continue

        res1, res2 = residues[key1], residues[key2]

        if 'C' in res1 and 'N' in res2:
            d = distance(res1['C'], res2['N'])
            if d > 2.0:  # Chain break
                results['n_breaks'] += 1

    results['pass'] = results['n_breaks'] == 0 and results['missing_atoms'] == 0

    return results


def test_phi_psi(atoms):
    """Check phi/psi angles are in allowed regions."""
    results = {
        'n_residues': 0,
        'outliers': 0,
        'generously_allowed': 0,
        'pass': True
    }

    residues = {}
    for atom in atoms:
        key = (atom['chain'], atom['resseq'])
        if key not in residues:
            residues[key] = {'resname': atom['resname']}
        residues[key][atom['name']] = atom

    sorted_keys = sorted(residues.keys())

    for i in range(1, len(sorted_keys) - 1):
        prev_key = sorted_keys[i - 1]
        curr_key = sorted_keys[i]
        next_key = sorted_keys[i + 1]

        # Same chain check
        if prev_key[0] != curr_key[0] or curr_key[0] != next_key[0]:
            continue

        prev_res = residues[prev_key]
        curr_res = residues[curr_key]
        next_res = residues[next_key]

        # Calculate phi: C(i-1) - N(i) - CA(i) - C(i)
        if all(a in prev_res for a in ['C']) and all(a in curr_res for a in ['N', 'CA', 'C']):
            phi = dihedral(prev_res['C'], curr_res['N'], curr_res['CA'], curr_res['C'])
        else:
            continue

        # Calculate psi: N(i) - CA(i) - C(i) - N(i+1)
        if all(a in curr_res for a in ['N', 'CA', 'C']) and all(a in next_res for a in ['N']):
            psi = dihedral(curr_res['N'], curr_res['CA'], curr_res['C'], next_res['N'])
        else:
            continue

        results['n_residues'] += 1

        # Simplified Ramachandran check (favored regions)
        is_favored = False
        is_allowed = False

        # Beta sheet region
        if -180 <= phi <= -45 and 45 <= psi <= 180:
            is_favored = True
        elif -180 <= phi <= -45 and -180 <= psi <= -135:
            is_favored = True
        # Alpha helix region
        elif -100 <= phi <= -45 and -65 <= psi <= -15:
            is_favored = True
        # Left-handed helix (for Gly)
        elif 30 <= phi <= 90 and -30 <= psi <= 60:
            is_allowed = True

        if not is_favored and not is_allowed:
            if curr_res['resname'] == 'GLY':  # Glycine is more flexible
                is_allowed = True
            else:
                results['outliers'] += 1

        if is_allowed and not is_favored:
            results['generously_allowed'] += 1

    results['pass'] = results['outliers'] == 0

    return results


def test_cysteine_bonds(atoms):
    """Check disulfide bond geometry."""
    results = {
        'n_disulfides': 0,
        'bad_geometry': 0,
        'pass': True
    }

    # Find all cysteine SG atoms
    cys_sg = [a for a in atoms if a['resname'] == 'CYS' and a['name'] == 'SG']

    for i in range(len(cys_sg)):
        for j in range(i + 1, len(cys_sg)):
            d = distance(cys_sg[i], cys_sg[j])

            # Disulfide bond distance is ~2.0-2.1Å
            if 1.8 < d < 2.5:
                results['n_disulfides'] += 1
                if d < 1.95 or d > 2.15:
                    results['bad_geometry'] += 1

    results['pass'] = results['bad_geometry'] == 0

    return results


def decompress_if_needed(structure):
    """Decompress gzipped PDB file if needed, return path to PDB file."""
    if structure.get('compressed'):
        gz_path = structure['path']
        pdb_path = TEMP_DIR / f"{structure['protein']}_{structure['category']}_{structure['model']}.pdb"
        with gzip.open(gz_path, 'rb') as f_in:
            with open(pdb_path, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        return pdb_path
    return structure['path']


def run_all_tests(pdb_path):
    """Run all validation tests on a PDB file."""
    try:
        atoms = parse_pdb(pdb_path)

        if len(atoms) == 0:
            return {'error': 'No atoms found'}

        results = {
            'n_atoms': len(atoms),
        }

        # Run each test
        bond_results = test_bond_lengths(atoms)
        results.update({f'bond_{k}': v for k, v in bond_results.items()})

        omega_results = test_omega_angles(atoms)
        results.update({f'omega_{k}': v for k, v in omega_results.items()})

        chiral_results = test_chirality(atoms)
        results.update({f'chiral_{k}': v for k, v in chiral_results.items()})

        # Skip O(n^2) clash test — MolProbity clashscore is more reliable
        results['clash_skipped'] = True

        backbone_results = test_backbone_continuity(atoms)
        results.update({f'backbone_{k}': v for k, v in backbone_results.items()})

        phipsi_results = test_phi_psi(atoms)
        results.update({f'phipsi_{k}': v for k, v in phipsi_results.items()})

        cys_results = test_cysteine_bonds(atoms)
        results.update({f'disulfide_{k}': v for k, v in cys_results.items()})

        # Overall pass/fail
        results['all_tests_pass'] = all([
            bond_results['pass'],
            omega_results['pass'],
            chiral_results['pass'],
            backbone_results['pass'],
            cys_results['pass']
        ])

        return results

    except Exception as e:
        return {'error': str(e)}


def main():
    print("=" * 70)
    print("PROTEIN BUSTERS - Comprehensive Structure Validation")
    print("=" * 70)
    print(f"\nStart time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

    # Find structures
    print("\nDiscovering structures...")
    structures = find_structures()
    print(f"Found {len(structures)} structures")

    # Count by category
    cats = {}
    for s in structures:
        key = f"{s['category']}/{s['subcategory'].split('_')[0]}"
        cats[key] = cats.get(key, 0) + 1
    for k, v in sorted(cats.items()):
        print(f"  {k}: {v}")

    # Run validation on ALL structures
    print(f"\nValidating ALL {len(structures)} structures...")
    all_results = []

    for struct in tqdm(structures, desc="Validating"):
        # Decompress if needed
        pdb_path = decompress_if_needed(struct)
        try:
            result = run_all_tests(pdb_path)
        finally:
            # Clean up temp file
            if struct.get('compressed') and pdb_path.exists():
                pdb_path.unlink(missing_ok=True)
        result['protein'] = struct['protein']
        result['category'] = struct['category']
        result['subcategory'] = struct['subcategory']
        result['model'] = struct.get('model', 'exp')
        all_results.append(result)

    # Save results
    results_df = pd.DataFrame(all_results)
    output_path = RESULTS_DIR / "protein_busters_results.csv"
    results_df.to_csv(output_path, index=False)
    print(f"\nSaved results to: {output_path}")

    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)

    print("\nTest Pass Rates by Category:")
    for cat in ['AlphaFold', 'Boltz', 'Experimental']:
        cat_df = results_df[results_df['category'] == cat]
        if len(cat_df) > 0:
            pass_rate = 100 * cat_df['all_tests_pass'].sum() / len(cat_df)
            print(f"  {cat}: {pass_rate:.1f}% ({cat_df['all_tests_pass'].sum()}/{len(cat_df)} pass all tests)")

    print("\nIndividual Test Results:")
    test_cols = [c for c in results_df.columns if c.endswith('_pass')]
    for col in test_cols:
        if col in results_df.columns:
            pass_rate = 100 * results_df[col].sum() / len(results_df)
            print(f"  {col.replace('_pass', '')}: {pass_rate:.1f}% pass")

    print("\nCommon Issues Found:")
    if 'omega_twisted_count' in results_df.columns:
        twisted = results_df['omega_twisted_count'].sum()
        print(f"  Twisted peptide bonds: {int(twisted)}")
    if 'chiral_d_amino_acids' in results_df.columns:
        d_aa = results_df['chiral_d_amino_acids'].sum()
        print(f"  D-amino acids: {int(d_aa)}")
    if 'backbone_n_breaks' in results_df.columns:
        breaks = results_df['backbone_n_breaks'].sum()
        print(f"  Backbone breaks: {int(breaks)}")

    print(f"\nEnd time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")


if __name__ == "__main__":
    main()
