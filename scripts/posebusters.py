#!/usr/bin/env python3
"""
posebusters.py
Protein structure validity checks inspired by PoseBusters (adapted for proteins).

Tests:
1. structure_loaded - Can the structure be parsed?
2. valid_residues - Are all residues recognized amino acids?
3. backbone_connected - Is the backbone continuous?
4. bond_lengths - Are peptide bond lengths within tolerance?
5. bond_angles - Are backbone bond angles within tolerance?
6. steric_clashes - Are there severe steric clashes?
7. aromatic_flatness - Are aromatic rings planar?
8. peptide_planarity - Are peptide bonds planar (omega angle)?
9. chirality - Are amino acids L-configuration?
10. complete_residues - Are all backbone atoms present?
11. internal_energy - Rosetta energy score (optional)
"""

import argparse
import gzip
import shutil
import subprocess
import tempfile
import os
import sys
import numpy as np
import pandas as pd
from pathlib import Path
from datetime import datetime
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm
import warnings
warnings.filterwarnings('ignore')

PROJECT_DIR = Path(__file__).parent.parent
PROTEINS_DIR = PROJECT_DIR / "proteins"
OUTPUT_DIR = PROJECT_DIR / "validation_results"
PER_PROTEIN_DIR = OUTPUT_DIR / "per_protein"
TEMP_DIR = OUTPUT_DIR / "temp"

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
PER_PROTEIN_DIR.mkdir(parents=True, exist_ok=True)
TEMP_DIR.mkdir(parents=True, exist_ok=True)

STANDARD_AA = {
    'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
    'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'
}

AROMATIC_RESIDUES = {'PHE', 'TYR', 'TRP', 'HIS'}

ROSETTA_BINS = [
    '/Users/muditagar/Rosetta/main/source/bin/score_jd2.default.macosclangrelease',
    'score_jd2.default.macosclangrelease',
    'score_jd2.default.linuxgccrelease',
    'score_jd2',
]


def find_rosetta():
    """Find Rosetta score_jd2 binary."""
    for bin_path in ROSETTA_BINS:
        if Path(bin_path).exists():
            return bin_path
        if shutil.which(bin_path):
            return bin_path
    return None


def find_structures():
    """Find ALL PDB structure files."""
    structures = []

    for protein_dir in sorted(PROTEINS_DIR.iterdir()):
        if not protein_dir.is_dir():
            continue

        protein_id = protein_dir.name

        # Experimental original
        exp_pdb = protein_dir / f"{protein_id}.pdb"
        if exp_pdb.exists():
            structures.append({
                'path': str(exp_pdb),
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
                    'path': str(pdb),
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
                    'path': str(pdb),
                    'protein': protein_id,
                    'category': 'Boltz',
                    'subcategory': 'raw',
                    'model': f"model{model_num}"
                })

        # Experimental relaxed
        for relax_type in ['cartesian_beta', 'cartesian_ref15', 'dualspace_beta',
                           'dualspace_ref15', 'normal_beta', 'normal_ref15']:
            relax_dir = protein_dir / relax_type
            if relax_dir.exists():
                for pdb_gz in sorted(relax_dir.glob(f"{protein_id}_r*.pdb.gz")):
                    rep = pdb_gz.stem.replace('.pdb', '').split('_r')[-1]
                    structures.append({
                        'path': str(pdb_gz),
                        'protein': protein_id,
                        'category': 'Experimental',
                        'subcategory': f'relaxed_{relax_type}',
                        'model': f"r{rep}",
                        'compressed': True
                    })

        # AF/Boltz relaxed
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
                                    'path': str(pdb_gz),
                                    'protein': protein_id,
                                    'category': category,
                                    'subcategory': f'relaxed_{relax_type}',
                                    'model': f"{model_name}_r{rep}",
                                    'compressed': True
                                })

    return structures


def parse_pdb(pdb_path: str) -> list:
    """Parse PDB file and extract atom coordinates."""
    atoms = []
    try:
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
                        if not atom['element']:
                            atom['element'] = atom['name'][0]
                        atoms.append(atom)
                    except (ValueError, IndexError):
                        continue
    except Exception:
        return []
    return atoms


def distance(a1, a2) -> float:
    """Calculate distance between two atoms."""
    return np.sqrt((a1['x']-a2['x'])**2 + (a1['y']-a2['y'])**2 + (a1['z']-a2['z'])**2)


def angle(a1, a2, a3) -> float:
    """Calculate angle between three atoms (a2 is vertex)."""
    v1 = np.array([a1['x']-a2['x'], a1['y']-a2['y'], a1['z']-a2['z']])
    v2 = np.array([a3['x']-a2['x'], a3['y']-a2['y'], a3['z']-a2['z']])
    n1, n2 = np.linalg.norm(v1), np.linalg.norm(v2)
    if n1 == 0 or n2 == 0:
        return 0.0
    cos_angle = np.dot(v1, v2) / (n1 * n2)
    cos_angle = np.clip(cos_angle, -1, 1)
    return np.degrees(np.arccos(cos_angle))


def dihedral(a1, a2, a3, a4) -> float:
    """Calculate dihedral angle between four atoms."""
    b1 = np.array([a2['x']-a1['x'], a2['y']-a1['y'], a2['z']-a1['z']])
    b2 = np.array([a3['x']-a2['x'], a3['y']-a2['y'], a3['z']-a2['z']])
    b3 = np.array([a4['x']-a3['x'], a4['y']-a3['y'], a4['z']-a3['z']])

    n1 = np.cross(b1, b2)
    n2 = np.cross(b2, b3)

    n1_norm = np.linalg.norm(n1)
    n2_norm = np.linalg.norm(n2)
    b2_norm = np.linalg.norm(b2)

    if n1_norm == 0 or n2_norm == 0 or b2_norm == 0:
        return 0.0

    m1 = np.cross(n1, b2 / b2_norm)
    x = np.dot(n1, n2)
    y = np.dot(m1, n2)

    return np.degrees(np.arctan2(y, x))


def group_by_residue(atoms: list) -> dict:
    """Group atoms by (chain, resseq, resname)."""
    residues = {}
    for atom in atoms:
        key = (atom['chain'], atom['resseq'], atom['resname'])
        if key not in residues:
            residues[key] = {}
        residues[key][atom['name']] = atom
    return residues


# =============================================================================
# TEST FUNCTIONS - Return both pass/fail AND raw values
# =============================================================================

def test_structure_loaded(atoms: list) -> dict:
    """Test 1: Can the structure be parsed?"""
    n_atoms = len(atoms)
    return {
        'structure_loaded': n_atoms > 0,
        'raw_n_atoms': n_atoms
    }


def test_valid_residues(atoms: list) -> dict:
    """Test 2: Are all residues recognized amino acids?"""
    residues = set(a['resname'] for a in atoms)
    non_standard = residues - STANDARD_AA
    n_total = len(residues)
    n_valid = len(residues & STANDARD_AA)
    return {
        'valid_residues': len(non_standard) == 0,
        'raw_n_residue_types': n_total,
        'raw_n_valid_residue_types': n_valid,
        'raw_non_standard_residues': ','.join(sorted(non_standard)) if non_standard else ''
    }


def test_backbone_connected(atoms: list) -> dict:
    """Test 3: Is the backbone continuous?"""
    residues = {}
    for atom in atoms:
        key = (atom['chain'], atom['resseq'])
        if key not in residues:
            residues[key] = {}
        residues[key][atom['name']] = atom

    sorted_keys = sorted(residues.keys())
    n_breaks = 0
    max_break_dist = 0.0

    for i in range(len(sorted_keys) - 1):
        key1, key2 = sorted_keys[i], sorted_keys[i + 1]
        if key1[0] != key2[0]:  # Different chains
            continue

        res1, res2 = residues[key1], residues[key2]
        if 'C' in res1 and 'N' in res2:
            d = distance(res1['C'], res2['N'])
            if d > 2.0:  # Chain break threshold
                n_breaks += 1
                max_break_dist = max(max_break_dist, d)

    return {
        'backbone_connected': n_breaks == 0,
        'raw_n_backbone_breaks': n_breaks,
        'raw_max_break_distance': round(max_break_dist, 3)
    }


def test_bond_lengths(atoms: list) -> dict:
    """Test 4: Are peptide bond lengths within tolerance?"""
    residues = {}
    for atom in atoms:
        key = (atom['chain'], atom['resseq'])
        if key not in residues:
            residues[key] = {}
        residues[key][atom['name']] = atom

    bond_lengths = []
    n_outliers = 0
    sorted_keys = sorted(residues.keys())

    for i in range(len(sorted_keys) - 1):
        key1, key2 = sorted_keys[i], sorted_keys[i + 1]
        if key1[0] != key2[0]:
            continue
        if key2[1] - key1[1] != 1:
            continue

        res1, res2 = residues[key1], residues[key2]
        if 'C' in res1 and 'N' in res2:
            d = distance(res1['C'], res2['N'])
            bond_lengths.append(d)
            # Ideal C-N peptide bond: 1.33 Å, tolerance ±0.15 Å
            if d < 1.18 or d > 1.48:
                n_outliers += 1

    mean_len = np.mean(bond_lengths) if bond_lengths else 0
    std_len = np.std(bond_lengths) if bond_lengths else 0

    return {
        'bond_lengths': n_outliers == 0,
        'raw_n_peptide_bonds': len(bond_lengths),
        'raw_n_bond_outliers': n_outliers,
        'raw_mean_bond_length': round(mean_len, 4),
        'raw_std_bond_length': round(std_len, 4)
    }


def test_bond_angles(atoms: list) -> dict:
    """Test 5: Are backbone bond angles within tolerance?"""
    residues = group_by_residue(atoms)
    sorted_keys = sorted(residues.keys())

    angles_list = []
    n_outliers = 0

    for key in sorted_keys:
        res = residues[key]
        # N-CA-C angle (ideal ~111°)
        if all(a in res for a in ['N', 'CA', 'C']):
            ang = angle(res['N'], res['CA'], res['C'])
            angles_list.append(ang)
            if ang < 100 or ang > 120:
                n_outliers += 1

    mean_ang = np.mean(angles_list) if angles_list else 0
    std_ang = np.std(angles_list) if angles_list else 0

    return {
        'bond_angles': n_outliers == 0,
        'raw_n_backbone_angles': len(angles_list),
        'raw_n_angle_outliers': n_outliers,
        'raw_mean_backbone_angle': round(mean_ang, 2),
        'raw_std_backbone_angle': round(std_ang, 2)
    }


def test_steric_clashes(atoms: list) -> dict:
    """Test 6: Are there severe steric clashes?"""
    vdw = {'C': 1.7, 'N': 1.55, 'O': 1.52, 'S': 1.8, 'H': 1.2}
    clash_threshold = 0.5  # Å overlap for severe clash

    n_atoms = len(atoms)
    n_clashes = 0
    worst_clash = 0.0

    # Sample atoms if too many (O(n^2) is slow)
    if n_atoms > 1000:
        indices = np.random.choice(n_atoms, 1000, replace=False)
        sample_atoms = [atoms[i] for i in indices]
    else:
        sample_atoms = atoms

    for i in range(len(sample_atoms)):
        for j in range(i + 1, len(sample_atoms)):
            a1, a2 = sample_atoms[i], sample_atoms[j]

            # Skip same residue
            if a1['chain'] == a2['chain'] and a1['resseq'] == a2['resseq']:
                continue
            # Skip adjacent residues
            if a1['chain'] == a2['chain'] and abs(a1['resseq'] - a2['resseq']) <= 1:
                continue

            d = distance(a1, a2)
            r1 = vdw.get(a1['element'], 1.7)
            r2 = vdw.get(a2['element'], 1.7)
            overlap = r1 + r2 - d

            if overlap > clash_threshold:
                n_clashes += 1
                worst_clash = max(worst_clash, overlap)

    return {
        'steric_clashes': n_clashes < 5,  # Allow a few minor clashes
        'raw_n_clashes': n_clashes,
        'raw_worst_clash': round(worst_clash, 3)
    }


def test_aromatic_flatness(atoms: list) -> dict:
    """Test 7: Are aromatic rings planar?"""
    residues = group_by_residue(atoms)

    n_aromatic = 0
    n_nonplanar = 0
    max_deviation = 0.0

    for key, res in residues.items():
        if key[2] not in AROMATIC_RESIDUES:
            continue

        # Get ring atoms based on residue type
        if key[2] == 'PHE':
            ring_atoms = ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ']
        elif key[2] == 'TYR':
            ring_atoms = ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ']
        elif key[2] == 'TRP':
            ring_atoms = ['CG', 'CD1', 'CD2', 'NE1', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2']
        elif key[2] == 'HIS':
            ring_atoms = ['CG', 'ND1', 'CD2', 'CE1', 'NE2']
        else:
            continue

        coords = []
        for atom_name in ring_atoms:
            if atom_name in res:
                a = res[atom_name]
                coords.append([a['x'], a['y'], a['z']])

        if len(coords) < 4:
            continue

        n_aromatic += 1
        coords = np.array(coords)

        # Fit plane and calculate RMSD
        centroid = coords.mean(axis=0)
        centered = coords - centroid
        _, _, vh = np.linalg.svd(centered)
        normal = vh[2]

        distances = np.abs(np.dot(centered, normal))
        rmsd = np.sqrt(np.mean(distances**2))
        max_deviation = max(max_deviation, rmsd)

        if rmsd > 0.1:  # 0.1 Å threshold
            n_nonplanar += 1

    return {
        'aromatic_flatness': n_nonplanar == 0,
        'raw_n_aromatic_rings': n_aromatic,
        'raw_n_nonplanar_rings': n_nonplanar,
        'raw_max_ring_deviation': round(max_deviation, 4)
    }


def test_peptide_planarity(atoms: list) -> dict:
    """Test 8: Are peptide bonds planar (omega angle)?"""
    residues = group_by_residue(atoms)
    sorted_keys = sorted(residues.keys())

    n_omega = 0
    n_cis = 0
    n_trans = 0
    n_twisted = 0
    omega_values = []

    for i in range(len(sorted_keys) - 1):
        key1, key2 = sorted_keys[i], sorted_keys[i + 1]
        if key1[0] != key2[0]:
            continue

        res1, res2 = residues[key1], residues[key2]

        if all(a in res1 for a in ['CA', 'C']) and all(a in res2 for a in ['N', 'CA']):
            omega = dihedral(res1['CA'], res1['C'], res2['N'], res2['CA'])
            omega_values.append(omega)
            n_omega += 1

            abs_omega = abs(omega)
            if abs_omega < 30:
                n_cis += 1
            elif abs_omega > 150:
                n_trans += 1
            else:
                n_twisted += 1

    return {
        'peptide_planarity': n_twisted == 0,
        'raw_n_omega_angles': n_omega,
        'raw_n_cis': n_cis,
        'raw_n_trans': n_trans,
        'raw_n_twisted': n_twisted,
        'raw_mean_abs_omega': round(np.mean(np.abs(omega_values)), 2) if omega_values else 0
    }


def test_chirality(atoms: list) -> dict:
    """Test 9: Are amino acids L-configuration?"""
    residues = group_by_residue(atoms)

    n_checked = 0
    n_d_amino = 0

    for key, res in residues.items():
        if key[2] == 'GLY':  # Glycine is achiral
            continue

        if all(a in res for a in ['N', 'CA', 'C', 'CB']):
            n_checked += 1
            chi = dihedral(res['N'], res['CA'], res['C'], res['CB'])
            if chi < -30:  # D-amino acid
                n_d_amino += 1

    return {
        'chirality': n_d_amino == 0,
        'raw_n_chiral_centers': n_checked,
        'raw_n_d_amino_acids': n_d_amino
    }


def test_complete_residues(atoms: list) -> dict:
    """Test 10: Are all backbone atoms present?"""
    residues = group_by_residue(atoms)
    backbone = ['N', 'CA', 'C', 'O']

    n_residues = len(residues)
    n_incomplete = 0
    missing_atoms = 0

    for key, res in residues.items():
        missing = [a for a in backbone if a not in res]
        if missing:
            n_incomplete += 1
            missing_atoms += len(missing)

    return {
        'complete_residues': n_incomplete == 0,
        'raw_n_residues': n_residues,
        'raw_n_incomplete_residues': n_incomplete,
        'raw_n_missing_backbone_atoms': missing_atoms
    }


def test_internal_energy(pdb_path: str, rosetta_bin: str) -> dict:
    """Test 11: Rosetta internal energy score."""
    if not rosetta_bin:
        return {
            'internal_energy': None,
            'raw_rosetta_score': None
        }

    try:
        with tempfile.TemporaryDirectory() as tmpdir:
            score_file = Path(tmpdir) / 'score.sc'
            cmd = [
                rosetta_bin,
                '-in:file:s', str(pdb_path),
                '-out:file:scorefile', str(score_file),
                '-ignore_unrecognized_res',
                '-mute', 'all'
            ]
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)

            if result.returncode != 0:
                return {
                    'internal_energy': None,
                    'raw_rosetta_score': None
                }

            if not score_file.exists():
                return {
                    'internal_energy': None,
                    'raw_rosetta_score': None
                }

            content = score_file.read_text()
            for line in content.split('\n'):
                if line.startswith('SCORE:') and 'total_score' not in line:
                    parts = line.split()
                    if len(parts) > 1:
                        try:
                            score = float(parts[1])
                            return {
                                'internal_energy': score < 0,
                                'raw_rosetta_score': round(score, 2)
                            }
                        except ValueError:
                            pass

            return {
                'internal_energy': None,
                'raw_rosetta_score': None
            }
    except subprocess.TimeoutExpired:
        return {
            'internal_energy': None,
            'raw_rosetta_score': None
        }
    except Exception:
        return {
            'internal_energy': None,
            'raw_rosetta_score': None
        }


def decompress_pdb(gz_path: str) -> str:
    """Decompress gzipped PDB file to unique temp location."""
    import uuid
    unique_id = uuid.uuid4().hex[:8]
    pdb_name = Path(gz_path).stem + f'_{unique_id}'
    pdb_path = TEMP_DIR / pdb_name
    with gzip.open(gz_path, 'rb') as f_in:
        with open(pdb_path, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    return str(pdb_path)


def validate_structure(args) -> dict:
    """Validate a single structure. Returns dict with all results."""
    struct, rosetta_bin = args

    result = {
        'protein': struct['protein'],
        'category': struct['category'],
        'subcategory': struct['subcategory'],
        'model': struct.get('model', 'exp'),
    }

    # Get PDB path (decompress if needed)
    pdb_path = struct['path']
    temp_pdb = None

    try:
        if struct.get('compressed'):
            pdb_path = decompress_pdb(pdb_path)
            temp_pdb = pdb_path

        atoms = parse_pdb(pdb_path)

        # Run all tests
        tests = [
            test_structure_loaded(atoms),
            test_valid_residues(atoms),
            test_backbone_connected(atoms),
            test_bond_lengths(atoms),
            test_bond_angles(atoms),
            test_steric_clashes(atoms),
            test_aromatic_flatness(atoms),
            test_peptide_planarity(atoms),
            test_chirality(atoms),
            test_complete_residues(atoms),
        ]

        for t in tests:
            result.update(t)

        # Rosetta energy (optional)
        if rosetta_bin:
            energy_result = test_internal_energy(pdb_path, rosetta_bin)
            result.update(energy_result)
        else:
            result['internal_energy'] = None
            result['raw_rosetta_score'] = None

        # Calculate all_pass and n_pass
        pass_cols = [
            'structure_loaded', 'valid_residues', 'backbone_connected',
            'bond_lengths', 'bond_angles', 'steric_clashes',
            'aromatic_flatness', 'peptide_planarity', 'chirality', 'complete_residues'
        ]
        n_pass = sum(1 for c in pass_cols if result.get(c) is True)
        all_pass = all(result.get(c) is True for c in pass_cols)

        if result.get('internal_energy') is not None:
            if result['internal_energy']:
                n_pass += 1
            else:
                all_pass = False

        result['all_pass'] = all_pass
        result['n_pass'] = n_pass

    except Exception as e:
        result['error'] = str(e)
        result['all_pass'] = False
        result['n_pass'] = 0

    finally:
        # Clean up temp file
        if temp_pdb and Path(temp_pdb).exists():
            try:
                Path(temp_pdb).unlink()
            except Exception:
                pass

    return result


def get_pass_cols():
    return ['protein', 'category', 'subcategory', 'model',
            'structure_loaded', 'valid_residues', 'backbone_connected',
            'bond_lengths', 'bond_angles', 'steric_clashes',
            'aromatic_flatness', 'peptide_planarity', 'chirality',
            'complete_residues', 'internal_energy', 'all_pass', 'n_pass']


def get_raw_cols(df):
    return ['protein', 'category', 'subcategory', 'model'] + \
           [c for c in df.columns if c.startswith('raw_')]


def save_results(results: list, output_name: str, raw_name: str):
    """Save pass/fail results and raw values to separate files."""
    df = pd.DataFrame(results)

    pass_cols = [c for c in get_pass_cols() if c in df.columns]
    pass_df = df[pass_cols]
    pass_df.to_csv(OUTPUT_DIR / output_name, index=False)

    raw_cols = [c for c in get_raw_cols(df) if c in df.columns]
    raw_df = df[raw_cols]
    raw_df.to_csv(OUTPUT_DIR / raw_name, index=False)

    return pass_df, raw_df


def sort_results(df: pd.DataFrame) -> pd.DataFrame:
    """Sort results in logical order: Experimental > AlphaFold > Boltz, original/raw > relaxed."""
    cat_order = {'Experimental': 0, 'AlphaFold': 1, 'Boltz': 2}
    sub_order = {
        'original': 0, 'raw': 1,
        'relaxed_cartesian_beta': 2, 'relaxed_cartesian_ref15': 3,
        'relaxed_dualspace_beta': 4, 'relaxed_dualspace_ref15': 5,
        'relaxed_normal_beta': 6, 'relaxed_normal_ref15': 7
    }
    df = df.copy()
    df['_cat_order'] = df['category'].map(cat_order).fillna(99)
    df['_sub_order'] = df['subcategory'].map(sub_order).fillna(99)
    df = df.sort_values(['_cat_order', '_sub_order', 'model'])
    df = df.drop(columns=['_cat_order', '_sub_order'])
    return df.reset_index(drop=True)


def save_per_protein_incremental(results: list, protein: str):
    """Save results for a single protein to proteins/{protein}/analysis/."""
    df = pd.DataFrame(results)
    df = sort_results(df)

    # Create analysis directory
    analysis_dir = PROTEINS_DIR / protein / "analysis"
    analysis_dir.mkdir(parents=True, exist_ok=True)

    # Pass/fail results
    pass_cols = ['category', 'subcategory', 'model',
                 'structure_loaded', 'valid_residues', 'backbone_connected',
                 'bond_lengths', 'bond_angles', 'steric_clashes',
                 'aromatic_flatness', 'peptide_planarity', 'chirality',
                 'complete_residues', 'internal_energy', 'all_pass', 'n_pass']
    pass_cols = [c for c in pass_cols if c in df.columns]
    pass_df = df[pass_cols]
    pass_df.to_csv(analysis_dir / "posebusters_results.csv", index=False)

    # Raw values
    raw_cols = ['category', 'subcategory', 'model'] + \
               [c for c in df.columns if c.startswith('raw_')]
    raw_cols = [c for c in raw_cols if c in df.columns]
    raw_df = df[raw_cols]
    raw_df.to_csv(analysis_dir / "posebusters_raw.csv", index=False)

    return len(results)


def main():
    parser = argparse.ArgumentParser(description='PoseBusters for proteins')
    parser.add_argument('--no-energy', action='store_true', help='Skip Rosetta energy calculation')
    parser.add_argument('--limit', type=int, help='Limit number of structures')
    parser.add_argument('-j', '--workers', type=int, default=os.cpu_count(), help='Number of parallel workers')
    parser.add_argument('--suffix', type=str, default='', help='Suffix for output files')
    args = parser.parse_args()

    print("=" * 70)
    print("POSEBUSTERS - Protein Structure Validity Checks")
    print("=" * 70)
    print(f"\nStart: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Workers: {args.workers}")

    # Find Rosetta
    rosetta_bin = None if args.no_energy else find_rosetta()
    if rosetta_bin:
        print(f"Rosetta: {rosetta_bin}")
    else:
        print("(Rosetta energy calculation disabled)")

    # Find structures
    structures = find_structures()
    print(f"Found {len(structures)} structures")

    if args.limit:
        structures = structures[:args.limit]
        print(f"Limited to {len(structures)} structures")

    # Group structures by protein for incremental saving
    by_protein = {}
    for s in structures:
        prot = s['protein']
        if prot not in by_protein:
            by_protein[prot] = []
        by_protein[prot].append(s)

    proteins = sorted(by_protein.keys())
    print(f"Processing {len(proteins)} proteins incrementally")

    all_results = []
    completed_proteins = 0

    # Process protein by protein
    for protein in proteins:
        protein_structures = by_protein[protein]
        tasks = [(s, rosetta_bin) for s in protein_structures]

        protein_results = []
        with ProcessPoolExecutor(max_workers=args.workers) as executor:
            futures = {executor.submit(validate_structure, t): t[0] for t in tasks}
            for future in tqdm(as_completed(futures), total=len(futures),
                              desc=f"{protein} ({completed_proteins+1}/{len(proteins)})",
                              leave=False):
                result = future.result()
                protein_results.append(result)
                all_results.append(result)

        # Save per-protein results immediately
        n_saved = save_per_protein_incremental(protein_results, protein)
        completed_proteins += 1
        print(f"[{completed_proteins}/{len(proteins)}] {protein}: {n_saved} structures saved to proteins/{protein}/analysis/", flush=True)

        # Also append to compiled results file (incremental)
        suffix_str = f'_{args.suffix}' if args.suffix else ''
        compiled_path = OUTPUT_DIR / f"posebusters_results{suffix_str}.csv"
        compiled_raw_path = OUTPUT_DIR / f"posebusters_raw{suffix_str}.csv"

        df = pd.DataFrame(protein_results)
        df = sort_results(df)
        pass_cols = [c for c in get_pass_cols() if c in df.columns]
        raw_cols = [c for c in get_raw_cols(df) if c in df.columns]

        # Append to compiled files
        write_header = not compiled_path.exists()
        df[pass_cols].to_csv(compiled_path, mode='a', header=write_header, index=False)
        df[raw_cols].to_csv(compiled_raw_path, mode='a', header=write_header, index=False)

    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)

    pass_df = pd.DataFrame(all_results)
    print(f"\nTotal structures: {len(all_results)}")
    print(f"All tests pass: {pass_df['all_pass'].sum()} ({100*pass_df['all_pass'].mean():.1f}%)")

    print("\nBy Category:")
    for cat in ['Experimental', 'AlphaFold', 'Boltz']:
        cat_df = pass_df[pass_df['category'] == cat]
        if len(cat_df) > 0:
            pct = 100 * cat_df['all_pass'].sum() / len(cat_df)
            print(f"  {cat}: {cat_df['all_pass'].sum()}/{len(cat_df)} pass ({pct:.1f}%)")

    print("\nIndividual Test Pass Rates:")
    test_cols = ['structure_loaded', 'valid_residues', 'backbone_connected',
                 'bond_lengths', 'bond_angles', 'steric_clashes',
                 'aromatic_flatness', 'peptide_planarity', 'chirality', 'complete_residues']
    for col in test_cols:
        if col in pass_df.columns:
            pct = 100 * pass_df[col].sum() / len(pass_df)
            print(f"  {col}: {pct:.1f}%")

    if 'internal_energy' in pass_df.columns:
        valid = pass_df['internal_energy'].notna()
        if valid.any():
            pct = 100 * pass_df.loc[valid, 'internal_energy'].sum() / valid.sum()
            print(f"  internal_energy: {pct:.1f}%")

    print(f"\nResults saved to:")
    print(f"  Per-protein: proteins/{{PROTEIN}}/analysis/posebusters_{{results,raw}}.csv")
    print(f"  Compiled: {OUTPUT_DIR}/posebusters_{{results,raw}}.csv")
    print(f"\nEnd: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")


if __name__ == "__main__":
    main()
