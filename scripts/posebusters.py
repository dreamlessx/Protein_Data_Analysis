#!/usr/bin/env python3
"""
Protein structure validity checks adapted from PoseBusters methodology.
Validates geometry, connectivity, and energetics of protein structures.
"""

import argparse
import gzip
import shutil
import subprocess
import tempfile
import os
import uuid
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
TEMP_DIR = OUTPUT_DIR / "temp"

for d in [OUTPUT_DIR, TEMP_DIR]:
    d.mkdir(parents=True, exist_ok=True)

STANDARD_AA = {
    'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
    'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'
}

AROMATIC_RESIDUES = {'PHE', 'TYR', 'TRP', 'HIS'}

ROSETTA_BINS = [
    'score_jd2.default.macosclangrelease',
    'score_jd2.default.linuxgccrelease',
    'score_jd2',
]

VDW_RADII = {'C': 1.7, 'N': 1.55, 'O': 1.52, 'S': 1.8, 'H': 1.2}

RING_ATOMS = {
    'PHE': ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'],
    'TYR': ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'],
    'TRP': ['CG', 'CD1', 'CD2', 'NE1', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2'],
    'HIS': ['CG', 'ND1', 'CD2', 'CE1', 'NE2']
}


def find_rosetta():
    for path in ROSETTA_BINS:
        if Path(path).exists() or shutil.which(path):
            return path
    return None


def find_structures():
    """Enumerate all PDB files across experimental, AlphaFold, and Boltz predictions."""
    structures = []

    for protein_dir in sorted(PROTEINS_DIR.iterdir()):
        if not protein_dir.is_dir():
            continue

        pid = protein_dir.name

        # Experimental
        exp_pdb = protein_dir / f"{pid}.pdb"
        if exp_pdb.exists():
            structures.append({
                'path': str(exp_pdb), 'protein': pid,
                'category': 'Experimental', 'subcategory': 'original', 'model': 'exp'
            })

        # AlphaFold raw
        af_dir = protein_dir / "AF"
        if af_dir.exists():
            for pdb in sorted(af_dir.glob("ranked_*.pdb")):
                structures.append({
                    'path': str(pdb), 'protein': pid,
                    'category': 'AlphaFold', 'subcategory': 'raw',
                    'model': f"model{pdb.stem.split('_')[1]}"
                })

        # Boltz raw
        boltz_dir = protein_dir / "Boltz"
        if boltz_dir.exists():
            for pdb in sorted(boltz_dir.glob("boltz_input_model_*.pdb")):
                structures.append({
                    'path': str(pdb), 'protein': pid,
                    'category': 'Boltz', 'subcategory': 'raw',
                    'model': f"model{pdb.stem.split('_')[-1]}"
                })

        # Relaxed structures (experimental)
        for protocol in ['cartesian_beta', 'cartesian_ref15', 'dualspace_beta',
                         'dualspace_ref15', 'normal_beta', 'normal_ref15']:
            relax_dir = protein_dir / protocol
            if relax_dir.exists():
                for pdb_gz in sorted(relax_dir.glob(f"{pid}_r*.pdb.gz")):
                    rep = pdb_gz.stem.replace('.pdb', '').split('_r')[-1]
                    structures.append({
                        'path': str(pdb_gz), 'protein': pid,
                        'category': 'Experimental', 'subcategory': f'relaxed_{protocol}',
                        'model': f"r{rep}", 'compressed': True
                    })

        # Relaxed structures (AF/Boltz)
        relax_base = protein_dir / "relax"
        if relax_base.exists():
            for cat_name in ['AF', 'Boltz']:
                cat_dir = relax_base / cat_name
                if not cat_dir.exists():
                    continue
                for model_dir in sorted(cat_dir.iterdir()):
                    if not model_dir.is_dir():
                        continue
                    for protocol_dir in sorted(model_dir.iterdir()):
                        if not protocol_dir.is_dir() or protocol_dir.name == 'log':
                            continue
                        for pdb_gz in sorted(protocol_dir.glob("*.pdb.gz")):
                            stem = pdb_gz.stem.replace('.pdb', '')
                            rep = stem.split('_r')[-1] if '_r' in stem else stem
                            structures.append({
                                'path': str(pdb_gz), 'protein': pid,
                                'category': 'AlphaFold' if cat_name == 'AF' else 'Boltz',
                                'subcategory': f'relaxed_{protocol_dir.name}',
                                'model': f"{model_dir.name}_r{rep}", 'compressed': True
                            })

    return structures


def parse_pdb(pdb_path: str) -> list:
    """Parse ATOM/HETATM records from PDB file."""
    atoms = []
    try:
        with open(pdb_path, 'r') as f:
            for line in f:
                if not (line.startswith('ATOM') or line.startswith('HETATM')):
                    continue
                try:
                    atom = {
                        'name': line[12:16].strip(),
                        'resname': line[17:20].strip(),
                        'chain': line[21].strip(),
                        'resseq': int(line[22:26]),
                        'x': float(line[30:38]),
                        'y': float(line[38:46]),
                        'z': float(line[46:54]),
                        'element': line[76:78].strip() if len(line) > 76 else line[12:16].strip()[0]
                    }
                    atoms.append(atom)
                except (ValueError, IndexError):
                    continue
    except Exception:
        return []
    return atoms


def dist(a1, a2) -> float:
    return np.sqrt((a1['x']-a2['x'])**2 + (a1['y']-a2['y'])**2 + (a1['z']-a2['z'])**2)


def angle(a1, a2, a3) -> float:
    """Angle at a2 between a1-a2-a3."""
    v1 = np.array([a1['x']-a2['x'], a1['y']-a2['y'], a1['z']-a2['z']])
    v2 = np.array([a3['x']-a2['x'], a3['y']-a2['y'], a3['z']-a2['z']])
    n1, n2 = np.linalg.norm(v1), np.linalg.norm(v2)
    if n1 == 0 or n2 == 0:
        return 0.0
    cos_ang = np.clip(np.dot(v1, v2) / (n1 * n2), -1, 1)
    return np.degrees(np.arccos(cos_ang))


def dihedral(a1, a2, a3, a4) -> float:
    """Dihedral angle defined by four atoms."""
    b1 = np.array([a2['x']-a1['x'], a2['y']-a1['y'], a2['z']-a1['z']])
    b2 = np.array([a3['x']-a2['x'], a3['y']-a2['y'], a3['z']-a2['z']])
    b3 = np.array([a4['x']-a3['x'], a4['y']-a3['y'], a4['z']-a3['z']])

    n1, n2 = np.cross(b1, b2), np.cross(b2, b3)
    n1_norm, n2_norm, b2_norm = np.linalg.norm(n1), np.linalg.norm(n2), np.linalg.norm(b2)

    if n1_norm == 0 or n2_norm == 0 or b2_norm == 0:
        return 0.0

    m1 = np.cross(n1, b2 / b2_norm)
    return np.degrees(np.arctan2(np.dot(m1, n2), np.dot(n1, n2)))


def group_by_residue(atoms: list) -> dict:
    residues = {}
    for atom in atoms:
        key = (atom['chain'], atom['resseq'], atom['resname'])
        if key not in residues:
            residues[key] = {}
        residues[key][atom['name']] = atom
    return residues


def test_structure_loaded(atoms):
    n = len(atoms)
    return {'structure_loaded': n > 0, 'raw_n_atoms': n}


def test_valid_residues(atoms):
    residues = set(a['resname'] for a in atoms)
    non_std = residues - STANDARD_AA
    return {
        'valid_residues': len(non_std) == 0,
        'raw_n_residue_types': len(residues),
        'raw_n_valid_residue_types': len(residues & STANDARD_AA),
        'raw_non_standard_residues': ','.join(sorted(non_std)) if non_std else ''
    }


def test_backbone_connected(atoms):
    residues = {}
    for a in atoms:
        key = (a['chain'], a['resseq'])
        if key not in residues:
            residues[key] = {}
        residues[key][a['name']] = a

    keys = sorted(residues.keys())
    n_breaks, max_break = 0, 0.0

    for i in range(len(keys) - 1):
        k1, k2 = keys[i], keys[i+1]
        if k1[0] != k2[0]:
            continue
        r1, r2 = residues[k1], residues[k2]
        if 'C' in r1 and 'N' in r2:
            d = dist(r1['C'], r2['N'])
            if d > 2.0:
                n_breaks += 1
                max_break = max(max_break, d)

    return {
        'backbone_connected': n_breaks == 0,
        'raw_n_backbone_breaks': n_breaks,
        'raw_max_break_distance': round(max_break, 3)
    }


def test_bond_lengths(atoms):
    residues = {}
    for a in atoms:
        key = (a['chain'], a['resseq'])
        if key not in residues:
            residues[key] = {}
        residues[key][a['name']] = a

    keys = sorted(residues.keys())
    lengths, n_outliers = [], 0

    for i in range(len(keys) - 1):
        k1, k2 = keys[i], keys[i+1]
        if k1[0] != k2[0] or k2[1] - k1[1] != 1:
            continue
        r1, r2 = residues[k1], residues[k2]
        if 'C' in r1 and 'N' in r2:
            d = dist(r1['C'], r2['N'])
            lengths.append(d)
            if d < 1.18 or d > 1.48:
                n_outliers += 1

    return {
        'bond_lengths': n_outliers == 0,
        'raw_n_peptide_bonds': len(lengths),
        'raw_n_bond_outliers': n_outliers,
        'raw_mean_bond_length': round(np.mean(lengths), 4) if lengths else 0,
        'raw_std_bond_length': round(np.std(lengths), 4) if lengths else 0
    }


def test_bond_angles(atoms):
    residues = group_by_residue(atoms)
    angles_list, n_outliers = [], 0

    for res in residues.values():
        if all(a in res for a in ['N', 'CA', 'C']):
            ang = angle(res['N'], res['CA'], res['C'])
            angles_list.append(ang)
            if ang < 100 or ang > 120:
                n_outliers += 1

    return {
        'bond_angles': n_outliers == 0,
        'raw_n_backbone_angles': len(angles_list),
        'raw_n_angle_outliers': n_outliers,
        'raw_mean_backbone_angle': round(np.mean(angles_list), 2) if angles_list else 0,
        'raw_std_backbone_angle': round(np.std(angles_list), 2) if angles_list else 0
    }


def test_steric_clashes(atoms):
    n_atoms = len(atoms)
    n_clashes, worst = 0, 0.0

    sample = atoms if n_atoms <= 1000 else [atoms[i] for i in np.random.choice(n_atoms, 1000, replace=False)]

    for i in range(len(sample)):
        for j in range(i + 1, len(sample)):
            a1, a2 = sample[i], sample[j]
            if a1['chain'] == a2['chain'] and abs(a1['resseq'] - a2['resseq']) <= 1:
                continue

            d = dist(a1, a2)
            r1, r2 = VDW_RADII.get(a1['element'], 1.7), VDW_RADII.get(a2['element'], 1.7)
            overlap = r1 + r2 - d

            if overlap > 0.5:
                n_clashes += 1
                worst = max(worst, overlap)

    return {
        'steric_clashes': n_clashes < 5,
        'raw_n_clashes': n_clashes,
        'raw_worst_clash': round(worst, 3)
    }


def test_aromatic_flatness(atoms):
    residues = group_by_residue(atoms)
    n_aromatic, n_nonplanar, max_dev = 0, 0, 0.0

    for key, res in residues.items():
        if key[2] not in RING_ATOMS:
            continue

        coords = []
        for name in RING_ATOMS[key[2]]:
            if name in res:
                a = res[name]
                coords.append([a['x'], a['y'], a['z']])

        if len(coords) < 4:
            continue

        n_aromatic += 1
        coords = np.array(coords)
        centered = coords - coords.mean(axis=0)
        _, _, vh = np.linalg.svd(centered)
        rmsd = np.sqrt(np.mean(np.dot(centered, vh[2])**2))
        max_dev = max(max_dev, rmsd)

        if rmsd > 0.1:
            n_nonplanar += 1

    return {
        'aromatic_flatness': n_nonplanar == 0,
        'raw_n_aromatic_rings': n_aromatic,
        'raw_n_nonplanar_rings': n_nonplanar,
        'raw_max_ring_deviation': round(max_dev, 4)
    }


def test_peptide_planarity(atoms):
    residues = group_by_residue(atoms)
    keys = sorted(residues.keys())

    n_omega, n_cis, n_trans, n_twisted = 0, 0, 0, 0
    omega_vals = []

    for i in range(len(keys) - 1):
        k1, k2 = keys[i], keys[i+1]
        if k1[0] != k2[0]:
            continue

        r1, r2 = residues[k1], residues[k2]
        if all(a in r1 for a in ['CA', 'C']) and all(a in r2 for a in ['N', 'CA']):
            omega = dihedral(r1['CA'], r1['C'], r2['N'], r2['CA'])
            omega_vals.append(omega)
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
        'raw_mean_abs_omega': round(np.mean(np.abs(omega_vals)), 2) if omega_vals else 0
    }


def test_chirality(atoms):
    residues = group_by_residue(atoms)
    n_checked, n_d = 0, 0

    for key, res in residues.items():
        if key[2] == 'GLY':
            continue
        if all(a in res for a in ['N', 'CA', 'C', 'CB']):
            n_checked += 1
            if dihedral(res['N'], res['CA'], res['C'], res['CB']) < -30:
                n_d += 1

    return {
        'chirality': n_d == 0,
        'raw_n_chiral_centers': n_checked,
        'raw_n_d_amino_acids': n_d
    }


def test_complete_residues(atoms):
    residues = group_by_residue(atoms)
    backbone = ['N', 'CA', 'C', 'O']

    n_incomplete, missing = 0, 0
    for res in residues.values():
        m = sum(1 for a in backbone if a not in res)
        if m > 0:
            n_incomplete += 1
            missing += m

    return {
        'complete_residues': n_incomplete == 0,
        'raw_n_residues': len(residues),
        'raw_n_incomplete_residues': n_incomplete,
        'raw_n_missing_backbone_atoms': missing
    }


def test_internal_energy(pdb_path: str, rosetta_bin: str) -> dict:
    if not rosetta_bin:
        return {'internal_energy': None, 'raw_rosetta_score': None}

    try:
        with tempfile.TemporaryDirectory() as tmpdir:
            score_file = Path(tmpdir) / 'score.sc'
            cmd = [
                rosetta_bin,
                '-in:file:s', str(pdb_path),
                '-out:file:scorefile', str(score_file),
                '-ignore_unrecognized_res', '-mute', 'all'
            ]
            subprocess.run(cmd, capture_output=True, timeout=120)

            if not score_file.exists():
                return {'internal_energy': None, 'raw_rosetta_score': None}

            for line in score_file.read_text().split('\n'):
                if line.startswith('SCORE:') and 'total_score' not in line:
                    parts = line.split()
                    if len(parts) > 1:
                        try:
                            score = float(parts[1])
                            return {'internal_energy': score < 0, 'raw_rosetta_score': round(score, 2)}
                        except ValueError:
                            pass
    except (subprocess.TimeoutExpired, Exception):
        pass

    return {'internal_energy': None, 'raw_rosetta_score': None}


def decompress_pdb(gz_path: str) -> str:
    pdb_name = Path(gz_path).stem + f'_{uuid.uuid4().hex[:8]}'
    pdb_path = TEMP_DIR / pdb_name
    with gzip.open(gz_path, 'rb') as f_in:
        with open(pdb_path, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    return str(pdb_path)


def validate_structure(args) -> dict:
    struct, rosetta_bin = args

    result = {
        'protein': struct['protein'],
        'category': struct['category'],
        'subcategory': struct['subcategory'],
        'model': struct.get('model', 'exp'),
    }

    pdb_path = struct['path']
    temp_pdb = None

    try:
        if struct.get('compressed'):
            pdb_path = decompress_pdb(pdb_path)
            temp_pdb = pdb_path

        atoms = parse_pdb(pdb_path)

        for test_fn in [test_structure_loaded, test_valid_residues, test_backbone_connected,
                        test_bond_lengths, test_bond_angles, test_steric_clashes,
                        test_aromatic_flatness, test_peptide_planarity, test_chirality,
                        test_complete_residues]:
            result.update(test_fn(atoms))

        if rosetta_bin:
            result.update(test_internal_energy(pdb_path, rosetta_bin))
        else:
            result['internal_energy'] = None
            result['raw_rosetta_score'] = None

        pass_cols = ['structure_loaded', 'valid_residues', 'backbone_connected', 'bond_lengths',
                     'bond_angles', 'steric_clashes', 'aromatic_flatness', 'peptide_planarity',
                     'chirality', 'complete_residues']
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
        if temp_pdb and Path(temp_pdb).exists():
            try:
                Path(temp_pdb).unlink()
            except Exception:
                pass

    return result


def sort_results(df: pd.DataFrame) -> pd.DataFrame:
    cat_order = {'Experimental': 0, 'AlphaFold': 1, 'Boltz': 2}
    sub_order = {
        'original': 0, 'raw': 1,
        'relaxed_cartesian_beta': 2, 'relaxed_cartesian_ref15': 3,
        'relaxed_dualspace_beta': 4, 'relaxed_dualspace_ref15': 5,
        'relaxed_normal_beta': 6, 'relaxed_normal_ref15': 7
    }
    df = df.copy()
    df['_cat'] = df['category'].map(cat_order).fillna(99)
    df['_sub'] = df['subcategory'].map(sub_order).fillna(99)
    df = df.sort_values(['_cat', '_sub', 'model']).drop(columns=['_cat', '_sub'])
    return df.reset_index(drop=True)


def save_per_protein(results: list, protein: str):
    df = sort_results(pd.DataFrame(results))
    out_dir = PROTEINS_DIR / protein / "analysis"
    out_dir.mkdir(parents=True, exist_ok=True)

    pass_cols = ['category', 'subcategory', 'model', 'structure_loaded', 'valid_residues',
                 'backbone_connected', 'bond_lengths', 'bond_angles', 'steric_clashes',
                 'aromatic_flatness', 'peptide_planarity', 'chirality', 'complete_residues',
                 'internal_energy', 'all_pass', 'n_pass']
    df[[c for c in pass_cols if c in df.columns]].to_csv(out_dir / "posebusters_results.csv", index=False)

    raw_cols = ['category', 'subcategory', 'model'] + [c for c in df.columns if c.startswith('raw_')]
    df[[c for c in raw_cols if c in df.columns]].to_csv(out_dir / "posebusters_raw.csv", index=False)

    return len(results)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--no-energy', action='store_true')
    parser.add_argument('--limit', type=int)
    parser.add_argument('-j', '--workers', type=int, default=os.cpu_count())
    args = parser.parse_args()

    print("=" * 70)
    print("POSEBUSTERS - Protein Structure Validity Checks")
    print("=" * 70)
    print(f"\nStart: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Workers: {args.workers}")

    rosetta_bin = None if args.no_energy else find_rosetta()
    print(f"Rosetta: {rosetta_bin}" if rosetta_bin else "(Rosetta disabled)")

    structures = find_structures()
    if args.limit:
        structures = structures[:args.limit]
    print(f"Found {len(structures)} structures")

    by_protein = {}
    for s in structures:
        by_protein.setdefault(s['protein'], []).append(s)

    proteins = sorted(by_protein.keys())
    print(f"Processing {len(proteins)} proteins")

    all_results = []

    for idx, protein in enumerate(proteins, 1):
        tasks = [(s, rosetta_bin) for s in by_protein[protein]]
        protein_results = []

        with ProcessPoolExecutor(max_workers=args.workers) as executor:
            futures = {executor.submit(validate_structure, t): t for t in tasks}
            for future in tqdm(as_completed(futures), total=len(futures),
                              desc=f"{protein} ({idx}/{len(proteins)})", leave=False):
                protein_results.append(future.result())
                all_results.append(future.result())

        n = save_per_protein(protein_results, protein)
        print(f"[{idx}/{len(proteins)}] {protein}: {n} structures saved", flush=True)

        # Append to compiled
        df = sort_results(pd.DataFrame(protein_results))
        pass_cols = ['protein', 'category', 'subcategory', 'model', 'structure_loaded', 'valid_residues',
                     'backbone_connected', 'bond_lengths', 'bond_angles', 'steric_clashes',
                     'aromatic_flatness', 'peptide_planarity', 'chirality', 'complete_residues',
                     'internal_energy', 'all_pass', 'n_pass']
        raw_cols = ['protein', 'category', 'subcategory', 'model'] + [c for c in df.columns if c.startswith('raw_')]

        compiled = OUTPUT_DIR / "posebusters_results.csv"
        compiled_raw = OUTPUT_DIR / "posebusters_raw.csv"
        header = not compiled.exists()
        df[[c for c in pass_cols if c in df.columns]].to_csv(compiled, mode='a', header=header, index=False)
        df[[c for c in raw_cols if c in df.columns]].to_csv(compiled_raw, mode='a', header=header, index=False)

    # Summary
    print("\n" + "=" * 70)
    df = pd.DataFrame(all_results)
    print(f"Total: {len(df)} structures")
    print(f"All pass: {df['all_pass'].sum()} ({100*df['all_pass'].mean():.1f}%)")

    print("\nBy category:")
    for cat in ['Experimental', 'AlphaFold', 'Boltz']:
        sub = df[df['category'] == cat]
        if len(sub):
            print(f"  {cat}: {sub['all_pass'].sum()}/{len(sub)} ({100*sub['all_pass'].mean():.1f}%)")

    print("\nTest pass rates:")
    for col in ['structure_loaded', 'valid_residues', 'backbone_connected', 'bond_lengths',
                'bond_angles', 'steric_clashes', 'aromatic_flatness', 'peptide_planarity',
                'chirality', 'complete_residues', 'internal_energy']:
        if col in df.columns:
            valid = df[col].notna()
            if valid.any():
                print(f"  {col}: {100*df.loc[valid, col].mean():.1f}%")

    print(f"\nEnd: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")


if __name__ == "__main__":
    main()
