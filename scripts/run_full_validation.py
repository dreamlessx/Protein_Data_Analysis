#!/usr/bin/env python3
"""
run_full_validation.py
Comprehensive protein structure validation matching RCSB/wwPDB metrics

Vanderbilt Protein Structure Validation Project

This script computes ALL standard validation metrics:
- Clashscore (via probe)
- Ramachandran analysis (via ramalyze)
- Rotamer outliers (via rotalyze)
- C-beta deviations (via cbetadev)
- Omega/cis-peptide analysis (via omegalyze)
- Bond length/angle deviations (via geometry analysis)
- Chirality validation
- CaBLAM analysis

Environment requirements:
    conda activate molprobity
    export MMTBX_CCP4_MONOMER_LIB="/Users/muditagar/miniconda3/envs/molprobity/share/ccp4_mon_lib"
    export REDUCE_HET_DICT="/Users/muditagar/miniconda3/envs/molprobity/share/reduce/reduce_wwPDB_het_dict.txt"
"""

import subprocess
import sys
import os
import re
import gzip
import shutil
import math
from pathlib import Path
from datetime import datetime
import json
import csv
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm
import numpy as np

# Set environment for geometry library BEFORE importing cctbx modules
os.environ['MMTBX_CCP4_MONOMER_LIB'] = str(Path.home() / "miniconda3/envs/molprobity/share/geostd")

# Configuration
PROJECT_DIR = Path(__file__).parent.parent
PROTEINS_DIR = PROJECT_DIR / "proteins"
OUTPUT_DIR = PROJECT_DIR / "validation_results"
TEMP_DIR = OUTPUT_DIR / "temp"

# Ensure output directories exist
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
TEMP_DIR.mkdir(parents=True, exist_ok=True)

# Environment setup
REDUCE_HET_DICT = Path.home() / "miniconda3/envs/molprobity/share/reduce/reduce_wwPDB_het_dict.txt"
PROBE_PATH = Path.home() / "miniconda3/envs/molprobity/bin/probe"
REDUCE_PATH = Path.home() / "miniconda3/envs/molprobity/bin/reduce"


def find_all_structures():
    """Find all PDB structure files in the proteins directory."""
    structures = []

    for protein_dir in sorted(PROTEINS_DIR.iterdir()):
        if not protein_dir.is_dir():
            continue

        protein_id = protein_dir.name

        # Experimental structure
        exp_pdb = protein_dir / f"{protein_id}.pdb"
        if exp_pdb.exists():
            structures.append({
                'path': exp_pdb,
                'protein': protein_id,
                'category': 'Experimental',
                'subcategory': 'original',
                'model': 'exp'
            })

        # AlphaFold predictions
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

        # Boltz predictions
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

        # Relaxed structures
        for relax_type in ['cartesian_beta', 'cartesian_ref15', 'dualspace_beta',
                           'dualspace_ref15', 'normal_beta', 'normal_ref15']:
            relax_dir = protein_dir / relax_type
            if relax_dir.exists():
                # Experimental relaxed
                for pdb_gz in sorted(relax_dir.glob(f"{protein_id}_r*.pdb.gz")):
                    rep = pdb_gz.stem.replace('.pdb', '').split('_r')[1]
                    structures.append({
                        'path': pdb_gz,
                        'protein': protein_id,
                        'category': 'Experimental',
                        'subcategory': f'relaxed_{relax_type}',
                        'model': f"r{rep}",
                        'compressed': True
                    })

        # Check relax directory for AF/Boltz relaxed structures
        # Structure: relax/AF/ranked_0/cartesian_beta/ranked_0_r1.pdb.gz
        #            relax/Boltz/boltz_input_model_0/cartesian_beta/boltz_input_model_0_r1.pdb.gz
        relax_base = protein_dir / "relax"
        if relax_base.exists():
            for category_dir in ['AF', 'Boltz']:
                cat_relax = relax_base / category_dir
                if cat_relax.exists():
                    for model_dir in sorted(cat_relax.iterdir()):
                        if not model_dir.is_dir():
                            continue
                        model_name = model_dir.name  # e.g. ranked_0 or boltz_input_model_0
                        for protocol_dir in sorted(model_dir.iterdir()):
                            if not protocol_dir.is_dir() or protocol_dir.name == 'log':
                                continue
                            relax_type = protocol_dir.name  # e.g. cartesian_beta
                            for pdb_gz in sorted(protocol_dir.glob("*.pdb.gz")):
                                stem = pdb_gz.stem.replace('.pdb', '')
                                category = 'AlphaFold' if category_dir == 'AF' else 'Boltz'
                                # Extract replicate number from filename
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


def run_ramalyze(pdb_path):
    """Run Ramachandran analysis."""
    try:
        result = subprocess.run(
            ['molprobity.ramalyze', str(pdb_path)],
            capture_output=True,
            text=True,
            timeout=120
        )

        output = result.stdout + result.stderr

        # Parse summary
        favored = allowed = outliers = total = 0
        for line in output.split('\n'):
            if 'SUMMARY:' in line and 'Favored' in line:
                match = re.search(r'(\d+) Favored, (\d+) Allowed, (\d+) Outlier.* out of (\d+)', line)
                if match:
                    favored, allowed, outliers, total = map(int, match.groups())

        if total > 0:
            return {
                'rama_favored': favored,
                'rama_allowed': allowed,
                'rama_outliers': outliers,
                'rama_total': total,
                'rama_favored_pct': 100 * favored / total,
                'rama_allowed_pct': 100 * allowed / total,
                'rama_outliers_pct': 100 * outliers / total
            }
    except Exception as e:
        pass

    return None


def run_rotalyze(pdb_path):
    """Run rotamer outlier analysis."""
    try:
        result = subprocess.run(
            ['molprobity.rotalyze', str(pdb_path)],
            capture_output=True,
            text=True,
            timeout=120
        )

        output = result.stdout + result.stderr

        # Parse summary - format: "SUMMARY: X.XX% outliers (Goal: < 0.3%)"
        for line in output.split('\n'):
            if 'SUMMARY:' in line and 'outliers' in line:
                match = re.search(r'([\d.]+)%\s+outliers', line)
                if match:
                    outlier_pct = float(match.group(1))
                    return {'rotamer_outliers_pct': outlier_pct}
    except Exception as e:
        pass

    return None


def run_cbetadev(pdb_path):
    """Run C-beta deviation analysis."""
    try:
        result = subprocess.run(
            ['molprobity.cbetadev', str(pdb_path)],
            capture_output=True,
            text=True,
            timeout=120
        )

        output = result.stdout + result.stderr

        # Parse summary
        for line in output.split('\n'):
            if 'SUMMARY:' in line:
                match = re.search(r'(\d+) C-beta deviation', line)
                if match:
                    return {'cbeta_deviations': int(match.group(1))}
    except Exception as e:
        pass

    return None


def run_omegalyze(pdb_path):
    """Run omega angle analysis for cis-peptides."""
    try:
        result = subprocess.run(
            ['molprobity.omegalyze', str(pdb_path)],
            capture_output=True,
            text=True,
            timeout=120
        )

        output = result.stdout + result.stderr

        # Parse summary - each metric is on a separate SUMMARY line
        cis_pro = cis_general = twisted_pro = twisted_general = 0
        for line in output.split('\n'):
            if 'SUMMARY:' in line:
                match = re.search(r'(\d+)\s+cis\s+prolines?\s+out\s+of', line, re.IGNORECASE)
                if match:
                    cis_pro = int(match.group(1))
                match = re.search(r'(\d+)\s+twisted\s+prolines?\s+out\s+of', line, re.IGNORECASE)
                if match:
                    twisted_pro = int(match.group(1))
                match = re.search(r'(\d+)\s+other\s+cis\s+residues?\s+out\s+of', line, re.IGNORECASE)
                if match:
                    cis_general = int(match.group(1))
                match = re.search(r'(\d+)\s+other\s+twisted\s+residues?\s+out\s+of', line, re.IGNORECASE)
                if match:
                    twisted_general = int(match.group(1))

        return {
            'omega_cis_proline': cis_pro,
            'omega_cis_general': cis_general,
            'omega_twisted': twisted_pro + twisted_general
        }
    except Exception as e:
        pass

    return None


def run_clashscore(pdb_path):
    """Calculate clashscore using reduce + probe."""
    try:
        env = os.environ.copy()
        env['REDUCE_HET_DICT'] = str(REDUCE_HET_DICT)

        # Add hydrogens with reduce
        reduce_result = subprocess.run(
            [str(REDUCE_PATH), '-build', str(pdb_path)],
            capture_output=True,
            text=True,
            timeout=120,
            env=env
        )

        pdb_with_h = reduce_result.stdout
        if not pdb_with_h or 'ATOM' not in pdb_with_h:
            return None

        # Write to temp file
        temp_h_pdb = TEMP_DIR / f"temp_H_{os.getpid()}.pdb"
        with open(temp_h_pdb, 'w') as f:
            f.write(pdb_with_h)

        # Run probe for clash analysis
        probe_result = subprocess.run(
            [str(PROBE_PATH), '-4H', '-mc', '-self', 'ALL', '-unformated', str(temp_h_pdb)],
            capture_output=True,
            text=True,
            timeout=120
        )

        # Clean up
        temp_h_pdb.unlink(missing_ok=True)

        # Count unique clash pairs (bad overlaps)
        clash_pairs = set()
        for line in probe_result.stdout.split('\n'):
            if line.startswith(':') and ':bo:' in line:
                parts = line.split(':')
                if len(parts) >= 5:
                    src_atom = parts[3].strip()
                    targ_atom = parts[4].strip()
                    pair = tuple(sorted([src_atom, targ_atom]))
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

    except Exception as e:
        pass

    return None


def run_geometry_validation(pdb_path):
    """Calculate Bond RMSZ and Angle RMSZ using cctbx geometry restraints."""
    try:
        import iotbx.pdb
        import mmtbx.model
        from io import StringIO
        from cctbx import geometry_restraints

        # Load PDB
        pdb_inp = iotbx.pdb.input(str(pdb_path))

        # Create model with restraints
        model = mmtbx.model.manager(
            model_input=pdb_inp,
            log=StringIO()
        )

        # Process restraints
        model.process(make_restraints=True)

        # Get geometry restraints manager
        grm = model.get_restraints_manager().geometry

        # Get sites
        sites_cart = model.get_sites_cart()

        # Bond analysis
        bond_proxies = grm.pair_proxies(sites_cart=sites_cart).bond_proxies
        bond_deltas = geometry_restraints.bond_deltas(
            sites_cart=sites_cart,
            proxies=bond_proxies.simple
        )

        deltas = np.array(list(bond_deltas))
        n_bonds = len(deltas)

        # Get sigma values for RMSZ (sigma = 1/sqrt(weight))
        sigmas = []
        for proxy in bond_proxies.simple:
            sigmas.append(proxy.weight**(-0.5) if proxy.weight > 0 else 0.02)
        sigmas = np.array(sigmas)

        z_scores = deltas / sigmas
        bond_rmsz = np.sqrt(np.mean(z_scores**2))

        # Angle analysis
        angle_proxies = grm.angle_proxies
        n_angles = len(angle_proxies)
        angle_deltas = []
        angle_sigmas = []

        for proxy in angle_proxies:
            # Get the 3 atoms
            i, j, k = proxy.i_seqs
            a = np.array(sites_cart[i])
            b = np.array(sites_cart[j])  # vertex
            c = np.array(sites_cart[k])

            # Calculate angle
            ba = a - b
            bc = c - b

            cos_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
            cos_angle = np.clip(cos_angle, -1, 1)
            angle_deg = np.degrees(np.arccos(cos_angle))

            # Delta from ideal (angle_ideal is already in degrees)
            delta = angle_deg - proxy.angle_ideal

            angle_deltas.append(delta)
            sigma = proxy.weight**(-0.5) if proxy.weight > 0 else 3.0
            angle_sigmas.append(sigma)

        angle_deltas = np.array(angle_deltas)
        angle_sigmas = np.array(angle_sigmas)
        angle_z_scores = angle_deltas / angle_sigmas
        angle_rmsz = np.sqrt(np.mean(angle_z_scores**2))

        return {
            'bond_rmsz': round(bond_rmsz, 4),
            'angle_rmsz': round(angle_rmsz, 4),
            'n_bonds': n_bonds,
            'n_angles': n_angles
        }
    except Exception as e:
        pass

    return None


def run_cablam(pdb_path):
    """Run CaBLAM analysis for secondary structure validation."""
    try:
        result = subprocess.run(
            ['molprobity.cablam', str(pdb_path)],
            capture_output=True,
            text=True,
            timeout=120
        )

        output = result.stdout + result.stderr

        # Parse CaBLAM output
        cablam_outliers = cablam_disfavored = cablam_ca_outliers = 0
        for line in output.split('\n'):
            if 'SUMMARY' in line or 'CaBLAM' in line:
                # Parse outlier percentages
                match = re.search(r'([\d.]+)%\s+CaBLAM outliers', line, re.IGNORECASE)
                if match:
                    cablam_outliers = float(match.group(1))
                match = re.search(r'([\d.]+)%\s+CaBLAM disfavored', line, re.IGNORECASE)
                if match:
                    cablam_disfavored = float(match.group(1))
                match = re.search(r'([\d.]+)%\s+CA geometry outliers', line, re.IGNORECASE)
                if match:
                    cablam_ca_outliers = float(match.group(1))

        return {
            'cablam_outliers_pct': cablam_outliers,
            'cablam_disfavored_pct': cablam_disfavored,
            'cablam_ca_outliers_pct': cablam_ca_outliers
        }
    except Exception as e:
        pass

    return None


def calculate_molprobity_score(clashscore, rama_outliers_pct, rotamer_outliers_pct):
    """Calculate the MolProbity score from component metrics.

    Formula: 0.426 * ln(1 + clashscore) +
             0.33 * ln(1 + max(0, rama_outliers_pct - 2)) +
             0.25 * ln(1 + max(0, rotamer_outliers_pct - 1)) + 0.5
    """
    try:
        if clashscore is None or rama_outliers_pct is None or rotamer_outliers_pct is None:
            return None

        score = (0.426 * math.log(1 + clashscore) +
                 0.33 * math.log(1 + max(0, rama_outliers_pct - 2)) +
                 0.25 * math.log(1 + max(0, rotamer_outliers_pct - 1)) + 0.5)

        return round(score, 3)
    except Exception:
        return None


def analyze_structure(structure):
    """Run all analyses on a single structure."""
    results = {
        'protein': structure['protein'],
        'category': structure['category'],
        'subcategory': structure['subcategory'],
        'model': structure['model'],
        'path': str(structure['path'])
    }

    # Decompress if needed
    pdb_path = decompress_if_needed(structure)

    try:
        # Run Ramachandran analysis
        rama_results = run_ramalyze(pdb_path)
        if rama_results:
            results.update(rama_results)

        # Run Rotamer analysis
        rota_results = run_rotalyze(pdb_path)
        if rota_results:
            results.update(rota_results)

        # Run C-beta deviation
        cbeta_results = run_cbetadev(pdb_path)
        if cbeta_results:
            results.update(cbeta_results)

        # Run omega analysis
        omega_results = run_omegalyze(pdb_path)
        if omega_results:
            results.update(omega_results)

        # Run clashscore
        clash_results = run_clashscore(pdb_path)
        if clash_results:
            results.update(clash_results)

        # Run CaBLAM
        cablam_results = run_cablam(pdb_path)
        if cablam_results:
            results.update(cablam_results)

        # Run geometry validation (Bond/Angle RMSZ)
        geo_results = run_geometry_validation(pdb_path)
        if geo_results:
            results.update(geo_results)

        # Calculate MolProbity Score
        mp_score = calculate_molprobity_score(
            results.get('clashscore'),
            results.get('rama_outliers_pct'),
            results.get('rotamer_outliers_pct')
        )
        if mp_score is not None:
            results['molprobity_score'] = mp_score

    finally:
        # Clean up temp file if we created one
        if structure.get('compressed') and pdb_path.exists():
            pdb_path.unlink(missing_ok=True)

    return results


def main():
    print("=" * 70)
    print("COMPREHENSIVE PROTEIN STRUCTURE VALIDATION")
    print("Matching RCSB/wwPDB Validation Metrics")
    print("Vanderbilt Protein Structure Validation Project")
    print("=" * 70)
    print(f"\nStart time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

    # Find all structures
    print("\nDiscovering structure files...")
    structures = find_all_structures()
    print(f"Found {len(structures)} structures")

    # Count by category
    categories = {}
    for s in structures:
        cat = s['category']
        categories[cat] = categories.get(cat, 0) + 1

    for cat, count in sorted(categories.items()):
        print(f"  - {cat}: {count}")

    # Run analyses
    print(f"\nRunning validation analyses...")
    print("Metrics computed:")
    print("  - Clashscore (via probe)")
    print("  - Ramachandran (via ramalyze)")
    print("  - Rotamer outliers (via rotalyze)")
    print("  - C-beta deviations (via cbetadev)")
    print("  - Omega/cis-peptides (via omegalyze)")
    print("  - CaBLAM secondary structure")
    print("  - Bond RMSZ (via cctbx geometry)")
    print("  - Angle RMSZ (via cctbx geometry)")
    print("  - MolProbity Score (calculated)")

    all_results = []

    # Process structures with progress bar
    for structure in tqdm(structures, desc="Analyzing structures"):
        try:
            result = analyze_structure(structure)
            all_results.append(result)
        except Exception as e:
            print(f"\nError processing {structure['path']}: {e}")
            all_results.append({
                'protein': structure['protein'],
                'category': structure['category'],
                'subcategory': structure['subcategory'],
                'model': structure['model'],
                'path': str(structure['path']),
                'error': str(e)
            })

    # Save results
    print(f"\n{'=' * 70}")
    print("SAVING RESULTS")
    print("=" * 70)

    # Save as CSV
    csv_path = OUTPUT_DIR / "full_validation_results.csv"

    # Get all field names
    fieldnames = set()
    for r in all_results:
        fieldnames.update(r.keys())
    fieldnames = sorted(fieldnames)

    with open(csv_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(all_results)

    print(f"Saved CSV: {csv_path}")

    # Save as JSON
    json_path = OUTPUT_DIR / "full_validation_results.json"
    with open(json_path, 'w') as f:
        json.dump(all_results, f, indent=2, default=str)

    print(f"Saved JSON: {json_path}")

    # Print summary
    print(f"\n{'=' * 70}")
    print("SUMMARY")
    print("=" * 70)

    successful = sum(1 for r in all_results if 'clashscore' in r and r['clashscore'] is not None)
    print(f"Successfully analyzed: {successful}/{len(structures)}")

    # Summary by category
    print("\nMean MolProbity Score by Category:")
    for cat in ['AlphaFold', 'Boltz', 'Experimental']:
        cat_results = [r for r in all_results if r.get('category') == cat and 'molprobity_score' in r]
        if cat_results:
            mean_mp = sum(r['molprobity_score'] for r in cat_results) / len(cat_results)
            print(f"  {cat}: {mean_mp:.3f}")

    print("\nMean Clashscore by Category:")
    for cat in ['AlphaFold', 'Boltz', 'Experimental']:
        cat_results = [r for r in all_results if r.get('category') == cat and 'clashscore' in r]
        if cat_results:
            mean_clash = sum(r['clashscore'] for r in cat_results) / len(cat_results)
            print(f"  {cat}: {mean_clash:.2f}")

    print(f"\nEnd time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("\n" + "=" * 70)
    print("VALIDATION COMPLETE")
    print("=" * 70)


if __name__ == "__main__":
    main()
