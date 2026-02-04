#!/usr/bin/env python3
"""
run_molprobity.py
Run MolProbity validation on all protein structures

Vanderbilt Protein Structure Validation Project

This script runs available MolProbity tools on all structures:
- Ramachandran analysis (ramalyze)
- C-beta deviations (cbetadev)
- Omega angle analysis (omegalyze)
- Clashscore (via probe)

Usage:
    conda activate molprobity
    python scripts/run_molprobity.py
"""

import subprocess
import sys
import os
import re
import gzip
import shutil
from pathlib import Path
from datetime import datetime
import json
import csv
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm

# Configuration
PROJECT_DIR = Path(__file__).parent.parent
PROTEINS_DIR = PROJECT_DIR / "proteins"
OUTPUT_DIR = PROJECT_DIR / "molprobity_results"
TEMP_DIR = OUTPUT_DIR / "temp"

# Ensure output directories exist
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
TEMP_DIR.mkdir(parents=True, exist_ok=True)

# Environment setup for reduce
REDUCE_HET_DICT = Path.home() / "miniconda3/envs/molprobity/share/reduce/reduce_wwPDB_het_dict.txt"


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

                # Check for AF/Boltz relaxed in relax subdirectory
                af_relax = relax_dir if not (protein_dir / "relax").exists() else None

        # Check relax directory for AF/Boltz relaxed structures
        relax_base = protein_dir / "relax"
        if relax_base.exists():
            for category_dir in ['AF', 'Boltz']:
                cat_relax = relax_base / category_dir
                if cat_relax.exists():
                    for relax_type_dir in cat_relax.iterdir():
                        if relax_type_dir.is_dir():
                            relax_type = relax_type_dir.name
                            for pdb_gz in sorted(relax_type_dir.glob("*.pdb.gz")):
                                # Parse filename to get model and replicate
                                stem = pdb_gz.stem.replace('.pdb', '')
                                category = 'AlphaFold' if category_dir == 'AF' else 'Boltz'
                                structures.append({
                                    'path': pdb_gz,
                                    'protein': protein_id,
                                    'category': category,
                                    'subcategory': f'relaxed_{relax_type}',
                                    'model': stem,
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
            elif 'SUMMARY:' in line and '% outliers' in line:
                match = re.search(r'([\d.]+)% outliers', line)
                if match:
                    outliers_pct = float(match.group(1))
            elif 'SUMMARY:' in line and '% favored' in line:
                match = re.search(r'([\d.]+)% favored', line)
                if match:
                    favored_pct = float(match.group(1))

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
                # Parse cis prolines: "SUMMARY: X cis prolines out of Y PRO"
                match = re.search(r'(\d+)\s+cis\s+prolines?\s+out\s+of', line, re.IGNORECASE)
                if match:
                    cis_pro = int(match.group(1))
                # Parse twisted prolines: "SUMMARY: X twisted prolines out of Y PRO"
                match = re.search(r'(\d+)\s+twisted\s+prolines?\s+out\s+of', line, re.IGNORECASE)
                if match:
                    twisted_pro = int(match.group(1))
                # Parse other cis: "SUMMARY: X other cis residues out of Y nonPRO"
                match = re.search(r'(\d+)\s+other\s+cis\s+residues?\s+out\s+of', line, re.IGNORECASE)
                if match:
                    cis_general = int(match.group(1))
                # Parse other twisted: "SUMMARY: X other twisted residues out of Y nonPRO"
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


def run_clashscore_via_probe(pdb_path):
    """Calculate clashscore using reduce + probe directly."""
    try:
        # First add hydrogens with reduce
        env = os.environ.copy()
        env['REDUCE_HET_DICT'] = str(REDUCE_HET_DICT)

        reduce_result = subprocess.run(
            ['reduce', '-build', str(pdb_path)],
            capture_output=True,
            text=True,
            timeout=120,
            env=env
        )

        # Get the PDB with hydrogens from stdout
        pdb_with_h = reduce_result.stdout

        if not pdb_with_h or 'ATOM' not in pdb_with_h:
            return None

        # Write to temp file
        temp_h_pdb = TEMP_DIR / f"temp_H_{os.getpid()}.pdb"
        with open(temp_h_pdb, 'w') as f:
            f.write(pdb_with_h)

        # Run probe for clash analysis with correct flags
        # -4H: extend bond chain dot removal to 4 for H
        # -mc: include mainchain-mainchain interactions
        # -self: self intersection
        # ALL: all atoms
        # -unformated: raw output for parsing
        probe_result = subprocess.run(
            ['/private/tmp/probe/probe', '-4H', '-mc', '-self', 'ALL', '-unformated', str(temp_h_pdb)],
            capture_output=True,
            text=True,
            timeout=120
        )

        # Clean up
        temp_h_pdb.unlink(missing_ok=True)

        # Parse probe output - count unique atom pairs with bad overlaps (:bo:)
        # Format: :name:type:srcAtom:targAtom:...
        # Type "bo" = bad overlap (clash)
        clash_pairs = set()
        for line in probe_result.stdout.split('\n'):
            if line.startswith(':') and ':bo:' in line:
                parts = line.split(':')
                if len(parts) >= 5:
                    # Create unique pair identifier (sorted to avoid counting A-B and B-A separately)
                    src_atom = parts[3].strip()
                    targ_atom = parts[4].strip()
                    pair = tuple(sorted([src_atom, targ_atom]))
                    clash_pairs.add(pair)

        clash_count = len(clash_pairs)

        # Count atoms for normalization
        atom_count = pdb_with_h.count('\nATOM ') + pdb_with_h.count('\nHETATM ')

        # Calculate clashscore (clashes per 1000 atoms)
        if atom_count > 0:
            clashscore = (clash_count * 1000) / atom_count
            return {'clashscore': round(clashscore, 2), 'clash_count': clash_count, 'atom_count': atom_count}

    except Exception as e:
        pass

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

        # Run C-beta deviation
        cbeta_results = run_cbetadev(pdb_path)
        if cbeta_results:
            results.update(cbeta_results)

        # Run omega analysis
        omega_results = run_omegalyze(pdb_path)
        if omega_results:
            results.update(omega_results)

        # Run clashscore via probe
        clash_results = run_clashscore_via_probe(pdb_path)
        if clash_results:
            results.update(clash_results)

    finally:
        # Clean up temp file if we created one
        if structure.get('compressed') and pdb_path.exists():
            pdb_path.unlink(missing_ok=True)

    return results


def main():
    print("=" * 70)
    print("MOLPROBITY VALIDATION ANALYSIS")
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
    print(f"\nRunning MolProbity analyses...")
    print("This may take a while for 577+ structures...")

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
    csv_path = OUTPUT_DIR / "molprobity_validation_results.csv"

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
    json_path = OUTPUT_DIR / "molprobity_validation_results.json"
    with open(json_path, 'w') as f:
        json.dump(all_results, f, indent=2, default=str)

    print(f"Saved JSON: {json_path}")

    # Print summary
    print(f"\n{'=' * 70}")
    print("SUMMARY")
    print("=" * 70)

    successful = sum(1 for r in all_results if 'rama_total' in r)
    print(f"Successfully analyzed: {successful}/{len(structures)}")

    # Summary by category
    print("\nMean Ramachandran Favored % by Category:")
    for cat in ['AlphaFold', 'Boltz', 'Experimental']:
        cat_results = [r for r in all_results if r.get('category') == cat and 'rama_favored_pct' in r]
        if cat_results:
            mean_favored = sum(r['rama_favored_pct'] for r in cat_results) / len(cat_results)
            print(f"  {cat}: {mean_favored:.2f}%")

    print(f"\nEnd time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("\n" + "=" * 70)
    print("MOLPROBITY ANALYSIS COMPLETE")
    print("=" * 70)


if __name__ == "__main__":
    main()
