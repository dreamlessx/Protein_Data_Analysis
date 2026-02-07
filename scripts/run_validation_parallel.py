#!/usr/bin/env python3
"""
run_validation_parallel.py
Parallelized MolProbity validation pipeline with per-protein organization

Uses multiprocessing for 10x+ speedup over sequential processing.

Output: proteins/{PDB}/analysis/
    molprobity_results.csv
    VALIDATION_SUMMARY.md

Usage:
    conda activate molprobity
    python scripts/run_validation_parallel.py
"""

import subprocess
import os
import re
import gzip
import math
import pandas as pd
from pathlib import Path
from datetime import datetime
from concurrent.futures import ProcessPoolExecutor, as_completed
from multiprocessing import cpu_count
import warnings

warnings.filterwarnings('ignore')

PROJECT_DIR = Path(__file__).parent.parent
PROTEINS_DIR = PROJECT_DIR / "proteins"
REDUCE_HET_DICT = str(Path.home() / "miniconda3/envs/molprobity/share/reduce/reduce_wwPDB_het_dict.txt")
CLIBD_MON = str(Path.home() / "miniconda3/envs/molprobity/chem_data/mon_lib")

# Set environment for MolProbity tools
os.environ['CLIBD_MON'] = CLIBD_MON

N_WORKERS = 12


def find_structures_for_protein(protein_id):
    """Find all structures for a protein."""
    protein_dir = PROTEINS_DIR / protein_id
    structures = []

    relax_protocols = ['cartesian_beta', 'cartesian_ref15', 'dualspace_beta',
                       'dualspace_ref15', 'normal_beta', 'normal_ref15']

    # Experimental original
    exp_pdb = protein_dir / f"{protein_id}.pdb"
    if exp_pdb.exists():
        structures.append({
            'path': str(exp_pdb), 'protein': protein_id,
            'category': 'Experimental', 'subcategory': 'original',
            'model': 'exp', 'compressed': False
        })

    # AlphaFold raw
    af_dir = protein_dir / "AF"
    if af_dir.exists():
        for pdb in af_dir.glob("ranked_*.pdb"):
            structures.append({
                'path': str(pdb), 'protein': protein_id,
                'category': 'AlphaFold', 'subcategory': 'raw',
                'model': pdb.stem, 'compressed': False
            })

    # Boltz raw
    boltz_dir = protein_dir / "Boltz"
    if boltz_dir.exists():
        for pdb in boltz_dir.glob("boltz_input_model_*.pdb"):
            structures.append({
                'path': str(pdb), 'protein': protein_id,
                'category': 'Boltz', 'subcategory': 'raw',
                'model': pdb.stem, 'compressed': False
            })

    # Experimental relaxed
    for protocol in relax_protocols:
        relax_dir = protein_dir / protocol
        if relax_dir.exists():
            for pdb_gz in relax_dir.glob(f"{protein_id}_r*.pdb.gz"):
                structures.append({
                    'path': str(pdb_gz), 'protein': protein_id,
                    'category': 'Experimental', 'subcategory': f'relaxed_{protocol}',
                    'model': pdb_gz.stem.replace('.pdb', ''), 'compressed': True
                })

    # AF/Boltz relaxed
    relax_base = protein_dir / "relax"
    if relax_base.exists():
        for cat_name in ['AF', 'Boltz']:
            cat_dir = relax_base / cat_name
            if not cat_dir.exists():
                continue
            category = 'AlphaFold' if cat_name == 'AF' else 'Boltz'

            for pdb_gz in cat_dir.rglob("*.pdb.gz"):
                rel_path = pdb_gz.relative_to(cat_dir)
                parts = list(rel_path.parts)
                protocol = 'unknown'
                for p in parts:
                    if p in relax_protocols:
                        protocol = p
                        break

                structures.append({
                    'path': str(pdb_gz), 'protein': protein_id,
                    'category': category, 'subcategory': f'relaxed_{protocol}',
                    'model': pdb_gz.stem.replace('.pdb', ''), 'compressed': True
                })

    return structures


def validate_structure(struct_dict):
    """Validate a single structure - runs in worker process."""
    import tempfile

    path = struct_dict['path']
    compressed = struct_dict['compressed']

    results = {
        'protein': struct_dict['protein'],
        'category': struct_dict['category'],
        'subcategory': struct_dict['subcategory'],
        'model': struct_dict['model']
    }

    temp_pdb = None
    try:
        # Decompress if needed
        if compressed:
            temp_pdb = tempfile.NamedTemporaryFile(suffix='.pdb', delete=False)
            with gzip.open(path, 'rb') as f_in:
                temp_pdb.write(f_in.read())
            temp_pdb.close()
            pdb_path = temp_pdb.name
        else:
            pdb_path = path

        # Run MolProbity tools
        results.update(run_molprobity_tools(pdb_path))

    except Exception as e:
        results['error'] = str(e)
    finally:
        if temp_pdb and os.path.exists(temp_pdb.name):
            os.unlink(temp_pdb.name)

    return results


def run_molprobity_tools(pdb_path):
    """Run all MolProbity analyses."""
    results = {}

    # Ramalyze
    try:
        r = subprocess.run(['molprobity.ramalyze', pdb_path],
                          capture_output=True, text=True, timeout=60)
        out = r.stdout + r.stderr
        for line in out.split('\n'):
            if 'SUMMARY:' in line:
                m = re.search(r'(\d+) Favored, (\d+) Allowed, (\d+) Outlier.* out of (\d+)', line)
                if m:
                    fav, allow, out_n, tot = map(int, m.groups())
                    results['rama_favored'] = fav
                    results['rama_allowed'] = allow
                    results['rama_outliers'] = out_n
                    results['rama_total'] = tot
                    results['rama_favored_pct'] = round(100*fav/tot, 2) if tot else 0
                    results['rama_outliers_pct'] = round(100*out_n/tot, 2) if tot else 0
    except Exception:
        pass

    # Rotalyze (rotamer analysis)
    try:
        r = subprocess.run(['molprobity.rotalyze', pdb_path],
                          capture_output=True, text=True, timeout=60)
        out = r.stdout + r.stderr
        for line in out.split('\n'):
            if 'SUMMARY:' in line:
                m = re.search(r'([0-9.]+)%\s+outliers', line)
                if m:
                    results['rota_outliers_pct'] = float(m.group(1))
        # Count rotamers from output lines
        rota_total = 0
        rota_outliers = 0
        rota_favored = 0
        for line in out.split('\n'):
            if line.strip() and not line.startswith('SUMMARY') and ':' in line:
                parts = line.split(':')
                if len(parts) >= 7:
                    rota_total += 1
                    line_upper = line.upper()
                    if 'OUTLIER' in line_upper:
                        rota_outliers += 1
                    elif 'FAVORED' in line_upper:
                        rota_favored += 1
        if rota_total > 0:
            results['rota_total'] = rota_total
            results['rota_favored'] = rota_favored
            results['rota_outliers'] = rota_outliers
            results['rota_favored_pct'] = round(100*rota_favored/rota_total, 2)
    except Exception:
        pass

    # Cbetadev
    try:
        r = subprocess.run(['molprobity.cbetadev', pdb_path],
                          capture_output=True, text=True, timeout=60)
        for line in (r.stdout + r.stderr).split('\n'):
            if 'SUMMARY:' in line:
                m = re.search(r'(\d+) C-beta deviation', line)
                if m:
                    results['cbeta_deviations'] = int(m.group(1))
    except Exception:
        pass

    # Omegalyze
    try:
        r = subprocess.run(['molprobity.omegalyze', pdb_path],
                          capture_output=True, text=True, timeout=60)
        out = r.stdout + r.stderr
        cis_pro = cis_gen = twisted = 0
        for line in out.split('\n'):
            if 'SUMMARY:' in line:
                m = re.search(r'(\d+)\s+cis\s+prolines?', line, re.I)
                if m: cis_pro = int(m.group(1))
                m = re.search(r'(\d+)\s+twisted\s+prolines?', line, re.I)
                if m: twisted += int(m.group(1))
                m = re.search(r'(\d+)\s+other\s+cis', line, re.I)
                if m: cis_gen = int(m.group(1))
                m = re.search(r'(\d+)\s+other\s+twisted', line, re.I)
                if m: twisted += int(m.group(1))
        results['omega_cis_proline'] = cis_pro
        results['omega_cis_general'] = cis_gen
        results['omega_twisted'] = twisted
    except Exception:
        pass

    # Clashscore via reduce+probe
    try:
        env = os.environ.copy()
        env['REDUCE_HET_DICT'] = REDUCE_HET_DICT
        r = subprocess.run(['reduce', '-build', pdb_path],
                          capture_output=True, text=True, timeout=60, env=env)
        pdb_h = r.stdout
        if pdb_h and 'ATOM' in pdb_h:
            import tempfile
            with tempfile.NamedTemporaryFile(suffix='.pdb', mode='w', delete=False) as f:
                f.write(pdb_h)
                temp_h = f.name

            probe_path = '/private/tmp/probe/probe' if os.path.exists('/private/tmp/probe/probe') else 'probe'
            r = subprocess.run([probe_path, '-4H', '-mc', '-self', 'ALL', '-unformated', temp_h],
                              capture_output=True, text=True, timeout=60)
            os.unlink(temp_h)

            clashes = set()
            for line in r.stdout.split('\n'):
                if line.startswith(':') and ':bo:' in line:
                    parts = line.split(':')
                    if len(parts) >= 5:
                        clashes.add(tuple(sorted([parts[3].strip(), parts[4].strip()])))

            atom_count = pdb_h.count('\nATOM ') + pdb_h.count('\nHETATM ')
            if atom_count > 0:
                results['clashscore'] = round(len(clashes)*1000/atom_count, 2)
                results['clash_count'] = len(clashes)
                results['atom_count'] = atom_count
    except Exception:
        pass

    # Calculate MolProbity score
    # Formula: 0.426*ln(1+clashscore) + 0.33*ln(1+max(0,rota_out-1)) + 0.25*ln(1+max(0,rama_out-2))
    try:
        clashscore = results.get('clashscore', 0)
        rota_out = results.get('rota_outliers_pct', 0)
        rama_out = results.get('rama_outliers_pct', 0)

        mp_score = (0.426 * math.log(1 + clashscore) +
                   0.33 * math.log(1 + max(0, rota_out - 1)) +
                   0.25 * math.log(1 + max(0, rama_out - 2)))
        results['molprobity_score'] = round(mp_score, 2)
    except Exception:
        pass

    return results


def generate_summary(protein_id, df, output_path):
    """Generate markdown validation summary."""

    # Define consistent ordering
    CATEGORY_ORDER = ['Experimental', 'AlphaFold', 'Boltz']
    SUBCATEGORY_ORDER = [
        'original', 'raw',
        'relaxed_normal_ref15', 'relaxed_normal_beta',
        'relaxed_cartesian_ref15', 'relaxed_cartesian_beta',
        'relaxed_dualspace_ref15', 'relaxed_dualspace_beta'
    ]

    def sort_subcategories(subs):
        """Sort subcategories in defined order."""
        return sorted(subs, key=lambda x: SUBCATEGORY_ORDER.index(x) if x in SUBCATEGORY_ORDER else 999)

    lines = [
        f"# Validation Summary: {protein_id}",
        f"",
        f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        f"",
        f"## Overview",
        f"",
        f"| Metric | Value |",
        f"|--------|-------|",
        f"| Total Structures | {len(df)} |",
    ]

    for cat in CATEGORY_ORDER:
        if cat in df['category'].values:
            lines.append(f"| {cat} | {len(df[df['category']==cat])} |")

    if 'rama_favored_pct' in df.columns:
        lines.extend([
            f"",
            f"## Ramachandran Quality",
            f"",
            f"| Category | Subcategory | Avg Favored % | Count |",
            f"|----------|-------------|---------------|-------|"
        ])
        for cat in CATEGORY_ORDER:
            cat_df = df[df['category'] == cat]
            if cat_df.empty:
                continue
            for sub in sort_subcategories(cat_df['subcategory'].unique()):
                sub_df = cat_df[cat_df['subcategory'] == sub]
                if 'rama_favored_pct' in sub_df.columns:
                    avg = sub_df['rama_favored_pct'].mean()
                    lines.append(f"| {cat} | {sub} | {avg:.2f} | {len(sub_df)} |")

    if 'rota_favored_pct' in df.columns:
        lines.extend([
            f"",
            f"## Rotamer Quality",
            f"",
            f"| Category | Subcategory | Avg Favored % | Avg Outliers % | Count |",
            f"|----------|-------------|---------------|----------------|-------|"
        ])
        for cat in CATEGORY_ORDER:
            cat_df = df[df['category'] == cat]
            if cat_df.empty:
                continue
            for sub in sort_subcategories(cat_df['subcategory'].unique()):
                sub_df = cat_df[cat_df['subcategory'] == sub]
                if 'rota_favored_pct' in sub_df.columns and sub_df['rota_favored_pct'].notna().any():
                    avg_fav = sub_df['rota_favored_pct'].mean()
                    avg_out = sub_df['rota_outliers_pct'].mean() if 'rota_outliers_pct' in sub_df else 0
                    lines.append(f"| {cat} | {sub} | {avg_fav:.2f} | {avg_out:.2f} | {len(sub_df)} |")

    if 'clashscore' in df.columns:
        lines.extend([
            f"",
            f"## Clashscore (lower is better)",
            f"",
            f"| Category | Subcategory | Avg Clashscore | Count |",
            f"|----------|-------------|----------------|-------|"
        ])
        for cat in CATEGORY_ORDER:
            cat_df = df[df['category'] == cat]
            if cat_df.empty:
                continue
            for sub in sort_subcategories(cat_df['subcategory'].unique()):
                sub_df = cat_df[cat_df['subcategory'] == sub]
                if 'clashscore' in sub_df.columns and sub_df['clashscore'].notna().any():
                    avg = sub_df['clashscore'].mean()
                    lines.append(f"| {cat} | {sub} | {avg:.2f} | {len(sub_df)} |")

    if 'molprobity_score' in df.columns:
        lines.extend([
            f"",
            f"## MolProbity Score (lower is better)",
            f"",
            f"| Category | Subcategory | Avg Score | Count |",
            f"|----------|-------------|-----------|-------|"
        ])
        for cat in CATEGORY_ORDER:
            cat_df = df[df['category'] == cat]
            if cat_df.empty:
                continue
            for sub in sort_subcategories(cat_df['subcategory'].unique()):
                sub_df = cat_df[cat_df['subcategory'] == sub]
                if 'molprobity_score' in sub_df.columns and sub_df['molprobity_score'].notna().any():
                    avg = sub_df['molprobity_score'].mean()
                    lines.append(f"| {cat} | {sub} | {avg:.2f} | {len(sub_df)} |")

    with open(output_path, 'w') as f:
        f.write('\n'.join(lines))


def process_protein(protein_id, skip_existing=True):
    """Process all structures for one protein."""
    protein_dir = PROTEINS_DIR / protein_id
    analysis_dir = protein_dir / "analysis"

    # Skip if already completed
    if skip_existing and (analysis_dir / "molprobity_results.csv").exists():
        return protein_id, 0, True  # skipped

    analysis_dir.mkdir(exist_ok=True)

    structures = find_structures_for_protein(protein_id)
    if not structures:
        return protein_id, 0, False

    # Process structures in parallel
    results = []
    with ProcessPoolExecutor(max_workers=N_WORKERS) as executor:
        futures = {executor.submit(validate_structure, s): s for s in structures}
        for future in as_completed(futures):
            try:
                results.append(future.result())
            except Exception as e:
                results.append({'error': str(e), **futures[future]})

    df = pd.DataFrame(results)

    # MolProbity columns only
    mp_cols = ['protein', 'category', 'subcategory', 'model',
               'rama_favored', 'rama_allowed', 'rama_outliers', 'rama_total',
               'rama_favored_pct', 'rama_outliers_pct',
               'rota_favored', 'rota_outliers', 'rota_total',
               'rota_favored_pct', 'rota_outliers_pct',
               'cbeta_deviations', 'omega_cis_proline', 'omega_cis_general', 'omega_twisted',
               'clashscore', 'clash_count', 'atom_count', 'molprobity_score']
    mp_cols = [c for c in mp_cols if c in df.columns]

    df[mp_cols].to_csv(analysis_dir / "molprobity_results.csv", index=False)
    generate_summary(protein_id, df, analysis_dir / "VALIDATION_SUMMARY.md")

    return protein_id, len(structures), False  # not skipped


def main():
    print("=" * 70)
    print("MOLPROBITY VALIDATION PIPELINE")
    print(f"Workers: {N_WORKERS}")
    print("=" * 70)
    print(f"Start: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

    proteins = sorted([d.name for d in PROTEINS_DIR.iterdir() if d.is_dir()])
    print(f"Found {len(proteins)} proteins")

    total = 0
    skipped = 0
    for protein_id in proteins:
        pid, count, was_skipped = process_protein(protein_id)
        if was_skipped:
            skipped += 1
            print(f"  {pid}: skipped (already done)")
        else:
            total += count
            print(f"  {pid}: {count} structures")

    print(f"\nTotal: {total} structures validated")
    print(f"End: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")


if __name__ == "__main__":
    main()
