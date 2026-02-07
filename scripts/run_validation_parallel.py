#!/usr/bin/env python3
"""
MolProbity validation pipeline for protein structure analysis.
Processes experimental, AlphaFold, and Boltz structures with relaxation variants.

Output: proteins/{PDB}/analysis/{molprobity_results.csv, VALIDATION_SUMMARY.md}
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
import warnings

warnings.filterwarnings('ignore')

ROOT = Path(__file__).parent.parent
PROTEINS = ROOT / "proteins"
REDUCE_DICT = str(Path.home() / "miniconda3/envs/molprobity/share/reduce/reduce_wwPDB_het_dict.txt")
os.environ['CLIBD_MON'] = str(Path.home() / "miniconda3/envs/molprobity/chem_data/mon_lib")

WORKERS = 12
RELAX_PROTOCOLS = ['cartesian_beta', 'cartesian_ref15', 'dualspace_beta',
                   'dualspace_ref15', 'normal_beta', 'normal_ref15']

CAT_ORDER = ['Experimental', 'AlphaFold', 'Boltz']
SUB_ORDER = ['original', 'raw'] + [f'relaxed_{p}' for p in
             ['normal_ref15', 'normal_beta', 'cartesian_ref15',
              'cartesian_beta', 'dualspace_ref15', 'dualspace_beta']]


def find_structures(pdb_id):
    pdir = PROTEINS / pdb_id
    structs = []

    exp = pdir / f"{pdb_id}.pdb"
    if exp.exists():
        structs.append({'path': str(exp), 'protein': pdb_id,
                       'category': 'Experimental', 'subcategory': 'original',
                       'model': 'exp', 'gz': False})

    for name, cat, pattern in [('AF', 'AlphaFold', 'ranked_*.pdb'),
                                ('Boltz', 'Boltz', 'boltz_input_model_*.pdb')]:
        d = pdir / name
        if d.exists():
            for f in d.glob(pattern):
                structs.append({'path': str(f), 'protein': pdb_id,
                               'category': cat, 'subcategory': 'raw',
                               'model': f.stem, 'gz': False})

    for proto in RELAX_PROTOCOLS:
        d = pdir / proto
        if d.exists():
            for f in d.glob(f"{pdb_id}_r*.pdb.gz"):
                structs.append({'path': str(f), 'protein': pdb_id,
                               'category': 'Experimental',
                               'subcategory': f'relaxed_{proto}',
                               'model': f.stem.replace('.pdb', ''), 'gz': True})

    relax = pdir / "relax"
    if relax.exists():
        for name, cat in [('AF', 'AlphaFold'), ('Boltz', 'Boltz')]:
            d = relax / name
            if not d.exists():
                continue
            for f in d.rglob("*.pdb.gz"):
                proto = next((p for p in f.parts if p in RELAX_PROTOCOLS), 'unknown')
                structs.append({'path': str(f), 'protein': pdb_id,
                               'category': cat, 'subcategory': f'relaxed_{proto}',
                               'model': f.stem.replace('.pdb', ''), 'gz': True})

    return structs


def run_mp_tools(pdb_path):
    out = {}

    try:
        r = subprocess.run(['molprobity.ramalyze', pdb_path],
                          capture_output=True, text=True, timeout=60)
        for line in (r.stdout + r.stderr).split('\n'):
            if 'SUMMARY:' in line:
                m = re.search(r'(\d+) Favored, (\d+) Allowed, (\d+) Outlier.* out of (\d+)', line)
                if m:
                    fav, allow, outl, tot = map(int, m.groups())
                    out.update({'rama_favored': fav, 'rama_allowed': allow,
                               'rama_outliers': outl, 'rama_total': tot,
                               'rama_favored_pct': round(100*fav/tot, 2),
                               'rama_outliers_pct': round(100*outl/tot, 2)})
    except Exception:
        pass

    try:
        r = subprocess.run(['molprobity.rotalyze', pdb_path],
                          capture_output=True, text=True, timeout=60)
        txt = r.stdout + r.stderr
        total = fav = outl = 0
        for line in txt.split('\n'):
            if line.strip() and not line.startswith('SUMMARY') and ':' in line:
                if len(line.split(':')) >= 7:
                    total += 1
                    upper = line.upper()
                    if 'OUTLIER' in upper:
                        outl += 1
                    elif 'FAVORED' in upper:
                        fav += 1
        if total:
            out.update({'rota_total': total, 'rota_favored': fav, 'rota_outliers': outl,
                       'rota_favored_pct': round(100*fav/total, 2),
                       'rota_outliers_pct': round(100*outl/total, 2)})
    except Exception:
        pass

    try:
        r = subprocess.run(['molprobity.cbetadev', pdb_path],
                          capture_output=True, text=True, timeout=60)
        for line in (r.stdout + r.stderr).split('\n'):
            if 'SUMMARY:' in line:
                m = re.search(r'(\d+) C-beta deviation', line)
                if m:
                    out['cbeta_deviations'] = int(m.group(1))
    except Exception:
        pass

    try:
        r = subprocess.run(['molprobity.omegalyze', pdb_path],
                          capture_output=True, text=True, timeout=60)
        cis_pro = cis_gen = twist = 0
        for line in (r.stdout + r.stderr).split('\n'):
            if 'SUMMARY:' in line:
                m = re.search(r'(\d+)\s+cis\s+prolines?', line, re.I)
                if m: cis_pro = int(m.group(1))
                m = re.search(r'(\d+)\s+twisted\s+prolines?', line, re.I)
                if m: twist += int(m.group(1))
                m = re.search(r'(\d+)\s+other\s+cis', line, re.I)
                if m: cis_gen = int(m.group(1))
                m = re.search(r'(\d+)\s+other\s+twisted', line, re.I)
                if m: twist += int(m.group(1))
        out.update({'omega_cis_proline': cis_pro, 'omega_cis_general': cis_gen,
                   'omega_twisted': twist})
    except Exception:
        pass

    try:
        env = os.environ.copy()
        env['REDUCE_HET_DICT'] = REDUCE_DICT
        r = subprocess.run(['reduce', '-build', pdb_path],
                          capture_output=True, text=True, timeout=60, env=env)
        pdb_h = r.stdout
        if pdb_h and 'ATOM' in pdb_h:
            import tempfile
            with tempfile.NamedTemporaryFile(suffix='.pdb', mode='w', delete=False) as f:
                f.write(pdb_h)
                tmp = f.name

            probe = 'probe'
            r = subprocess.run([probe, '-4H', '-mc', '-self', 'ALL', '-unformated', tmp],
                              capture_output=True, text=True, timeout=60)
            os.unlink(tmp)

            clashes = set()
            for line in r.stdout.split('\n'):
                if line.startswith(':') and ':bo:' in line:
                    parts = line.split(':')
                    if len(parts) >= 5:
                        clashes.add(tuple(sorted([parts[3].strip(), parts[4].strip()])))

            atoms = pdb_h.count('\nATOM ') + pdb_h.count('\nHETATM ')
            if atoms:
                out.update({'clashscore': round(len(clashes)*1000/atoms, 2),
                           'clash_count': len(clashes), 'atom_count': atoms})
    except Exception:
        pass

    try:
        cs = out.get('clashscore', 0)
        ro = out.get('rota_outliers_pct', 0)
        ra = out.get('rama_outliers_pct', 0)
        score = (0.426 * math.log(1 + cs) +
                 0.33 * math.log(1 + max(0, ro - 1)) +
                 0.25 * math.log(1 + max(0, ra - 2)))
        out['molprobity_score'] = round(score, 2)
    except Exception:
        pass

    return out


def validate(s):
    import tempfile

    result = {k: s[k] for k in ['protein', 'category', 'subcategory', 'model']}
    tmp = None

    try:
        if s['gz']:
            tmp = tempfile.NamedTemporaryFile(suffix='.pdb', delete=False, mode='wb')
            with gzip.open(s['path'], 'rb') as f:
                tmp.write(f.read())
            tmp.close()
            path = tmp.name
        else:
            path = s['path']

        result.update(run_mp_tools(path))
    except Exception as e:
        result['error'] = str(e)
    finally:
        if tmp and os.path.exists(tmp.name):
            os.unlink(tmp.name)

    return result


def write_summary(pdb_id, df, path):
    def sort_sub(subs):
        return sorted(subs, key=lambda x: SUB_ORDER.index(x) if x in SUB_ORDER else 999)

    lines = [f"# {pdb_id} Validation", "",
             f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}",
             "", "## Counts", "",
             "| Metric | N |", "|--------|---|",
             f"| Total | {len(df)} |"]

    for cat in CAT_ORDER:
        n = len(df[df['category'] == cat])
        if n:
            lines.append(f"| {cat} | {n} |")

    metrics = [
        ('rama_favored_pct', 'Ramachandran', 'Favored %'),
        ('rota_favored_pct', 'Rotamer', 'Favored %'),
        ('clashscore', 'Clashscore', 'Score'),
        ('molprobity_score', 'MolProbity', 'Score'),
    ]

    for col, title, label in metrics:
        if col not in df.columns or not df[col].notna().any():
            continue
        lines.extend(["", f"## {title}", "",
                     f"| Category | Subcategory | {label} | N |",
                     "|----------|-------------|---------|---|"])
        for cat in CAT_ORDER:
            cdf = df[df['category'] == cat]
            if cdf.empty:
                continue
            for sub in sort_sub(cdf['subcategory'].unique()):
                sdf = cdf[cdf['subcategory'] == sub]
                if sdf[col].notna().any():
                    lines.append(f"| {cat} | {sub} | {sdf[col].mean():.2f} | {len(sdf)} |")

    path.write_text('\n'.join(lines))


def process(pdb_id, skip_done=True):
    analysis = PROTEINS / pdb_id / "analysis"

    if skip_done and (analysis / "molprobity_results.csv").exists():
        return pdb_id, 0, True

    analysis.mkdir(exist_ok=True)
    structs = find_structures(pdb_id)
    if not structs:
        return pdb_id, 0, False

    results = []
    with ProcessPoolExecutor(max_workers=WORKERS) as ex:
        futs = {ex.submit(validate, s): s for s in structs}
        for fut in as_completed(futs):
            try:
                results.append(fut.result())
            except Exception as e:
                results.append({'error': str(e), **futs[fut]})

    df = pd.DataFrame(results)
    cols = ['protein', 'category', 'subcategory', 'model',
            'rama_favored', 'rama_allowed', 'rama_outliers', 'rama_total',
            'rama_favored_pct', 'rama_outliers_pct',
            'rota_favored', 'rota_outliers', 'rota_total',
            'rota_favored_pct', 'rota_outliers_pct',
            'cbeta_deviations', 'omega_cis_proline', 'omega_cis_general', 'omega_twisted',
            'clashscore', 'clash_count', 'atom_count', 'molprobity_score']
    cols = [c for c in cols if c in df.columns]

    df[cols].to_csv(analysis / "molprobity_results.csv", index=False)
    write_summary(pdb_id, df, analysis / "VALIDATION_SUMMARY.md")

    return pdb_id, len(structs), False


def main():
    print("=" * 60)
    print("MolProbity Validation Pipeline")
    print("=" * 60)
    t0 = datetime.now()
    print(f"Start: {t0.strftime('%H:%M:%S')}")

    proteins = sorted(d.name for d in PROTEINS.iterdir() if d.is_dir())
    print(f"Proteins: {len(proteins)}")

    total = skipped = 0
    for pid in proteins:
        _, n, skip = process(pid)
        if skip:
            skipped += 1
            print(f"  {pid}: cached")
        else:
            total += n
            print(f"  {pid}: {n}")

    print(f"\nValidated: {total} structures")
    print(f"Skipped: {skipped} proteins (cached)")
    print(f"Done: {datetime.now().strftime('%H:%M:%S')}")


if __name__ == "__main__":
    main()
