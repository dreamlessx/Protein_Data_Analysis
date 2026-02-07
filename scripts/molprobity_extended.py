#!/usr/bin/env python3
"""
Extended MolProbity geometry metrics via Python API.
Adds C-beta deviation statistics and omega angle distributions.
"""

import os
import gzip
import tempfile
import pandas as pd
import numpy as np
from pathlib import Path
from datetime import datetime
from concurrent.futures import ProcessPoolExecutor, as_completed

os.environ['CLIBD_MON'] = os.path.expanduser('~/miniconda3/envs/molprobity/chem_data/mon_lib')

ROOT = Path(__file__).parent.parent
PROTEINS = ROOT / "proteins"
WORKERS = 12

RELAX_PROTOCOLS = ['cartesian_beta', 'cartesian_ref15', 'dualspace_beta',
                   'dualspace_ref15', 'normal_beta', 'normal_ref15']


def get_metrics(pdb_path):
    """Extract C-beta and omega metrics via mmtbx."""
    from mmtbx.validation import cbetadev, omegalyze
    from iotbx import pdb

    out = {}
    try:
        hierarchy = pdb.input(str(pdb_path)).construct_hierarchy()

        # C-beta deviations
        cb = cbetadev.cbetadev(pdb_hierarchy=hierarchy, outliers_only=False)
        devs = [r.deviation for r in cb.results if hasattr(r, 'deviation')]
        if devs:
            out['cbeta_mean'] = round(np.mean(devs), 4)
            out['cbeta_max'] = round(max(devs), 4)
            out['cbeta_std'] = round(np.std(devs), 4)
            out['cbeta_rmsz'] = round(np.sqrt(np.mean([d**2 for d in devs])) / 0.25, 3)
        out['cbeta_n'] = len(devs)
        out['cbeta_outliers'] = cb.n_outliers

        # Omega angles
        om = omegalyze.omegalyze(pdb_hierarchy=hierarchy, nontrans_only=False)
        vals = [r.omega for r in om.results if hasattr(r, 'omega') and r.omega is not None]
        if vals:
            out['omega_mean'] = round(np.mean(vals), 2)
            out['omega_std'] = round(np.std(vals), 2)
            out['omega_min'] = round(min(vals), 2)
            out['omega_max'] = round(max(vals), 2)
            out['omega_trans'] = sum(1 for o in vals if abs(abs(o) - 180) < 30)
            out['omega_cis'] = sum(1 for o in vals if abs(o) < 30)

    except Exception as e:
        out['error'] = str(e)[:50]

    return out


def process(info):
    """Process single structure."""
    path, protein, category, subcategory, model, gz = info
    result = {'protein': protein, 'category': category,
              'subcategory': subcategory, 'model': model}

    tmp = None
    try:
        if gz:
            tmp = tempfile.NamedTemporaryFile(suffix='.pdb', delete=False)
            with gzip.open(path, 'rb') as f:
                tmp.write(f.read())
            tmp.close()
            pdb_path = tmp.name
        else:
            pdb_path = path

        result.update(get_metrics(pdb_path))
    except Exception as e:
        result['error'] = str(e)[:50]
    finally:
        if tmp and os.path.exists(tmp.name):
            os.unlink(tmp.name)

    return result


def find_structures(pdb_id):
    """Find all structure variants."""
    pdir = PROTEINS / pdb_id
    structs = []

    exp = pdir / f"{pdb_id}.pdb"
    if exp.exists():
        structs.append((str(exp), pdb_id, 'Experimental', 'original', 'exp', False))

    for name, cat, pat in [('AF', 'AlphaFold', 'ranked_*.pdb'),
                            ('Boltz', 'Boltz', 'boltz_input_model_*.pdb')]:
        d = pdir / name
        if d.exists():
            for f in d.glob(pat):
                structs.append((str(f), pdb_id, cat, 'raw', f.stem, False))

    for proto in RELAX_PROTOCOLS:
        d = pdir / proto
        if d.exists():
            for f in d.glob(f"{pdb_id}_r*.pdb.gz"):
                structs.append((str(f), pdb_id, 'Experimental',
                               f'relaxed_{proto}', f.stem.replace('.pdb', ''), True))

    relax = pdir / "relax"
    if relax.exists():
        for name, cat in [('AF', 'AlphaFold'), ('Boltz', 'Boltz')]:
            d = relax / name
            if not d.exists():
                continue
            for f in d.rglob("*.pdb.gz"):
                proto = next((p for p in f.parts if p in RELAX_PROTOCOLS), 'unknown')
                structs.append((str(f), pdb_id, cat,
                               f'relaxed_{proto}', f.stem.replace('.pdb', ''), True))

    return structs


def main():
    print("=" * 60)
    print("MolProbity Extended Metrics")
    print("=" * 60)
    t0 = datetime.now()
    print(f"Start: {t0.strftime('%H:%M:%S')}")

    proteins = sorted(d.name for d in PROTEINS.iterdir() if d.is_dir())
    print(f"Proteins: {len(proteins)}")

    all_structs = []
    for pid in proteins:
        all_structs.extend(find_structures(pid))
    print(f"Structures: {len(all_structs)}")

    results = []
    with ProcessPoolExecutor(max_workers=WORKERS) as ex:
        futs = {ex.submit(process, s): s for s in all_structs}
        done = 0
        for fut in as_completed(futs):
            try:
                results.append(fut.result())
            except Exception:
                pass
            done += 1
            if done % 500 == 0:
                print(f"  {done}/{len(all_structs)}")

    df = pd.DataFrame(results)
    out = ROOT / "validation_results" / "molprobity_extended.csv"
    df.to_csv(out, index=False)
    print(f"\nSaved: {out}")

    # quick summary
    for cat in ['Experimental', 'AlphaFold', 'Boltz']:
        sub = df[(df['category'] == cat) & (df['subcategory'].isin(['original', 'raw']))]
        if len(sub) and 'cbeta_mean' in sub.columns and sub['cbeta_mean'].notna().any():
            print(f"{cat}: cbeta_mean={sub['cbeta_mean'].mean():.4f}A")

    print(f"Done: {datetime.now().strftime('%H:%M:%S')}")


if __name__ == "__main__":
    main()
