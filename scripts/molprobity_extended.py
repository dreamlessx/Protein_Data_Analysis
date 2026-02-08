#!/usr/bin/env python3
"""
Extended MolProbity geometry metrics via Python API.
Adds C-beta deviation, omega distributions, and bond/angle RMSZ.
"""

import os
import gzip
import tempfile
import math
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

# Engh & Huber ideal geometry values (used by wwPDB)
IDEAL_BONDS = {
    ('N', 'CA'): (1.458, 0.019),
    ('CA', 'C'): (1.525, 0.021),
    ('C', 'O'): (1.231, 0.020),
    ('CA', 'CB'): (1.530, 0.020),
}

IDEAL_ANGLES = {
    ('N', 'CA', 'C'): (111.2, 2.8),
    ('CA', 'C', 'O'): (120.8, 1.7),
    ('N', 'CA', 'CB'): (110.5, 1.7),
    ('CB', 'CA', 'C'): (110.1, 1.9),
}


def calc_angle(a1, a2, a3):
    """Calculate angle between three atoms in degrees."""
    v1 = np.array([a1.xyz[i] - a2.xyz[i] for i in range(3)])
    v2 = np.array([a3.xyz[i] - a2.xyz[i] for i in range(3)])
    cos_angle = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
    return math.degrees(math.acos(np.clip(cos_angle, -1, 1)))


def get_bond_angle_rmsz(hierarchy):
    """Calculate bond and angle RMSZ using Engh & Huber values."""
    bond_zs = []
    angle_zs = []

    for model in hierarchy.models():
        for chain in model.chains():
            for rg in chain.residue_groups():
                for conf in rg.conformers():
                    for res in conf.residues():
                        atoms = {a.name.strip(): a for a in res.atoms()}

                        for (a1, a2), (ideal, esd) in IDEAL_BONDS.items():
                            if a1 in atoms and a2 in atoms:
                                dist = atoms[a1].distance(atoms[a2])
                                bond_zs.append((dist - ideal) / esd)

                        for (a1, a2, a3), (ideal, esd) in IDEAL_ANGLES.items():
                            if a1 in atoms and a2 in atoms and a3 in atoms:
                                angle = calc_angle(atoms[a1], atoms[a2], atoms[a3])
                                angle_zs.append((angle - ideal) / esd)

    out = {}
    if bond_zs:
        out['bond_rmsz'] = round(math.sqrt(sum(z**2 for z in bond_zs) / len(bond_zs)), 3)
        out['bond_outliers'] = sum(1 for z in bond_zs if abs(z) > 4)
        out['bond_n'] = len(bond_zs)
    if angle_zs:
        out['angle_rmsz'] = round(math.sqrt(sum(z**2 for z in angle_zs) / len(angle_zs)), 3)
        out['angle_outliers'] = sum(1 for z in angle_zs if abs(z) > 4)
        out['angle_n'] = len(angle_zs)

    return out


def get_metrics(pdb_path):
    """Extract extended geometry metrics."""
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

        # Bond/angle RMSZ
        out.update(get_bond_angle_rmsz(hierarchy))

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
            tmp = tempfile.NamedTemporaryFile(suffix='.pdb', delete=False, mode='wb')
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

    # summary
    print("\n=== Summary (raw/original only) ===")
    for cat in ['Experimental', 'AlphaFold', 'Boltz']:
        sub = df[(df['category'] == cat) & (df['subcategory'].isin(['original', 'raw']))]
        if len(sub) == 0:
            continue
        print(f"\n{cat} (n={len(sub)}):")
        if 'cbeta_mean' in sub.columns and sub['cbeta_mean'].notna().any():
            print(f"  C-beta mean: {sub['cbeta_mean'].mean():.4f} A")
        if 'bond_rmsz' in sub.columns and sub['bond_rmsz'].notna().any():
            print(f"  Bond RMSZ: {sub['bond_rmsz'].mean():.3f}")
        if 'angle_rmsz' in sub.columns and sub['angle_rmsz'].notna().any():
            print(f"  Angle RMSZ: {sub['angle_rmsz'].mean():.3f}")

    print(f"\nDone: {datetime.now().strftime('%H:%M:%S')}")


if __name__ == "__main__":
    main()
