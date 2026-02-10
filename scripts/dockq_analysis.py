#!/usr/bin/env python3
"""
DockQ analysis pipeline for BM5.5 relaxation benchmark.

Computes DockQ scores for predicted structures against bound experimental references.
Analyzes whether relaxation improves docking accuracy (not just MolProbity metrics).

Usage:
    python dockq_analysis.py --predictions <dir> --references <dir> --output <csv>

Requirements:
    - DockQ installed (https://github.com/bjornwallner/DockQ)
    - Bound structures from BM5.5 benchmark
"""

import argparse
import subprocess
import pandas as pd
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
import json
import re


def run_dockq(model_pdb: Path, native_pdb: Path) -> dict:
    """Run DockQ on a single model-native pair."""
    try:
        result = subprocess.run(
            ["DockQ", str(model_pdb), str(native_pdb)],
            capture_output=True,
            text=True,
            timeout=120
        )

        output = result.stdout

        # Parse DockQ output
        dockq_match = re.search(r"DockQ\s+([\d.]+)", output)
        fnat_match = re.search(r"Fnat\s+([\d.]+)", output)
        irms_match = re.search(r"iRMS\s+([\d.]+)", output)
        lrms_match = re.search(r"LRMS\s+([\d.]+)", output)

        return {
            "model": model_pdb.name,
            "native": native_pdb.name,
            "dockq": float(dockq_match.group(1)) if dockq_match else None,
            "fnat": float(fnat_match.group(1)) if fnat_match else None,
            "irms": float(irms_match.group(1)) if irms_match else None,
            "lrms": float(lrms_match.group(1)) if lrms_match else None,
            "success": True,
            "error": None
        }
    except subprocess.TimeoutExpired:
        return {
            "model": model_pdb.name,
            "native": native_pdb.name,
            "dockq": None,
            "fnat": None,
            "irms": None,
            "lrms": None,
            "success": False,
            "error": "timeout"
        }
    except Exception as e:
        return {
            "model": model_pdb.name,
            "native": native_pdb.name,
            "dockq": None,
            "fnat": None,
            "irms": None,
            "lrms": None,
            "success": False,
            "error": str(e)
        }


def find_native_structure(target_id: str, references_dir: Path) -> Path:
    """Find the bound structure for a target."""
    # BM5.5 naming: <PDB>_l_b.pdb (ligand bound) or <PDB>_r_b.pdb (receptor bound)
    # For complex, we need the bound complex
    patterns = [
        f"{target_id}_bound.pdb",
        f"{target_id}.pdb",
        f"{target_id}_complex.pdb",
    ]

    for pattern in patterns:
        candidate = references_dir / pattern
        if candidate.exists():
            return candidate

    # Try subdirectory structure
    subdir = references_dir / target_id
    if subdir.exists():
        for pattern in patterns:
            candidate = subdir / pattern
            if candidate.exists():
                return candidate

    return None


def collect_predictions(predictions_dir: Path) -> list:
    """Collect all prediction files organized by target/source/protocol."""
    predictions = []

    for target_dir in predictions_dir.iterdir():
        if not target_dir.is_dir():
            continue

        target_id = target_dir.name

        # AlphaFold predictions
        af_dir = target_dir / "af_out"
        if af_dir.exists():
            for pdb in af_dir.glob("ranked_*.pdb"):
                predictions.append({
                    "target": target_id,
                    "source": "AlphaFold",
                    "protocol": "AMBER" if "ranked" in pdb.name else "raw",
                    "model": pdb.stem,
                    "path": pdb
                })

        # Unrelaxed AF
        af_unrelaxed = target_dir / "af_out_unrelaxed"
        if af_unrelaxed.exists():
            for pdb in af_unrelaxed.glob("*.pdb"):
                predictions.append({
                    "target": target_id,
                    "source": "AlphaFold",
                    "protocol": "unrelaxed",
                    "model": pdb.stem,
                    "path": pdb
                })

        # Boltz predictions
        boltz_dir = target_dir / "boltz_out"
        if boltz_dir.exists():
            for pdb in boltz_dir.glob("*.pdb"):
                predictions.append({
                    "target": target_id,
                    "source": "Boltz",
                    "protocol": "raw",
                    "model": pdb.stem,
                    "path": pdb
                })

        # Rosetta relaxed structures
        rosetta_dir = target_dir / "rosetta_out"
        if rosetta_dir.exists():
            for protocol_dir in rosetta_dir.iterdir():
                if not protocol_dir.is_dir():
                    continue
                protocol = protocol_dir.name
                for pdb in protocol_dir.glob("*.pdb"):
                    predictions.append({
                        "target": target_id,
                        "source": "Rosetta",
                        "protocol": protocol,
                        "model": pdb.stem,
                        "path": pdb
                    })

    return predictions


def analyze_dockq_results(results_df: pd.DataFrame) -> dict:
    """Compute summary statistics and comparisons."""
    summary = {}

    # Overall stats by source
    for source in results_df["source"].unique():
        source_df = results_df[results_df["source"] == source]
        summary[f"{source}_mean_dockq"] = source_df["dockq"].mean()
        summary[f"{source}_median_dockq"] = source_df["dockq"].median()
        summary[f"{source}_std_dockq"] = source_df["dockq"].std()

    # Protocol comparison
    for protocol in results_df["protocol"].unique():
        protocol_df = results_df[results_df["protocol"] == protocol]
        summary[f"protocol_{protocol}_mean_dockq"] = protocol_df["dockq"].mean()

    # Key question: Does relaxation improve DockQ?
    # Compare raw vs relaxed for each source

    return summary


def main():
    parser = argparse.ArgumentParser(description="DockQ analysis for BM5.5 benchmark")
    parser.add_argument("--predictions", type=Path, required=True,
                        help="Directory containing predicted structures")
    parser.add_argument("--references", type=Path, required=True,
                        help="Directory containing bound reference structures")
    parser.add_argument("--output", type=Path, default=Path("dockq_results.csv"),
                        help="Output CSV file")
    parser.add_argument("--workers", type=int, default=8,
                        help="Number of parallel workers")
    args = parser.parse_args()

    print(f"Collecting predictions from {args.predictions}")
    predictions = collect_predictions(args.predictions)
    print(f"Found {len(predictions)} prediction files")

    results = []

    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        futures = {}

        for pred in predictions:
            native = find_native_structure(pred["target"], args.references)
            if native is None:
                print(f"Warning: No native structure for {pred['target']}")
                continue

            future = executor.submit(run_dockq, pred["path"], native)
            futures[future] = pred

        for future in as_completed(futures):
            pred = futures[future]
            dockq_result = future.result()

            results.append({
                "target": pred["target"],
                "source": pred["source"],
                "protocol": pred["protocol"],
                "model": pred["model"],
                **dockq_result
            })

    # Create DataFrame and save
    df = pd.DataFrame(results)
    df.to_csv(args.output, index=False)
    print(f"Saved {len(df)} results to {args.output}")

    # Print summary
    if not df.empty and df["dockq"].notna().any():
        print("\n=== DockQ Summary ===")
        summary = analyze_dockq_results(df)
        for key, value in summary.items():
            print(f"  {key}: {value:.3f}")


if __name__ == "__main__":
    main()
