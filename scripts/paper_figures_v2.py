#!/usr/bin/env python3
"""
Paper figures v2 - aligned with God Claude's original spec.

Original spec:
- Fig 1: Heatmap — clashscore improvement by (protocol) x (prediction method)
- Fig 2: Paired scatter — pre vs post relaxation MolProbity scores per structure
- Fig 3: Convergence — RMSD between crystal-relaxed and AF-relaxed (needs structural data)
- Fig 4: Box plots — per-protocol delta_toward_bound (needs bound structures)
- Table 1: Summary statistics

This script adds what's missing from v1 and can run with current data.
"""

import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

RESULTS_DIR = Path(__file__).parent.parent / "validation_results"
FIGURES_DIR = Path(__file__).parent.parent / "figures"
FIGURES_DIR.mkdir(exist_ok=True)

CATEGORIES = ["Experimental", "AlphaFold", "Boltz"]
RELAXED_PROTOCOLS = [
    "relaxed_normal_beta",
    "relaxed_normal_ref15",
    "relaxed_cartesian_beta",
    "relaxed_cartesian_ref15",
    "relaxed_dualspace_beta",
    "relaxed_dualspace_ref15",
]

KEY_METRICS = ["clashscore", "molprobity_score", "rama_outliers_pct", "rota_outliers_pct"]


def load_data():
    return pd.read_csv(RESULTS_DIR / "molprobity_full.csv")


# =============================================================================
# Fig 1: Clashscore improvement heatmap (protocol x category)
# =============================================================================
def fig1_clashscore_heatmap(df):
    """Heatmap: mean clashscore improvement by protocol and category."""
    improvements = []

    for category in CATEGORIES:
        cat_df = df[df["category"] == category]

        # Baseline
        if category == "Experimental":
            baseline = cat_df[cat_df["subcategory"] == "original"]["clashscore"].mean()
        else:
            baseline = cat_df[cat_df["subcategory"] == "raw"]["clashscore"].mean()

        for protocol in RELAXED_PROTOCOLS:
            relaxed = cat_df[cat_df["subcategory"] == protocol]["clashscore"].mean()
            if pd.notna(baseline) and pd.notna(relaxed):
                improvement = baseline - relaxed  # Positive = improved
                improvements.append({
                    "category": category,
                    "protocol": protocol.replace("relaxed_", ""),
                    "improvement": improvement,
                    "baseline": baseline,
                    "relaxed": relaxed
                })

    if not improvements:
        print("  No data for clashscore heatmap")
        return

    imp_df = pd.DataFrame(improvements)
    pivot = imp_df.pivot(index="protocol", columns="category", values="improvement")

    fig, ax = plt.subplots(figsize=(10, 6))
    sns.heatmap(pivot, annot=True, fmt=".1f", cmap="RdYlGn", center=0, ax=ax,
                cbar_kws={"label": "Clashscore Reduction"})
    ax.set_title("Fig 1: Clashscore Improvement by Protocol and Source\n(positive = fewer clashes)")
    ax.set_ylabel("Relaxation Protocol")
    ax.set_xlabel("Structure Source")

    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "fig1_clashscore_heatmap.png", dpi=150)
    plt.close()
    print("  Saved fig1_clashscore_heatmap.png")

    return pivot


# =============================================================================
# Fig 2: Paired scatter (pre vs post per structure)
# =============================================================================
def fig2_paired_scatter(df):
    """Scatter plot: unrelaxed vs relaxed MolProbity score per structure."""
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    axes = axes.flatten()

    for i, protocol in enumerate(RELAXED_PROTOCOLS):
        ax = axes[i]

        for category in CATEGORIES:
            cat_df = df[df["category"] == category]

            # Get baseline
            if category == "Experimental":
                baseline_df = cat_df[cat_df["subcategory"] == "original"]
            else:
                baseline_df = cat_df[cat_df["subcategory"] == "raw"]

            relaxed_df = cat_df[cat_df["subcategory"] == protocol]

            # Match by protein
            for protein in baseline_df["protein"].unique():
                base_score = baseline_df[baseline_df["protein"] == protein]["molprobity_score"].mean()
                rel_score = relaxed_df[relaxed_df["protein"] == protein]["molprobity_score"].mean()

                if pd.notna(base_score) and pd.notna(rel_score):
                    color = {"Experimental": "green", "AlphaFold": "blue", "Boltz": "red"}[category]
                    ax.scatter(base_score, rel_score, c=color, alpha=0.6, s=30)

        # Diagonal line (no change)
        lims = [0, max(ax.get_xlim()[1], ax.get_ylim()[1])]
        ax.plot(lims, lims, 'k--', alpha=0.5, label='No change')
        ax.set_xlim(lims)
        ax.set_ylim(lims)

        ax.set_xlabel("Unrelaxed MolProbity Score")
        ax.set_ylabel("Relaxed MolProbity Score")
        ax.set_title(protocol.replace("relaxed_", ""))

    # Legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='green', label='Experimental'),
        Patch(facecolor='blue', label='AlphaFold'),
        Patch(facecolor='red', label='Boltz'),
    ]
    fig.legend(handles=legend_elements, loc='upper right', bbox_to_anchor=(0.98, 0.98))

    plt.suptitle("Fig 2: Pre vs Post Relaxation MolProbity Score\n(below diagonal = improved)", fontsize=14)
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "fig2_paired_scatter.png", dpi=150)
    plt.close()
    print("  Saved fig2_paired_scatter.png")


# =============================================================================
# Table 1: Summary statistics
# =============================================================================
def table1_summary_stats(df):
    """Generate summary statistics table."""
    results = []

    for category in CATEGORIES:
        cat_df = df[df["category"] == category]

        # Baseline stats
        if category == "Experimental":
            baseline = cat_df[cat_df["subcategory"] == "original"]
            baseline_name = "original"
        else:
            baseline = cat_df[cat_df["subcategory"] == "raw"]
            baseline_name = "raw"

        results.append({
            "category": category,
            "protocol": baseline_name,
            "n": len(baseline),
            "clashscore_mean": baseline["clashscore"].mean(),
            "clashscore_std": baseline["clashscore"].std(),
            "molprobity_mean": baseline["molprobity_score"].mean(),
            "molprobity_std": baseline["molprobity_score"].std(),
            "rama_outliers_mean": baseline["rama_outliers_pct"].mean(),
            "rota_outliers_mean": baseline["rota_outliers_pct"].mean(),
        })

        # Relaxed stats
        for protocol in RELAXED_PROTOCOLS:
            relaxed = cat_df[cat_df["subcategory"] == protocol]
            if relaxed.empty:
                continue

            results.append({
                "category": category,
                "protocol": protocol.replace("relaxed_", ""),
                "n": len(relaxed),
                "clashscore_mean": relaxed["clashscore"].mean(),
                "clashscore_std": relaxed["clashscore"].std(),
                "molprobity_mean": relaxed["molprobity_score"].mean(),
                "molprobity_std": relaxed["molprobity_score"].std(),
                "rama_outliers_mean": relaxed["rama_outliers_pct"].mean(),
                "rota_outliers_mean": relaxed["rota_outliers_pct"].mean(),
            })

    summary_df = pd.DataFrame(results)
    summary_df.to_csv(FIGURES_DIR / "table1_summary_stats.csv", index=False)
    print("  Saved table1_summary_stats.csv")

    return summary_df


# =============================================================================
# Statistical tests
# =============================================================================
def statistical_tests(df):
    """Run statistical tests comparing protocols."""
    print("\n=== Statistical Tests ===\n")

    results = []

    for category in CATEGORIES:
        cat_df = df[df["category"] == category]

        if category == "Experimental":
            baseline = cat_df[cat_df["subcategory"] == "original"]
        else:
            baseline = cat_df[cat_df["subcategory"] == "raw"]

        for protocol in RELAXED_PROTOCOLS:
            relaxed = cat_df[cat_df["subcategory"] == protocol]

            # Match by protein for paired test
            paired_baseline = []
            paired_relaxed = []

            for protein in baseline["protein"].unique():
                b = baseline[baseline["protein"] == protein]["molprobity_score"].mean()
                r = relaxed[relaxed["protein"] == protein]["molprobity_score"].mean()
                if pd.notna(b) and pd.notna(r):
                    paired_baseline.append(b)
                    paired_relaxed.append(r)

            if len(paired_baseline) >= 3:
                # Wilcoxon signed-rank (non-parametric paired test)
                stat, p_value = stats.wilcoxon(paired_baseline, paired_relaxed)
                mean_diff = np.mean(np.array(paired_baseline) - np.array(paired_relaxed))

                results.append({
                    "category": category,
                    "protocol": protocol.replace("relaxed_", ""),
                    "n_pairs": len(paired_baseline),
                    "mean_improvement": mean_diff,
                    "wilcoxon_stat": stat,
                    "p_value": p_value,
                    "significant": p_value < 0.05
                })

    stats_df = pd.DataFrame(results)
    stats_df.to_csv(FIGURES_DIR / "statistical_tests.csv", index=False)
    print(stats_df.to_string())
    print("\n  Saved statistical_tests.csv")

    return stats_df


# =============================================================================
# Outlier analysis
# =============================================================================
def outlier_analysis(df):
    """Which structures don't improve with relaxation?"""
    outliers = []

    for category in CATEGORIES:
        cat_df = df[df["category"] == category]

        if category == "Experimental":
            baseline = cat_df[cat_df["subcategory"] == "original"]
        else:
            baseline = cat_df[cat_df["subcategory"] == "raw"]

        best_protocol = "relaxed_normal_beta"  # Best from pilot
        relaxed = cat_df[cat_df["subcategory"] == best_protocol]

        for protein in baseline["protein"].unique():
            b = baseline[baseline["protein"] == protein]["molprobity_score"].mean()
            r = relaxed[relaxed["protein"] == protein]["molprobity_score"].mean()

            if pd.notna(b) and pd.notna(r):
                improvement = b - r
                if improvement < 0:  # Got worse
                    outliers.append({
                        "protein": protein,
                        "category": category,
                        "baseline_score": b,
                        "relaxed_score": r,
                        "degradation": -improvement
                    })

    if outliers:
        outlier_df = pd.DataFrame(outliers).sort_values("degradation", ascending=False)
        outlier_df.to_csv(FIGURES_DIR / "outliers_degraded.csv", index=False)
        print(f"\n  Found {len(outliers)} structures that got WORSE after relaxation:")
        print(outlier_df.head(10).to_string())
        print("\n  Saved outliers_degraded.csv")
    else:
        print("\n  No structures degraded after relaxation")


# =============================================================================
# Main
# =============================================================================
def main():
    print("Loading data...")
    df = load_data()
    print(f"Loaded {len(df)} rows\n")

    print("Generating Fig 1: Clashscore Heatmap...")
    pivot = fig1_clashscore_heatmap(df)
    if pivot is not None:
        print("\nClashscore improvement matrix:")
        print(pivot)

    print("\nGenerating Fig 2: Paired Scatter...")
    fig2_paired_scatter(df)

    print("\nGenerating Table 1: Summary Stats...")
    summary = table1_summary_stats(df)

    print("\nRunning Statistical Tests...")
    statistical_tests(df)

    print("\nRunning Outlier Analysis...")
    outlier_analysis(df)

    print(f"\n=== All outputs saved to {FIGURES_DIR} ===")


if __name__ == "__main__":
    main()
