#!/usr/bin/env python3
"""
Paper figure generation for BM5.5 relaxation benchmark.

5 Key Analyses:
1. Relaxation delta (MolProbity before/after by protocol)
2. Protocol ranking (which of 7 protocols wins?)
3. Convergence test (crystal→relax vs AF→relax)
4. MSA depth effect (full_dbs vs reduced_dbs)
5. Bound vs unbound (does relaxation help docking?)
"""

import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns

# Config
RESULTS_DIR = Path(__file__).parent.parent / "validation_results"
FIGURES_DIR = Path(__file__).parent.parent / "figures"
FIGURES_DIR.mkdir(exist_ok=True)

# Actual data structure from molprobity_full.csv:
# category: Experimental, AlphaFold, Boltz
# subcategory: original/raw (unrelaxed), relaxed_* (6 protocols)

CATEGORIES = ["Experimental", "AlphaFold", "Boltz"]

SUBCATEGORIES = [
    "original",  # Experimental unrelaxed
    "raw",  # AF/Boltz unrelaxed
    "relaxed_normal_beta",
    "relaxed_normal_ref15",
    "relaxed_cartesian_beta",
    "relaxed_cartesian_ref15",
    "relaxed_dualspace_beta",
    "relaxed_dualspace_ref15",
]

RELAXED_PROTOCOLS = [s for s in SUBCATEGORIES if s.startswith("relaxed_")]

METRICS = [
    "rama_favored_pct",
    "rama_outliers_pct",
    "rota_outliers_pct",
    "clashscore",
    "cbeta_outliers",
    "bond_rmsz",
    "angle_rmsz",
    "molprobity_score",
]

# Direction: higher is better (+1) or lower is better (-1)
METRIC_DIRECTION = {
    "rama_favored_pct": 1,
    "rama_outliers_pct": -1,
    "rota_outliers_pct": -1,
    "clashscore": -1,
    "cbeta_outliers": -1,
    "bond_rmsz": -1,
    "angle_rmsz": -1,
    "molprobity_score": -1,
}


def load_validation_data(csv_path: str = None) -> pd.DataFrame:
    """Load consolidated MolProbity results."""
    if csv_path is None:
        csv_path = RESULTS_DIR / "molprobity_full.csv"
    return pd.read_csv(csv_path)


# =============================================================================
# Figure 1: Relaxation Delta
# =============================================================================
def figure1_relaxation_delta(df: pd.DataFrame):
    """
    Box plots showing MolProbity metric changes after relaxation.
    Compares each relaxed protocol against its unrelaxed baseline.
    """
    fig, axes = plt.subplots(2, 4, figsize=(16, 8))
    axes = axes.flatten()

    for i, metric in enumerate(METRICS):
        ax = axes[i]
        deltas = []

        for category in CATEGORIES:
            cat_df = df[df["category"] == category]

            # Get unrelaxed baseline
            if category == "Experimental":
                unrelaxed = cat_df[cat_df["subcategory"] == "original"]
            else:
                unrelaxed = cat_df[cat_df["subcategory"] == "raw"]

            if unrelaxed.empty:
                continue

            # Calculate deltas for each relaxed protocol
            for protocol in RELAXED_PROTOCOLS:
                relaxed = cat_df[cat_df["subcategory"] == protocol]
                if relaxed.empty:
                    continue

                # Match by protein and model where possible
                for protein in unrelaxed["protein"].unique():
                    unrel_vals = unrelaxed[unrelaxed["protein"] == protein][metric].mean()
                    rel_vals = relaxed[relaxed["protein"] == protein][metric].mean()
                    if pd.notna(unrel_vals) and pd.notna(rel_vals):
                        delta = rel_vals - unrel_vals
                        # Flip sign so positive = improvement
                        if METRIC_DIRECTION[metric] == -1:
                            delta = -delta
                        deltas.append({
                            "protocol": protocol.replace("relaxed_", ""),
                            "category": category,
                            "delta": delta
                        })

        if deltas:
            delta_df = pd.DataFrame(deltas)
            sns.boxplot(data=delta_df, x="protocol", y="delta", hue="category", ax=ax)
            ax.axhline(0, color="red", linestyle="--", alpha=0.5)
            ax.legend(fontsize=6)
        ax.set_title(metric.replace("_pct", "%").replace("_", " "))
        ax.tick_params(axis="x", rotation=45)
        ax.set_xlabel("")

    plt.suptitle("Figure 1: Relaxation Improvement (positive = better)", fontsize=14)
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "fig1_relaxation_delta.png", dpi=150)
    plt.close()
    print("  Saved fig1_relaxation_delta.png")


# =============================================================================
# Figure 2: Protocol Ranking
# =============================================================================
def figure2_protocol_ranking(df: pd.DataFrame):
    """
    Heatmap showing mean metric values by subcategory (protocol).
    Protocols ranked by composite score.
    """
    # Aggregate by subcategory
    agg = df.groupby("subcategory")[METRICS].mean()

    # Normalize each metric to 0-1 scale where 1 = best
    normalized = agg.copy()
    for col in METRICS:
        if METRIC_DIRECTION[col] == 1:  # Higher is better
            normalized[col] = (agg[col] - agg[col].min()) / (agg[col].max() - agg[col].min())
        else:  # Lower is better
            normalized[col] = 1 - (agg[col] - agg[col].min()) / (agg[col].max() - agg[col].min())

    # Composite score
    normalized["composite"] = normalized[METRICS].mean(axis=1)
    normalized = normalized.sort_values("composite", ascending=False)

    # Create heatmap
    fig, axes = plt.subplots(1, 2, figsize=(16, 8), gridspec_kw={"width_ratios": [3, 1]})

    # Left: normalized scores heatmap
    display_cols = [m.replace("_pct", "%").replace("_", " ") for m in METRICS]
    heatmap_data = normalized[METRICS].copy()
    heatmap_data.columns = display_cols
    sns.heatmap(heatmap_data, annot=True, fmt=".2f", cmap="RdYlGn", ax=axes[0],
                cbar_kws={"label": "Score (1=best)"})
    axes[0].set_title("Protocol Ranking by Metric (normalized)")
    axes[0].set_ylabel("Protocol")

    # Right: composite bar chart
    axes[1].barh(normalized.index, normalized["composite"], color="steelblue")
    axes[1].set_xlabel("Composite Score")
    axes[1].set_title("Overall Ranking")
    axes[1].set_xlim(0, 1)

    plt.suptitle("Figure 2: Protocol Comparison", fontsize=14)
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "fig2_protocol_ranking.png", dpi=150)
    plt.close()

    print("  Saved fig2_protocol_ranking.png")
    return normalized


# =============================================================================
# Figure 3: Convergence Test
# =============================================================================
def figure3_convergence(df: pd.DataFrame):
    """
    Do crystal and AF/Boltz structures converge to same geometry after relaxation?
    Compare experimental→relax vs predicted→relax endpoint distributions.
    """
    fig, axes = plt.subplots(2, 4, figsize=(16, 8))
    axes = axes.flatten()

    # Filter to relaxed only
    relaxed = df[df["subcategory"].str.startswith("relaxed_")]

    for i, metric in enumerate(METRICS):
        ax = axes[i]

        exp_vals = relaxed[relaxed["category"] == "Experimental"][metric].dropna()
        af_vals = relaxed[relaxed["category"] == "AlphaFold"][metric].dropna()
        boltz_vals = relaxed[relaxed["category"] == "Boltz"][metric].dropna()

        if len(exp_vals) > 0:
            ax.hist(exp_vals, alpha=0.5, label=f"Exp (n={len(exp_vals)})", bins=30, density=True)
        if len(af_vals) > 0:
            ax.hist(af_vals, alpha=0.5, label=f"AF (n={len(af_vals)})", bins=30, density=True)
        if len(boltz_vals) > 0:
            ax.hist(boltz_vals, alpha=0.5, label=f"Boltz (n={len(boltz_vals)})", bins=30, density=True)

        ax.legend(fontsize=7)
        ax.set_title(metric.replace("_pct", "%").replace("_", " "))
        ax.set_ylabel("Density")

    plt.suptitle("Figure 3: Convergence - Do different sources converge after relaxation?", fontsize=14)
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "fig3_convergence.png", dpi=150)
    plt.close()
    print("  Saved fig3_convergence.png")


# =============================================================================
# Figure 4: Predictor Comparison (AlphaFold vs Boltz)
# =============================================================================
def figure4_predictor_comparison(df: pd.DataFrame):
    """
    Compare AlphaFold vs Boltz predictions before and after relaxation.
    (MSA depth analysis requires BM5.5 full run with both full_dbs and reduced_dbs)
    """
    fig, axes = plt.subplots(2, 4, figsize=(16, 8))
    axes = axes.flatten()

    for i, metric in enumerate(METRICS):
        ax = axes[i]

        data_to_plot = []
        labels = []

        for category in ["AlphaFold", "Boltz"]:
            cat_df = df[df["category"] == category]

            # Unrelaxed
            raw = cat_df[cat_df["subcategory"] == "raw"][metric].dropna()
            if len(raw) > 0:
                data_to_plot.append(raw)
                labels.append(f"{category[:2]}_raw")

            # Best relaxed (cartesian_ref15 as example)
            relaxed = cat_df[cat_df["subcategory"] == "relaxed_cartesian_ref15"][metric].dropna()
            if len(relaxed) > 0:
                data_to_plot.append(relaxed)
                labels.append(f"{category[:2]}_relax")

        if data_to_plot:
            bp = ax.boxplot(data_to_plot, tick_labels=labels, patch_artist=True)
            colors = ["#1f77b4", "#aec7e8", "#ff7f0e", "#ffbb78"]
            for patch, color in zip(bp["boxes"], colors[:len(bp["boxes"])]):
                patch.set_facecolor(color)

        ax.set_title(metric.replace("_pct", "%").replace("_", " "))
        ax.tick_params(axis="x", rotation=45)

    plt.suptitle("Figure 4: AlphaFold vs Boltz (raw vs relaxed)", fontsize=14)
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "fig4_predictor_comparison.png", dpi=150)
    plt.close()
    print("  Saved fig4_predictor_comparison.png")


# =============================================================================
# Figure 5: Protocol Improvement Summary
# =============================================================================
def figure5_protocol_improvement(df: pd.DataFrame):
    """
    Summary: which protocol improves structures the most?
    Shows improvement in MolProbity score (lower = better).

    Note: For BM5.5 docking analysis, will need RMSD to bound structures.
    """
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Left: MolProbity score by category and protocol
    ax1 = axes[0]
    plot_data = []
    for category in CATEGORIES:
        cat_df = df[df["category"] == category]

        # Get unrelaxed baseline
        if category == "Experimental":
            baseline = cat_df[cat_df["subcategory"] == "original"]["molprobity_score"].mean()
        else:
            baseline = cat_df[cat_df["subcategory"] == "raw"]["molprobity_score"].mean()

        for protocol in RELAXED_PROTOCOLS:
            relaxed_score = cat_df[cat_df["subcategory"] == protocol]["molprobity_score"].mean()
            if pd.notna(baseline) and pd.notna(relaxed_score):
                improvement = baseline - relaxed_score  # Positive = improved
                plot_data.append({
                    "category": category,
                    "protocol": protocol.replace("relaxed_", ""),
                    "improvement": improvement
                })

    if plot_data:
        plot_df = pd.DataFrame(plot_data)
        pivot = plot_df.pivot(index="protocol", columns="category", values="improvement")
        pivot.plot(kind="bar", ax=ax1, width=0.8)
        ax1.axhline(0, color="black", linestyle="-", linewidth=0.5)
        ax1.set_ylabel("MolProbity Score Improvement")
        ax1.set_xlabel("Protocol")
        ax1.set_title("Improvement vs Unrelaxed")
        ax1.tick_params(axis="x", rotation=45)
        ax1.legend(title="Source")

    # Right: Best protocol by category
    ax2 = axes[1]
    summary = df.groupby(["category", "subcategory"])["molprobity_score"].mean().reset_index()
    summary = summary[summary["subcategory"].str.startswith("relaxed_")]

    best = summary.loc[summary.groupby("category")["molprobity_score"].idxmin()]
    colors = {"Experimental": "#2ecc71", "AlphaFold": "#3498db", "Boltz": "#e74c3c"}

    bars = ax2.bar(best["category"], best["molprobity_score"],
                   color=[colors.get(c, "gray") for c in best["category"]])
    ax2.set_ylabel("Best MolProbity Score")
    ax2.set_title("Best Achievable Score by Source")

    for bar, protocol in zip(bars, best["subcategory"]):
        height = bar.get_height()
        ax2.annotate(protocol.replace("relaxed_", ""),
                     xy=(bar.get_x() + bar.get_width() / 2, height),
                     xytext=(0, 3), textcoords="offset points",
                     ha="center", va="bottom", fontsize=8, rotation=45)

    plt.suptitle("Figure 5: Protocol Improvement Summary", fontsize=14)
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "fig5_protocol_improvement.png", dpi=150)
    plt.close()
    print("  Saved fig5_protocol_improvement.png")


# =============================================================================
# Main
# =============================================================================
def generate_all_figures(csv_path: str = None):
    """Generate all 5 paper figures."""
    print("Loading data...")
    df = load_validation_data(csv_path)
    print(f"Loaded {len(df)} rows")
    print(f"Categories: {df['category'].unique().tolist()}")
    print(f"Subcategories: {df['subcategory'].unique().tolist()}")
    print()

    print("Generating Figure 1: Relaxation Delta...")
    figure1_relaxation_delta(df)

    print("Generating Figure 2: Protocol Ranking...")
    rankings = figure2_protocol_ranking(df)
    print("\nTop protocols by composite score:")
    print(rankings[["composite"]].head(5))
    print()

    print("Generating Figure 3: Convergence...")
    figure3_convergence(df)

    print("Generating Figure 4: Predictor Comparison...")
    figure4_predictor_comparison(df)

    print("Generating Figure 5: Protocol Improvement...")
    figure5_protocol_improvement(df)

    print(f"\nAll figures saved to {FIGURES_DIR}")


if __name__ == "__main__":
    generate_all_figures()
