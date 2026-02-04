# Comprehensive Validation Analysis Report

## Vanderbilt Protein Structure Validation Project
**Generated:** 2026-01-28 23:34:25

---

## Executive Summary

- **Total Structures Analyzed:** 820
- **Unique Proteins:** 20
- **Categories:** AlphaFold (100), Boltz (100), Experimental (620)
- **Validation Pipeline:** Deterministic (3 trials, all identical)

---

## Statistical Analysis Files

| File | Description |
|------|-------------|
| `stats_1_descriptive.csv` | Descriptive statistics for all metrics |
| `stats_2_category_comparisons.csv` | Mann-Whitney U tests between categories |
| `stats_3a_pearson_correlations.csv` | Pearson correlation matrix |
| `stats_3b_spearman_correlations.csv` | Spearman correlation matrix |
| `stats_3c_significant_correlations.csv` | Strong correlations (|r| > 0.3) |
| `stats_4_relaxation_effects.csv` | Before/after relaxation comparisons |
| `stats_5_outliers.csv` | Extreme values and outliers |
| `stats_6_protein_level.csv` | Per-protein analysis |
| `stats_8_anomalies.csv` | Unexpected findings |
| `stats_9_effect_sizes.csv` | Cohen's d, Cliff's delta |
| `stats_10_mp_decomposition.csv` | MolProbity score components |

---

## Key Statistical Findings

### 1. Category Differences Are Significant
- Kruskal-Wallis tests show significant differences (p < 0.001) between AlphaFold, Boltz, and Experimental for most metrics
- Boltz has significantly higher clashscores than both other categories (large effect size)
- Boltz has significantly better Ramachandran geometry (medium effect size)

### 2. Strong Correlations Found
- Clashscore ↔ MolProbity Score: r ≈ 0.7 (strong positive)
- Rama favored ↔ Rama outliers: r ≈ -0.9 (strong negative, expected)
- Bond RMSZ ↔ Angle RMSZ: moderate correlation
- Atom count has weak correlation with quality metrics

### 3. Relaxation Effects
- All relaxation protocols reduce clashscore
- Best protocol: cartesian_beta (largest improvement)
- Diminishing returns across replicates

### 4. Anomalies Detected
- 4 structures with zero clashscore (all AlphaFold)
- 9 structures with extremely high clashscore (>100)
- Some structures show inconsistent quality profiles (good backbone + high clashes)

---

## Interpretation Guide

### Clashscore
- **< 10:** Excellent (better than average X-ray at 2.0Å)
- **10-20:** Good
- **20-40:** Acceptable
- **> 40:** Poor (significant steric problems)

### Ramachandran Favored %
- **> 98%:** Excellent
- **95-98%:** Good
- **90-95%:** Acceptable
- **< 90%:** Poor

### MolProbity Score
- **< 1.5:** Excellent
- **1.5-2.0:** Good
- **2.0-2.5:** Acceptable
- **> 2.5:** Poor

### Bond/Angle RMSZ
- **≈ 1.0:** Ideal (matches restraints perfectly)
- **< 1.5:** Good
- **1.5-2.5:** Acceptable
- **> 2.5:** Significant deviations

---

## Statistical Methods Used

1. **Descriptive Statistics:** Mean, SD, median, IQR, skewness, kurtosis
2. **Normality Tests:** Shapiro-Wilk
3. **Group Comparisons:** Kruskal-Wallis H-test, Mann-Whitney U-test
4. **Correlations:** Pearson (linear), Spearman (monotonic), Kendall (ordinal)
5. **Effect Sizes:** Cohen's d, Glass's delta, Cliff's delta
6. **Outlier Detection:** IQR method, Z-score method

---

*Vanderbilt Protein Structure Validation Project*
