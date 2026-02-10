# Wednesday Meeting Prep — BM5.5 Relaxation Benchmark

## 1. Project Overview

**Goal:** Systematic evaluation of Rosetta relaxation protocols on the BM5.5 protein-protein docking benchmark using AlphaFold-Multimer and Boltz-1 predictions.

**Design:**
- Two independent pipelines (Blue + Green) for reproducibility verification
- 257 BM5.5 targets (253 standard + 4 non-standard alternate chain combinations)
- 7 relaxation protocols tested:
  - 6 Rosetta protocols: normal_beta, normal_ref15, cartesian_beta, cartesian_ref15, dualspace_beta, dualspace_ref15
  - 1 built-in: AMBER (AlphaFold's default relaxation)
- 5 replicates per Rosetta protocol
- Expected output: ~81,000 structures per pipeline

**Prediction Sources:**
- AlphaFold-Multimer v2.3 (full MSA database, --models_to_relax=all)
- Boltz-1 (5 diffusion samples per target)
- Experimental crystal structures (BM5.5 bound/unbound pairs)

**Validation Metrics:**
- MolProbity: clashscore, Ramachandran outliers, rotamer outliers, bond/angle RMSZ
- DockQ: interface accuracy relative to experimental bound structure (planned)

---

## 2. Key Findings (Pilot Data — 20 proteins, 6,820 structures)

### 2.1 The Clashscore ~14 Attractor (r = -0.997)

All structures converge to clashscore ~14 after Rosetta relaxation, regardless of starting point:

| Source | Initial Clashscore | Post-Relaxation |
|--------|-------------------|-----------------|
| Boltz | 38.7 | 13.7 |
| AF ranked_1-4 (unrelaxed) | 22.8 | 14.1 |
| Experimental | 15.7 | 14.0 |
| AF ranked_0 (AMBER-relaxed) | 0.93 | 13.8 |

**Interpretation:** ~14 is the Rosetta force field equilibrium. Structures above improve; structures below get worse. This is regression to force field minimum, not "fixing" structures.

### 2.2 Split AF Analysis — The Hidden Story

Initial pooled analysis showed AF improvement was NOT statistically significant (p=0.28). Split analysis revealed:

| AF Model | Mean Improvement | p-value | Direction |
|----------|-----------------|---------|-----------|
| ranked_0 (AMBER-relaxed) | -0.83 | 0.000002 | **WORSE** |
| ranked_1-4 (unrelaxed) | +0.33 | 0.000322 | **BETTER** |

**Key insight:** Rosetta relaxation helps unrefined predictions but harms already-refined structures. The p=0.28 was an artifact of averaging opposite effects.

### 2.3 Protocol Ranking

| Protocol | Composite Score |
|----------|-----------------|
| normal_beta | 0.95 |
| normal_ref15 | 0.93 |
| dualspace_beta | 0.77 |
| cartesian_ref15 | 0.73 |
| raw (unrelaxed) | 0.72 |

**Result:** Normal-mode relaxation dominates. Cartesian/dualspace underperform — aggressive minimization may overshoot for structures already near energy minima.

### 2.4 Ceiling Effect

22 structures degraded after relaxation. Characterization:
- NOT all AlphaFold (distributed across categories)
- Initial MolProbity: 0.99 (degraded) vs 1.66 (improved)
- Initial clashscore: 9.9 (degraded) vs 32.7 (improved)

**Conclusion:** Degraded structures started too good (below ~14 attractor). Relaxation pushed them toward equilibrium = got worse.

---

## 3. Pipeline Status

| Component | Blue | Green |
|-----------|------|-------|
| AF progress | 43/257 (17%) | 29/257 (11%) |
| Boltz progress | pending | 248/257 (96%) |
| Rosetta relaxation | pending (after AF) | pending (after AF) |
| Input verification | Clean (10 FASTAs fixed) | Clean (DNA stripped) |
| ETA to AF complete | ~12h | ~24-36h |

**Notes:**
- Blue: 1ATN AMBER failure (unrelaxed only), 10 FASTAs regenerated from RCSB
- Green: Boltz 96% complete, 9 targets permanently OOM (>3000 residues, AF-only)
- Both pipelines use --models_to_relax=all (enables paired AMBER analysis)
- 203 targets have reversed chain order between pipelines (robustness test)

---

## 4. Paper Outline

**Working Title:** "Systematic evaluation of relaxation protocols for AI-predicted protein complex structures on the BM5.5 docking benchmark"

### Planned Analyses:
1. **Relaxation delta** — improvement by protocol and source
2. **Protocol ranking** — which Rosetta protocol is best?
3. **Convergence test** — do all sources converge to same geometry?
4. **AMBER vs Rosetta** — paired analysis (same AF model, ±AMBER, then Rosetta)
5. **MSA depth effect** — full_dbs vs reduced_dbs fallback targets
6. **DockQ correlation** — does clashscore improvement correlate with docking accuracy?

### Key Questions:
1. Does relaxation improve or harm AI predictions?
2. Is the improvement cosmetic (MolProbity) or functional (DockQ)?
3. Which protocol should practitioners use?

### Target Journal: TBD

---

## 5. Figures (Pilot Data)

### Figure A: Clashscore Trajectory to ~14
![Clashscore Correlation](figures/clashscore_correlation.png)

All sources converge to ~14 from opposite directions. r = -0.997 correlation between initial clashscore and change.

### Figure B: Split AF Analysis
![Split AF Analysis](figures/split_af_analysis.png)

ranked_0 (AMBER-relaxed) gets worse; ranked_1-4 (unrelaxed) get better. Both converge to ~14.

### Figure C: Protocol Ranking Heatmap
![Protocol Ranking](figures/fig2_protocol_ranking.png)

normal_beta dominates across metrics.

### Figure D: Outlier Characterization
![Outlier Analysis](figures/outlier_characterization.png)

Ceiling effect confirmed — degraded structures started with better scores.

---

## 6. Quality Tiers for Analysis

### Boltz Tiers:
1. **Standard** (245 targets): 5 samples on L40S/H100
2. **Reduced** (2 targets: 1GXD, 3EO1): 1 sample on H100
3. **No Boltz** (9 targets): >3000 residues, OOM even on H100 80GB

### AF Tiers:
1. **Standard**: full_dbs MSA
2. **Fallback**: reduced_dbs MSA (HHblits titin issue)
3. **Partial**: 1ATN (AMBER failed, unrelaxed only)

---

## 7. Known Limitations

1. **9 largest complexes have no Boltz predictions** — systematic bias toward smaller structures in Boltz analysis
2. **DNA/RNA excluded** — some BM5.5 targets have nucleic acid-mediated interfaces (e.g., 3P57). DockQ may be reduced for these.
3. **Chain order differs between pipelines** — Blue uses RCSB order, Green uses BM5.5 receptor-first. Should not affect global metrics.

---

## 8. Next Steps

1. Complete AF runs (Blue ~12h, Green ~24-36h)
2. Run Rosetta relaxation (6 protocols × 5 replicates × ~81k structures)
3. Collect MolProbity metrics
4. Run DockQ against bound experimental structures
5. Statistical analysis and figure generation
6. Write manuscript

---

*Last updated: 2026-02-10 ~02:00 PST*
*Prepared by: Teal (Analysis Pipeline)*
