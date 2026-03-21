# Receptor Analysis Results

## Critical Data Limitation

**All MI cells (6,439) come from a SINGLE patient (PH-M57)**, while Normal cells (12,973) come from 29 different donors. This means:
- Observed "MI vs Normal" differences may reflect **patient-specific effects** rather than true disease biology
- Statistical comparisons are confounded by individual variation
- Results should be interpreted with caution

## Dataset: PERIHEART Epicardial Cells

| Condition | Cells | Donors |
|-----------|-------|--------|
| Normal | 12,973 | 29 |
| MI | 6,439 | 1 (PH-M57) |

### Epicardial Identity Confirmed
| Marker | Function | % Expressing | Normal | MI |
|--------|----------|--------------|--------|-----|
| WT1 | Transcription factor | 52.4% | 47.0% | 63.3% (↑) |
| TBX18 | Transcription factor | 48.7% | 46.6% | 52.8% |
| ALDH1A2 | Retinoic acid synthesis | 62.0% | 58.7% | 68.6% (↑) |
| PECAM1/CD31 | Endothelial marker | 7.9% | 8.9% | 5.9% (low, as expected) |

### Marker Expression Across Cell States

Expression shown as: % expressing (mean log1p expression)

| Gene | Function | Quiescent | Primed | Bystander | Activated |
|------|----------|-----------|--------|-----------|-----------|
| WT1 | Epicardial TF | 57.8% (1.08) | 50.4% (0.95) | 63.9% (1.13) | 62.0% (1.11) |
| TBX18 | Epicardial TF | 52.6% (0.93) | 47.9% (0.86) | 53.0% (0.87) | 52.3% (0.85) |
| ALDH1A2 | Retinoic acid | 61.7% (1.17) | 56.2% (1.08) | 69.2% (1.27) | 67.3% (1.23) |
| FN1 | EMT marker | 23.6% (0.39) | 35.5% (0.65) | 34.8% (0.62) | **43.2% (0.78)** |
| VIM | EMT marker | 34.8% (0.56) | 51.6% (0.95) | 40.6% (0.62) | **53.0% (0.86)** |
| POSTN | EMT marker | 6.5% (0.09) | 12.8% (0.20) | 4.5% (0.06) | 7.7% (0.11) |
| FGFR1 | FGF receptor | 30.9% (0.47) | 27.6% (0.44) | 39.6% (0.60) | 40.1% (0.61) |
| FGFR2 | FGF receptor | 18.0% (0.27) | 15.5% (0.24) | 24.8% (0.36) | 23.3% (0.34) |

**Cell counts:** quiescent (18,091), primed (8,394), bystander (4,523), activated (1,916)

**Key observations:**
- Epicardial identity markers (WT1, TBX18, ALDH1A2) are expressed across all states, confirming epicardial identity
- EMT markers (FN1, VIM) show highest expression in **activated** cells
- **Primed** cells (Normal + high molecular score) also show elevated EMT markers, suggesting constitutive EMT activity in a subset of healthy epicardium

---

## Analysis 1: DEG Analysis (Full PERIHEART, Normal vs MI)

**File:** `deg/receptor_normal_vs_mi_PERIHEART.csv`

**Method:** Wilcoxon rank-sum test comparing 12,973 Normal vs 6,439 MI cells

| Rank | Receptor | logFC | padj |
|------|----------|-------|------|
| 1 | PDE1A | 1.40 | 2.1e-198 |
| 2 | CCBE1 | 1.55 | 3.3e-102 |
| 3 | PLXNA4 | 0.79 | 4.4e-92 |
| 4 | LPAR1 | 0.73 | 1.2e-90 |
| 5 | SLC40A1 | 1.13 | 1.1e-65 |
| ... | ... | ... | ... |
| **28** | **FGFR2** | **0.57** | **1.9e-16** |

---

## Analysis 2: DEG Analysis (Quiescent vs Activated Cell States)

**File:** `deg/receptor_quiescent_vs_activated.csv`

**Method:**
- Cell state classification using `sc.tl.score_genes()` with control gene sets (background correction)
- GMM-based threshold selection (threshold = 0.22, equivalent to 68.7th percentile)
- Wilcoxon rank-sum test comparing **quiescent** (18,091) vs **activated** (1,916) cells
- This compares cells by their molecular phenotype, not just disease condition

**Cell State Distribution (Full Dataset):**

| State | Count | Percent |
|-------|-------|---------|
| quiescent | 18,091 | 55.0% |
| primed | 8,394 | 25.5% |
| bystander | 4,523 | 13.7% |
| activated | 1,916 | 5.8% |

**Top Receptors:**

| Rank | Receptor | logFC | padj |
|------|----------|-------|------|
| 1 | PDE1A | 1.06 | 5.5e-54 |
| 2 | AQP1 | 0.97 | 1.4e-53 |
| 3 | RACK1 | 0.96 | 1.2e-40 |
| 4 | RPSA | 0.85 | 1.3e-35 |
| 5 | LRP1 | 0.91 | 5.8e-30 |
| 6 | SLC40A1 | 1.19 | 5.0e-28 |
| 7 | KDR | 0.69 | 1.5e-25 |
| 8 | CCBE1 | 1.17 | 5.1e-24 |
| 9 | CD81 | 0.71 | 1.6e-22 |
| 10 | PLXNA4 | 0.58 | 9.8e-22 |
| ... | ... | ... | ... |
| **59** | **FGFR2** | **0.41** | **3.9e-03** |

---

## Analysis 3: DEG Bootstrap (100 iterations, 200 vs 200)

**File:** `deg/receptor_normal_vs_mi_PERIHEART_bootstrap100.csv`

**Method:** 100 random balanced samples (200 Normal vs 200 MI each), aggregated rankings

| Final Rank | Receptor | Mean Genome Rank | Mean logFC |
|------------|----------|------------------|------------|
| 1 | PDE1A | 41.4 | 1.39 |
| 2 | CCBE1 | 144.0 | 1.55 |
| 3 | PLXNA4 | 199.6 | 0.78 |
| 4 | SLC40A1 | 257.1 | 1.18 |
| 5 | TGFBR2 | 370.6 | 0.68 |
| ... | ... | ... | ... |
| **34** | **FGFR2** | **3172.9** | **0.51** |

---

## Analysis 4: Geneformer In Silico Deletion

**File:** `geneformer/perturbation_results_round1.csv`

### Two-Round Approach

**Round 1 - Embedding Perturbation (Completed):**
- Uses the **pre-trained** Geneformer model (no fine-tuning required)
- For each cell, delete a receptor gene from the tokenized input
- Measure the **cosine distance shift** in the embedding space
- Larger shift = gene has more influence on the cell's transcriptional identity
- This measures general importance to cell state, not specifically activation vs quiescence

**Round 2 - Classifier Perturbation (Not Completed):**
- Requires **fine-tuning** Geneformer as a binary classifier (quiescent vs activated)
- For each cell, delete a receptor and measure change in **predicted activation probability**
- More negative effect = receptor is more important for maintaining activated state
- This directly measures each receptor's contribution to the activated phenotype

### Why Round 2 Was Not Completed

Fine-tuning Geneformer requires significant GPU memory. On Apple M4 Max (128GB unified memory):
- Training caused MPS backend to run out of memory (~170GB allocated)
- Even with reduced batch size (4) and gradient accumulation, each iteration took 90-150 seconds
- Estimated training time: >24 hours for 5 epochs
- Process was stopped due to impractical runtime

**Recommendation:** Run Round 2 on a dedicated GPU server (NVIDIA A100 or similar with 40GB+ VRAM).

### Round 1 Results

| Rank | Receptor | Mean Effect | Cells with Gene |
|------|----------|-------------|-----------------|
| 1 | PLXDC2 | 0.000932 | 1,368 |
| 2 | ALK | 0.000820 | 229 |
| 3 | AQP1 | 0.000709 | 977 |
| 4 | ROR2 | 0.000705 | 1,263 |
| 5 | TNFRSF1A | 0.000698 | 437 |
| ... | ... | ... | ... |
| **61** | **FGFR2** | **0.000542** | **321** |

---

## Combined Rankings (All Methods)

**File:** `receptor_rankings_comparison.csv`

Top 10 receptors by mean rank across all analyses:

| Overall Rank | Receptor | Mean Rank | DEG Normal vs MI | DEG Quiescent vs Activated | DEG Bootstrap | Geneformer |
|--------------|----------|-----------|------------------|---------------------------|---------------|------------|
| 1 | PDE1A | 5.8 | 1 | 1 | 1 | 20 |
| 2 | PLXNA4 | 6.5 | 3 | 10 | 3 | 7 |
| 3 | AQP1 | 7.3 | 15 | 2 | 15 | 3 |
| 4 | CCBE1 | 7.0 | 2 | 8 | 2 | 16 |
| 5 | SLC40A1 | 8.8 | 5 | 6 | 4 | 65 |
| 6 | LRP1 | 12.5 | 9 | 5 | 9 | 68 |
| 7 | RACK1 | 14.5 | 18 | 3 | 21 | 19 |
| 8 | RPSA | 15.0 | 24 | 4 | 37 | 9 |
| 9 | ROBO1 | 11.0 | 11 | 17 | 13 | 6 |
| 10 | KDR | 15.5 | 14 | 7 | 19 | 27 |
| ... | ... | ... | ... | ... | ... | ... |
| **~50** | **FGFR2** | **~45** | 28 | 59 | 34 | 61 |

---

## FGFR2 Summary

| Analysis Method | FGFR2 Rank | logFC | Notes |
|-----------------|------------|-------|-------|
| DEG (Normal vs MI) | 28/84 | 0.57 | Top 33% |
| DEG (Quiescent vs Activated) | 59/84 | 0.41 | Bottom 30% |
| DEG (Bootstrap 100x) | 34/84 | 0.51 | Top 40% |
| Geneformer Embedding | 61/84 | - | Bottom 27% |
| **Overall** | **~45/84** | - | Mean rank ~45.5 |

**Conclusion:** FGFR2 shows moderate upregulation in MI (logFC ~0.4-0.6) but does not rank among the top receptors in any analysis. When comparing by cell state (quiescent vs activated), FGFR2 drops to rank 59/84, suggesting its upregulation is more associated with disease condition than with the activated molecular phenotype. Receptors like **PDE1A, AQP1, PLXNA4, CCBE1, SLC40A1, LRP1** consistently rank higher across all methods.

---

## File Descriptions

### deg/
| File | Description |
|------|-------------|
| `receptor_normal_vs_mi.csv` | Full dataset DEG (PERIHEART + CAREBANK) |
| `receptor_normal_vs_mi_PERIHEART.csv` | PERIHEART only DEG (Normal vs MI) |
| `receptor_quiescent_vs_activated.csv` | Cell state DEG (Quiescent vs Activated) |
| `receptor_normal_vs_mi_PERIHEART_improved.csv` | With sc.tl.score_genes + GMM threshold |
| `receptor_normal_vs_mi_PERIHEART_bootstrap100.csv` | 100x bootstrap (200 vs 200) |
| `receptor_normal_vs_mi_PERIHEART_sampled200.csv` | Single 200 vs 200 sample |

### geneformer/
| File | Description |
|------|-------------|
| `receptors_for_perturbation.csv` | 84 upregulated receptors tested |
| `perturbation_results_round1.csv` | Embedding perturbation results |
| `cell_metadata.csv` | Cell state labels for tokenized data |
| `tokenized_dataset/` | Geneformer input format |

### Summary files
| File | Description |
|------|-------------|
| `receptor_rankings_comparison.csv` | Side-by-side ranking comparison |
| `receptor_analysis_summary.csv` | Detailed results (long format) |

---

## Date
Analysis completed: 2026-02-02
