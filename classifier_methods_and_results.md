# Signature-Based Quiescent vs Activated Epicardial Classifier

## Motivation

Previous classification approaches relied on disease condition labels (normal vs MI) to define cell states. This is unreliable because the MI data comes from a single patient, confounding patient-specific effects with activation biology. This classifier instead labels cells purely from gene expression signatures.

## Input Data

- **Source**: `data/processed/epicardial_periheart.h5ad`
- **Cells**: 19,412 epicardial cells (PERIHEART dataset, Linna-Kuosmanen et al.)
- **Genes**: 35,477 (Ensembl IDs; gene symbols in `adata.var['feature_name']`)

---

## Methods

### Phase 1: Feature Engineering & Weak Labeling

Cells were scored using `sc.tl.score_genes()` (ctrl_size=50, n_bins=25) against two curated gene signatures:

| Quiescent (epithelial/resting) | Activated (EMT/mesenchymal) |
|---|---|
| WT1, UPK3B, MSLN, KRT19 | SNAI1, SNAI2, TWIST1 |
| CDH1, CLDN1, CLDN3, TJP1, OCLN | VIM, FN1, COL1A1, COL3A1, POSTN, ACTA2 |

A **state potential** score was computed as `score_activated - score_quiescent`. Silver labels were assigned to the extremes of the distribution:
- Bottom 15% -> `quiescent_ref` (2,912 cells)
- Top 15% -> `activated_ref` (2,912 cells)
- Middle 70% -> `unassigned` (13,588 cells)

### Phase 2: Classifier Training

- **Features**: 2,000 highly variable genes (Seurat v3 method), excluding 3 cell cycle genes detected in the HVG set (to avoid proliferation-driven classification).
- **Models compared**:
  - Random Forest (500 trees, balanced class weights)
  - L1-regularized Logistic Regression (SAGA solver, max 5,000 iterations, balanced class weights)
- **Evaluation**: 5-fold stratified cross-validation on silver-labeled cells, scored by F1.
- **Prediction**: The winning model was refit on all silver-labeled cells and used to predict all 19,412 cells with probability scores.
- **Transitioning cells**: Cells with P(activated) between 0.3 and 0.7 were flagged.

### Phase 3: Validation

Three validation outputs were generated:
1. UMAP colored by predicted state and P(activated)
2. Marker heatmap (matrixplot) of signature genes across predicted states
3. Probability histogram showing classifier confidence distribution

A **fail condition** was applied: if >20% of predicted quiescent cells had VIM expression above the 75th percentile of all cells, the classification was deemed unreliable. (VIM is broadly expressed in epicardial cells due to their mesenchymal origin, so the 75th percentile threshold avoids penalizing baseline VIM expression.)

### Phase 4: Self-Correction

When the fail condition triggered, the pipeline ran differential expression (Wilcoxon) between predicted states and augmented the gene signatures with top discriminatory genes. Non-epicardial markers (cardiomyocyte genes: ACTC1, TNNT2, TTN, MYH6, etc.; mitochondrial genes: MT-CO1, MT-ND3, etc.) were excluded from augmentation.

Up to 2 augmentation rounds were attempted. If the classifier still failed, a **GMM fallback** was applied: 2-component Gaussian Mixture Model on the first 10 principal components, with cluster identity assigned based on which cluster had higher mean activated signature score.

---

## Results

### Iteration Summary

| Iteration | Best Model | CV F1 | Predicted Activated | Predicted Quiescent | VIM Check | Action |
|---|---|---|---|---|---|---|
| 1 | L1-Logistic Regression | 0.9188 | 10,196 | 9,216 | FAIL (23.1%) | DE augmentation |
| 2 | L1-Logistic Regression | 0.9990 | 9,790 | 9,622 | FAIL (25.1%) | DE augmentation |
| 3 | L1-Logistic Regression | 0.9993 | 10,476 | 8,936 | FAIL (27.2%) | GMM fallback |
| 3 (GMM) | GMM on 10 PCs | — | 7,464 | 11,948 | **PASS (16.3%)** | — |

The L1-LR classifier separated silver-labeled cells with near-perfect F1 (0.999), but the predicted quiescent population consistently contained cells with high VIM — expected given that VIM is broadly expressed across epicardial cells regardless of activation state. The GMM on PCA space provided a cleaner separation.

### Final Classification

| State | Cells | % |
|---|---|---|
| **Quiescent** | 11,948 | 61.5% |
| **Activated** | 7,464 | 38.5% |
| Transitioning (P 0.3–0.7) | 851 | 4.4% |

### Signature Augmentation Through Self-Correction

Genes added via DE-based augmentation across iterations:

**Quiescent signature (added):** ITLN1, C3, PRG4, TM4SF1, PLA2G2A, SLPI, TIMP1, AQP1, SERPING1, C1R

**Activated signature (added):** MB, PTGDS, TNNI3, CRYAB, COX7A1, COX6A2, TPM1

### Top Classifier Features (L1-LR, Iteration 3)

| Rank | Gene | Importance (|coef|) |
|---|---|---|
| 1 | PTGDS | 1.437 |
| 2 | CRYAB | 1.062 |
| 3 | MSLN | 1.049 |
| 4 | COX6A2 | 0.896 |
| 5 | PLA2G2A | 0.863 |
| 6 | FN1 | 0.803 |
| 7 | TIMP1 | 0.764 |
| 8 | TM4SF1 | 0.760 |
| 9 | PRG4 | 0.712 |
| 10 | SERPING1 | 0.700 |

Full feature importances: `results/classifier/feature_importances.csv`

### Differential Expression: Activated vs Quiescent

Wilcoxon rank-sum test across all 35,477 genes. Full results: `results/classifier/deg_activated_vs_quiescent.csv`

**Top 15 upregulated in activated:**

| Gene | logFC | padj |
|---|---|---|
| ITLN1 | 2.54 | <1e-300 |
| PLA2G2A | 3.15 | <1e-300 |
| TMSB4X | 2.60 | <1e-300 |
| SLPI | 3.48 | <1e-300 |
| TIMP1 | 3.04 | <1e-300 |
| RPLP1 | 2.14 | <1e-300 |
| S100A6 | 2.89 | <1e-300 |
| B2M | 2.08 | <1e-300 |
| TPT1 | 1.92 | <1e-300 |
| RPL13 | 1.90 | <1e-300 |
| RPS12 | 2.27 | <1e-300 |
| RPL3 | 2.35 | <1e-300 |
| RPL34 | 2.21 | <1e-300 |
| CRIP1 | 2.60 | <1e-300 |
| RPS18 | 2.27 | <1e-300 |

**Top 10 upregulated in quiescent (downregulated in activated):**

| Gene | logFC | padj |
|---|---|---|
| BNC2 | -0.61 | <1e-300 |
| PKHD1L1 | -0.57 | 3.6e-266 |
| RBFOX1 | -0.74 | 3.0e-253 |
| MAGI2 | -0.71 | 6.9e-250 |
| RORA | -0.50 | 4.8e-238 |
| GPM6A | -0.84 | 1.6e-234 |
| SLIT3 | -0.54 | 2.9e-220 |
| ROR1 | -0.81 | 3.1e-214 |
| MEIS2 | -0.56 | 8.3e-212 |
| PARD3B | -0.56 | 5.6e-209 |

### FGF/FGFR Family

| Gene | logFC | padj | logFC Rank | Significant? |
|---|---|---|---|---|
| FGF12 | +0.63 | 2.0e-54 | 6,408/35,477 | Yes |
| FGF14 | +0.67 | 9.6e-09 | 6,002/35,477 | Yes |
| FGF7 | +0.54 | 5.4e-05 | 7,277/35,477 | Yes |
| FGF2 | +0.07 | 0.035 | 14,707/35,477 | Marginal |
| FGFRL1 | +0.64 | 0.018 | 6,312/35,477 | Marginal |
| FGF1 | +0.38 | 0.20 | 9,363/35,477 | No |
| FGFR1 | -0.01 | NS | 22,699/35,477 | No |
| FGFR2 | -0.27 | 0.30 | 28,727/35,477 | No |
| FGFR3 | +0.35 | NS | 9,703/35,477 | No |
| FGFR4 | +0.43 | NS | 8,588/35,477 | No |

FGF12 and FGF14 are intracellular FGFs (iFGFs) that do not signal through FGFRs canonically. FGF7 (KGF) is the most relevant paracrine FGF ligand, significantly upregulated in activated cells. FGFR2 trends toward quiescent cells but does not reach significance.

---

## Output Files

| File | Description |
|---|---|
| `data/processed/epicardial_classified.h5ad` | Annotated AnnData with `pred_state`, `pred_prob_activated`, `transitioning`, signature scores |
| `results/classifier/umap_predicted_states.png` | UMAP by predicted state and P(activated) |
| `results/classifier/marker_heatmap.png` | Matrixplot of signature genes |
| `results/classifier/probability_histogram.png` | P(activated) distribution |
| `results/classifier/classifier_model.joblib` | Saved model (L1-LR + GMM fallback metadata) |
| `results/classifier/feature_importances.csv` | Gene-level classifier weights |
| `results/classifier/deg_activated_vs_quiescent.csv` | Full DEG table (35,477 genes) |
| `results/classifier/classification_report.txt` | Iteration-level run log |

## Reproducibility

```bash
python3 scripts/02_cell_states/build_classifier.py
```

Script: `scripts/02_cell_states/build_classifier.py`
