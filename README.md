# AI-Driven Epicardial Target Discovery in Myocardial Infarction

> **Purpose**: Identify potential therapeutic targets on epicardial cells that promote proliferation and EMT for cardioprotection after myocardial infarction

**Version**: 2.3 (January 2026)

---

## Table of Contents

1. [Overview and Rationale](#1-overview-and-rationale)
2. [Current Progress](#2-current-progress)
3. [Phase 1: Data Acquisition and Preprocessing](#3-phase-1-data-acquisition-and-preprocessing)
4. [Phase 2: Define Target Phenotypes](#4-phase-2-define-target-phenotypes)
5. [Phase 3: Ligand-Receptor Analysis](#5-phase-3-ligand-receptor-analysis)
6. [Phase 4: Machine Learning Model (Geneformer)](#6-phase-4-machine-learning-model-geneformer)
7. [Phase 5: Model Interpretation and Prioritization](#7-phase-5-model-interpretation-and-prioritization)
8. [Tools and Resources](#8-tools-and-resources)
9. [Key References](#9-key-references)

---

## 1. Overview and Rationale

### Goal

Identify receptors on epicardial cells whose activation leads to:
- Epicardial cell proliferation
- Epithelial-to-mesenchymal transition (EMT)
- Cardioprotective paracrine secretion

### Biological Context

- After MI, epicardial cells can be activated to support cardiac repair
- FGF10 delivered via pericardial cavity activates epicardium through FGFR2b
- Epicardium-derived cells (EPDCs) contribute to repair through:
  - Paracrine factor secretion (VEGF, PDGF, FGFs)
  - Differentiation into cardiac fibroblasts and smooth muscle
  - Support of revascularization

### Sender/Receiver Framework

| Role | Cell Types |
|------|------------|
| **Sender cells** | Cardiomyocytes, macrophages, fibroblasts, endothelial cells |
| **Receiver cells** | Epicardial cells (target population) |

**Positive Control**: FGF10/FGFR2b axis (should rank highly if model is working)

---

## 2. Current Progress

### Project Structure

```
epicardial-target-discovery/
├── README.md
├── data/
│   ├── raw/
│   │   ├── kuppe/
│   │   │   └── 54d24dbe-*.h5ad          # Kuppe MI data (1.8GB, 191,795 cells)
│   │   └── linna_kuosmanen/
│   │       ├── 77df1b04-*.h5ad          # PERIHEART (2.6GB, 392,819 cells)
│   │       └── 8e1afde5-*.h5ad          # CAREBANK (1.5GB, 221,756 cells)
│   ├── processed/
│   │   ├── epicardial_all_merged.h5ad     # 59,431 cells (merged, log-normalized)
│   │   ├── epicardial_with_states.h5ad    # 59,431 cells (with cell state labels)
│   │   ├── epicardial_periheart_log.h5ad  # 28,851 cells (log-normalized)
│   │   ├── epicardial_carebank_log.h5ad   # 20,493 cells (log-normalized)
│   │   ├── epicardial_kuppe.h5ad          # 10,087 cells (log-normalized)
│   │   ├── kuppe_cell_metadata.csv
│   │   └── kuppe_gene_metadata.csv
│   └── references/
├── scripts/
│   ├── 01_preprocessing/
│   │   ├── analyze_kuppe.py
│   │   ├── extract_epicardial_simple.py
│   │   ├── extract_kuppe_epicardial.py
│   │   ├── compare_epicardial_scores.py
│   │   ├── normalize_and_compare.py
│   │   └── merge_epicardial_datasets.py
│   └── 02_cell_states/
│       └── classify_cell_states.py
└── papers/
```

### Completed

- [x] Download Kuppe et al. 2022 from CellxGene (MI heart, left ventricle)
- [x] Download Linna-Kuosmanen et al. 2024 from CellxGene (PERIHEART + CAREBANK)
- [x] Explore data structure (h5ad format, AnnData)
- [x] Identify epicardial cell populations:
  - Kuppe: ~523 adipocytes (epicardial fat); true epicardial may be in Fibroblast
  - Linna-Kuosmanen: **~49,000 mesothelial cells** (7-9% of total) - best source!
- [x] Verify epicardial markers present: WT1, TBX18, ALDH1A2, UPK3B, MSLN
- [x] Set up Python-based workflow (scanpy/anndata)
- [x] Extract epicardial cells from all datasets:
  - Linna-Kuosmanen: filtered by `cell_type == 'mesothelial cell'` (49,344 cells)
  - Kuppe: filtered by marker score (95th percentile) + adipocyte annotation (10,087 cells)
- [x] Harmonize normalization across datasets:
  - Discovered: Linna-Kuosmanen uses CPM, Kuppe uses log1p(CPM)
  - Applied log1p to Linna-Kuosmanen for consistency
  - Verified score comparability (ratio 1.40x)
- [x] Merge all datasets into `epicardial_all_merged.h5ad` (59,431 cells, 28,380 genes)
- [x] Phase 2: Define cell states with dual-method validation:
  - Method 1: Condition-based (disease/spatial zone annotations)
  - Method 2: Score-based (EMT + proliferation signatures)
  - Cross-validation agreement: 66.8%
  - Output: `epicardial_with_states.h5ad`

### Next Steps
- [ ] Run NicheNet/LIANA with sender cells → epicardium
- [ ] Verify FGF10 ranks in top ligands (positive control)
- [ ] Fine-tune Geneformer for quiescent vs activated classification
- [ ] Run in silico perturbation on candidate receptors
- [ ] Verify FGFR2 shows strong activation effect
- [ ] Prioritize novel receptor candidates

---

## 3. Phase 1: Data Acquisition and Preprocessing

### 3.1 Data Sources

**Data Portal**: https://cellxgene.cziscience.com/gene-expression

#### Kuppe et al. 2022 (Nature)

| Property | Value |
|----------|-------|
| **File** | `data/raw/kuppe/54d24dbe-d39a-4844-bb21-07b5f4e173ad.h5ad` |
| **Size** | 1.8 GB |
| **Tissue** | Heart left ventricle |
| **Cells** | 191,795 |
| **Genes** | 28,380 |
| **Conditions** | MI (ischemic, border, remote zones) vs normal |
| **Paper** | https://www.nature.com/articles/s41586-022-05060-x |

**Key findings:**
- Disease states: myocardial infarction (168,044 cells), normal (23,751 cells)
- Patient groups: myogenic, ischemic, fibrotic
- Epicardial markers present: WT1, TBX18, ALDH1A2, UPK3B (MSLN not annotated)
- Only ~523 "adipocyte" cells labeled; true epicardial cells likely mixed in Fibroblast population

#### Linna-Kuosmanen et al. 2024 (Cell Reports Medicine)

Two cohorts from Finnish cardiac surgery patients (right atrium tissue):

| Dataset | File | Cells | Genes | Mesothelial (Epicardial) |
|---------|------|-------|-------|--------------------------|
| **PERIHEART** | `77df1b04-*.h5ad` | 392,819 | 35,477 | 28,851 (7.3%) |
| **CAREBANK** | `8e1afde5-*.h5ad` | 221,756 | 35,477 | 20,493 (9.2%) |

**Key findings:**
- **Rich epicardial annotations**: "mesothelial cell" = epicardial cells (~49,000 total!)
- All 6 epicardial markers present including MSLN
- Multiple cell types with clear annotations
- Best source for epicardial cell analysis

### 3.2 Standard scRNA-seq Processing

**Quality Control:**
- Filter cells: 200 < nFeature_RNA < 6000
- Mitochondrial content: < 10-20% (tissue-dependent)
- Doublet removal: Scrublet

**Normalization:**
- Log-normalization: `log(raw_count / total_counts × 10000 + 1)`
- Or SCTransform equivalent in scanpy

**Batch correction (if integrating datasets):**
- Harmony
- scVI
- BBKNN

### 3.3 Epicardial Cell Identification

**Canonical epicardial markers:**
| Marker | Ensembl ID | Function |
|--------|------------|----------|
| WT1 | ENSG00000184937 | Master regulator |
| TBX18 | ENSG00000112837 | Transcription factor |
| ALDH1A2 | ENSG00000128918 | Retinoic acid synthesis |
| UPK3B | ENSG00000243566 | Mesothelial marker |
| MSLN | ENSG00000102854 | Mesothelin |

> **Note**: TCF21 excluded due to low expression in adult human epicardium

**Activated/EMT epicardial markers:**
- POSTN (Periostin)
- COL1A1, COL3A1 (Collagens)
- ACTA2 (alpha-SMA)
- VIM (Vimentin)

**Cell selection methods:**
- Linna-Kuosmanen: `cell_type == 'mesothelial cell'`
- Kuppe: `epicardial_score > 95th percentile` (using sc.tl.score_genes)

### 3.4 Normalization Harmonization

**Issue Discovered:** The two data sources use different normalization methods:

| Dataset | Original Format | Expression Range | Notes |
|---------|----------------|------------------|-------|
| **Linna-Kuosmanen** | CPM/TPM (linear) | 0 - 1353 | No log transformation |
| **Kuppe** | log1p(CPM) | 0 - 7.7 | `X_approximate_distribution: normal` |

**Verification:** Applying `log1p` to Linna-Kuosmanen data produces range (0 - 7.2), matching Kuppe.

**Solution:** Applied `log1p` transformation to Linna-Kuosmanen datasets for consistency.

**Epicardial Score Comparison After Normalization:**

| Dataset | Cells | Mean Score | Median | Std | 95th % |
|---------|-------|------------|--------|-----|--------|
| PERIHEART | 28,851 | 0.8664 | 0.8763 | 0.5003 | 1.6881 |
| CAREBANK | 20,493 | 0.9254 | 0.9537 | 0.4983 | 1.7108 |
| Kuppe_MI | 10,087 | 0.6590 | 0.6287 | 0.2071 | 1.0155 |

**Individual Marker Expression (Mean):**

| Dataset | WT1 | TBX18 | ALDH1A2 | UPK3B |
|---------|-----|-------|---------|-------|
| PERIHEART | 0.95 | 0.83 | 1.12 | 0.57 |
| CAREBANK | 1.18 | 1.00 | 1.16 | 0.36 |
| Kuppe_MI | 0.34 | 0.89 | 1.39 | 0.01 |

**Interpretation:**
- Max/Min score ratio = **1.40x** (acceptable for combined analysis)
- Linna-Kuosmanen datasets show higher WT1/UPK3B (true epicardial annotations)
- Kuppe shows higher ALDH1A2 (marker used for selection)
- UPK3B nearly absent in Kuppe (0.8% cells express vs 23-33% in Linna-Kuosmanen)

**Note on Marker Selection:**
- **MSLN not used**: Not present in Kuppe dataset (28,380 genes vs 35,477 in Linna-Kuosmanen)
- **UPK3B low in Kuppe**: Kuppe cells selected by marker score contain many EMT-derived fibroblast-like cells that lost UPK3B expression; Linna-Kuosmanen cells are true mesothelial cells with explicit annotations

**Normalized Files:**
- `data/processed/epicardial_periheart_log.h5ad`
- `data/processed/epicardial_carebank_log.h5ad`
- `data/processed/epicardial_kuppe.h5ad` (already log-normalized)

---

## 4. Phase 2: Define Target Phenotypes

### 4.1 Proliferation Signature

**Core proliferation genes:**
- MKI67 (Ki-67)
- TOP2A (Topoisomerase II alpha)
- PCNA (Proliferating cell nuclear antigen)
- CDK1 (Cyclin-dependent kinase 1)
- CCNB1, CCNB2 (Cyclin B1/B2)
- CCNA2 (Cyclin A2)
- MCM2, MCM6 (Minichromosome maintenance)
- AURKA, AURKB (Aurora kinases)

### 4.2 EMT Signature

**EMT-promoting transcription factors (upregulated):**
- SNAI1, SNAI2 (Snail family)
- TWIST1, TWIST2
- ZEB1, ZEB2
- PRRX1

**Mesenchymal markers (upregulated):**
- VIM (Vimentin)
- CDH2 (N-cadherin)
- FN1 (Fibronectin)
- ACTA2 (alpha-SMA)

**Epithelial markers (downregulated):**
- CDH1 (E-cadherin)
- CLDN1, CLDN3 (Claudins)
- TJP1 (ZO-1)
- OCLN (Occludin)

**EMT score**: `emt_up_score - emt_down_score`

### 4.3 Cell State Classification Results

**Dual-method approach:**

**Method 1 - Condition-based:**
- Quiescent: `disease == 'normal'`
- Activated: `disease == 'myocardial infarction'` or `'myocardial ischemia'`
- Kuppe only: further split by `major_labl` (IZ/BZ=high, FZ=medium, RZ=low activation)
- Other diseases (valve disorder, heart failure): marked as 'other'

**Method 2 - Score-based:**
- Calculate `proliferation_score` (mean of 9 genes) and `emt_score` (up - down)
- Combine into z-normalized `activation_score = (prolif_z + emt_z) / 2`
- Activated: `activation_score > 75th percentile`

**Cross-validation results:**

| Metric | Value |
|--------|-------|
| Agreement rate | 66.8% |
| Activated (condition) EMT score | +0.171 |
| Quiescent (condition) EMT score | -0.113 |

**Consensus state distribution:**

| State | Count | Description |
|-------|-------|-------------|
| quiescent | 22,951 | Both methods agree: normal |
| activated | 6,162 | Both methods agree: activated |
| early_activated | 8,684 | Condition=MI but low scores |
| pre_activated | 5,780 | Condition=normal but high scores |
| ambiguous | 15,854 | Other disease conditions |

**Per-dataset activation scores:**

| Dataset | Mean Score | Quiescent | Activated |
|---------|------------|-----------|-----------|
| Kuppe_MI | +0.56 | 818 | 7,841 |
| PERIHEART | -0.05 | 10,391 | 6,439 |
| CAREBANK | -0.20 | 11,742 | 566 |

### 4.4 Prepare Full Communication Dataset

Extract sender cells from raw data and merge with epicardial cells for L-R analysis.

**Sender cell types:**
- Cardiomyocytes
- Macrophages / Myeloid cells
- Fibroblasts
- Endothelial cells

**Data source selection:**

| Dataset | Tissue | Sender cells available | Notes |
|---------|--------|------------------------|-------|
| **Kuppe** | Left ventricle | ✅ All 4 types | MI context, has spatial zones |
| **Linna-Kuosmanen** | Right atrium | ✅ All 4 types | Cardiac surgery patients |

> Both datasets can provide sender cells. Kuppe is preferred for MI-specific analysis due to spatial zone annotations (IZ/BZ/FZ/RZ). Linna-Kuosmanen can be used for validation or combined analysis.

**Output**: `full_communication_dataset.h5ad`
- Contains both sender and receiver (epicardial) cells
- Used by LIANA and NicheNet in Phase 3

### 4.5 Temporal Labeling (Not Used)

> Temporal labeling was not used because the available datasets lack precise timepoint annotations. Cell states were classified using condition-based and score-based methods instead (see 4.3).

---

## 5. Phase 3: Ligand-Receptor Analysis

### 5.1 Workflow Overview

```
Step 1: DEG Analysis (pyDESeq2)
  └── Input: epicardial_with_states.h5ad
  └── Compare: activated vs quiescent
  └── Output: differentially expressed genes

Step 2: LIANA Analysis
  └── Input: full_communication_dataset.h5ad (prepared in Phase 2)
  └── No DEG required
  └── Output: L-R pair rankings

Step 3: NicheNet Analysis
  └── Input: DEG geneset + full dataset
  └── Predict ligands causing DEG upregulation
  └── Output: Ligand activity rankings
```

### 5.2 DEG Analysis (pyDESeq2)

Find genes differentially expressed between activated and quiescent epicardial cells.

**Input**: `epicardial_with_states.h5ad` (Phase 2 output)

**Comparison**: `final_state == 'activated'` vs `final_state == 'quiescent'`

**Output**: Upregulated genes in activated cells → used as geneset for NicheNet

### 5.3 LIANA Analysis

LIANA integrates multiple methods (CellPhoneDB, NATMI, Connectome, etc.) for consensus L-R rankings.

**Input**: `full_communication_dataset.h5ad` (Phase 2 output, contains sender + receiver cells)

**Key function**: `li.mt.rank_aggregate(adata, groupby='cell_type', resource_name='consensus', expr_prop=0.1)`

**Output**: Filter `adata.uns['liana_res']` for `target == 'mesothelial cell'`, rank by `magnitude_rank`

### 5.4 NicheNet Analysis

NicheNet predicts which ligands from sender cells best explain transcriptional changes in receiver cells.

**Input**:
- DEG geneset from 5.2
- `full_communication_dataset.h5ad` (Phase 2 output)

**Required data** (from Zenodo):
- `ligand_target_matrix.rds`
- `lr_network.rds`

**Key function**: `predict_ligand_activities(geneset, ligand_target_matrix, potential_ligands)`

### 5.5 Expected Results

**Priority ligands**: FGF10, PDGFA, TGFB1, WNT5A, HGF

**FGF10 receptors**: FGFR1, FGFR2 (specifically FGFR2b isoform), FGFR3

**Validation**: FGF10-FGFR2 should rank highly in both LIANA and NicheNet results

### 5.6 Backup Resources

Additional L-R databases (not used in main analysis):
- [CellTalkDB](http://tcm.zju.edu.cn/celltalkdb/)
- NATMI database
- Ramilowski et al. ligand-receptor pairs

**Filter criteria** (if using raw databases):
- Receptors expressed in epicardial cells (mean > 0.1, pct > 10%)
- Ligands expressed in sender populations
- Documented signaling activity

---

## 6. Phase 4: Machine Learning Model (Geneformer)

### 6.1 Feature Engineering

For each candidate ligand-receptor pair, construct features:

**A) Expression features:**
- Receptor expression level in epicardium (mean, max, pct.expressed)
- Ligand expression in sender cells
- Fold change post-MI vs sham

**B) NicheNet/LIANA scores:**
- Ligand activity score (Pearson correlation)
- Communication probability
- Ligand-target regulatory potential

**C) Correlation features:**
- Correlation of receptor expression with proliferation score
- Correlation with EMT score
- Temporal correlation with activation

**D) Pathway features:**
- Pathway membership (FGFR, PDGFR, TGFBR, Wnt, Notch)
- Downstream target enrichment

**E) Druggability features:**
- Existing approved agonists/antagonists
- Protein class (GPCR, RTK, cytokine receptor)

### 6.2 Label Construction

**Positive ligand-receptor pairs (known epicardial activators):**
- FGF10 → FGFR2b
- FGF2 → FGFR1/2
- PDGFA/B → PDGFRA/B
- TGFB1/2/3 → TGFBR1/2
- WNT ligands → FZD receptors
- HGF → MET
- IGF1 → IGF1R

**Negative pairs:**
- Ligand-receptor pairs not expressed in relevant cell types
- Pairs with no documented effect on epicardium

### 6.3 Geneformer Fine-tuning

Geneformer is a foundation model pretrained on ~95M single-cell transcriptomes.
Fine-tuning enables cell state classification and **in silico perturbation**.

**Requirements:**
- GPU with at least 16GB VRAM (A100, V100, or RTX 3090+)
- Python 3.8+
- PyTorch, Transformers, Datasets libraries

**Installation**: `pip install transformers datasets torch geneformer`

**Workflow:**

1. **Prepare Data**: Use `TranscriptomeTokenizer` to tokenize transcriptomes. Label cells as 0 (quiescent) or 1 (activated). Gene names must be Ensembl IDs.

2. **Fine-tune**: Load `ctheodoris/Geneformer`, fine-tune with `BertForSequenceClassification` (num_labels=2). Recommended: 5 epochs, batch_size=8, warmup_steps=500.

3. **In Silico Perturbation**: Use `InSilicoPerturber` with `perturb_type="delete"` to test candidate receptors (FGFR2, FGFR1, PDGFRA, PDGFRB, TGFBR1/2, MET, IGF1R, FZD1/2/7, BMPR1A/2).

4. **Rank Results**: Calculate `activation_effect = prob_class1_perturbed - prob_class1_original`. Most negative = most important for activation.

5. **Validate**: FGFR2 deletion should show strong negative effect (prevents activation).

### 6.4 Model Evaluation

**Evaluation metrics:**
- Accuracy on held-out test set
- F1 score (balanced metric for imbalanced classes)
- AUROC for cell state prediction

**Positive control validation:**
- FGFR2 deletion should significantly reduce activation probability
- FGF pathway receptors should cluster among top hits

---

## 7. Phase 5: Model Interpretation and Prioritization

### 7.1 Feature Importance

Extract and visualize:
- Attention weights for transformer models
- Perturbation effect sizes
- SHAP values (if using ensemble methods)

**Key questions:**
- Is the model learning biologically meaningful patterns?
- Which features drive predictions for top candidates?

### 7.2 Ranked Output

Generate final ranked list:
1. Rank all receptor candidates by perturbation effect
2. Annotate with: known biology, druggability status, expression specificity, pathway membership
3. Flag novel predictions (not in training positives)

### 7.3 Pathway Enrichment of Top Candidates

Check if top candidates converge on expected pathways:
- FGF signaling
- PDGF signaling
- TGF-beta superfamily
- Wnt signaling
- Receptor tyrosine kinases
- Hippo/YAP pathway

*If top hits are scattered randomly, revisit model design.*

---

## 8. Tools and Resources

### Software Tools

| Category | Tool | URL |
|----------|------|-----|
| **scRNA-seq** | Scanpy | https://scanpy.readthedocs.io/ |
| **scRNA-seq** | scVI | https://scvi-tools.org/ |
| **Cell communication** | LIANA | https://github.com/saezlab/liana |
| **Cell communication** | NicheNet | https://github.com/saeyslab/nichenetr |
| **Cell communication** | CellChat | https://github.com/sqjin/CellChat |
| **Cell communication** | CellPhoneDB | https://www.cellphonedb.org/ |
| **Deep learning** | Geneformer | https://huggingface.co/ctheodoris/Geneformer |
| **Deep learning** | TxGNN | https://github.com/mims-harvard/TxGNN |
| **Protein design** | RFdiffusion | https://github.com/RosettaCommons/RFdiffusion |
| **Protein design** | BindCraft | https://github.com/martinpacesa/BindCraft |

### Databases

| Category | Resource | URL |
|----------|----------|-----|
| **Single-cell data** | CellxGene | https://cellxgene.cziscience.com/ |
| **Cardiac atlas** | Heart Cell Atlas | https://www.heartcellatlas.org/ |
| **Ligand-receptor** | CellTalkDB | http://tcm.zju.edu.cn/celltalkdb/ |
| **Ligand-receptor** | CellPhoneDB | https://www.cellphonedb.org/ |
| **Drug-target** | DrugBank | https://go.drugbank.com/ |
| **Drug-target** | ChEMBL | https://www.ebi.ac.uk/chembl/ |

---

## 9. Key References

### Spatial/Single-Cell Cardiac Atlases

1. Kuppe C, et al. (2022) "Spatial multi-omic map of human myocardial infarction" *Nature* 608:766-777. [DOI](https://doi.org/10.1038/s41586-022-05060-x)

2. Linna-Kuosmanen S, et al. (2024) "Single-nucleus RNA sequencing of human heart" *Cell Reports Medicine*

3. Molenaar B, et al. (2021) "Single-cell transcriptomics following ischemic injury identifies a role for B2M in cardiac repair" *Commun Biol* 4:146

4. Yamada S, et al. (2022) "Spatiotemporal transcriptome analysis reveals critical roles for mechano-sensing genes at the border zone" *Nat Cardiovasc Res* 1:1039-1055

### AI/ML for Target Discovery

5. Theodoris CV, et al. (2023) "Transfer learning enables predictions in network biology" *Nature* 618:616-624 — **Geneformer**

6. Huang K, et al. (2024) "TxGNN: A foundation model for clinician-centered drug repurposing" *Nat Med*

7. Singh R, et al. (2023) "Contrastive learning in protein language space predicts interactions between drugs and protein targets" *PNAS* 120:e2220778120 — **ConPLex**

### Protein Design

8. Watson JL, et al. (2023) "De novo design of protein structure and function with RFdiffusion" *Nature* 620:1089-1100

9. Pacesa M, et al. (2025) "BindCraft: One-shot design of functional protein binders" *Nature*

### Cell-Cell Communication

10. Browaeys R, et al. (2020) "NicheNet: modeling intercellular communication by linking ligands to target genes" *Nat Methods* 17:159-162

11. Jin S, et al. (2021) "Inference and analysis of cell-cell communication using CellChat" *Nat Commun* 12:1088

12. Dimitrov D, et al. (2022) "Comparison of methods and resources for cell-cell communication inference from single-cell RNA-Seq data" *Nat Commun* 13:3224 — **LIANA**

13. Efremova M, et al. (2020) "CellPhoneDB: inferring cell–cell communication from combined expression of multi-subunit ligand–receptor complexes" *Nat Protoc* 15:1484-1506

### Epicardial Biology

14. Cao J, Poss KD (2018) "The epicardium as a hub for heart regeneration" *Nat Rev Cardiol* 15:631-647

15. Wang S, et al. (2015) "Epicardial regeneration is guided by cardiac outflow tract and Hedgehog signalling" *Nature* 522:226-230

---

## Summary Notes

- **FGF10/FGFR2b** serves as positive control throughout pipeline
- Cardiomyocytes are **SENDER** cells; epicardial cells are **RECEIVER** cells
- Epicardial markers: WT1, TBX18, ALDH1A2, UPK3B, MSLN (**TCF21 excluded** - low expression)
- Mesothelial cells in Linna-Kuosmanen = epicardial cells (~49,000 available!)
- Data format: `.h5ad` (AnnData/HDF5), use `backed='r'` for large files
- GPU with **16GB+ VRAM** required for Geneformer fine-tuning
- Start with **LIANA/NicheNet** analysis before Geneformer for initial candidate identification
