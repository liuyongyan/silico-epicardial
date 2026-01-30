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
│   │   ├── kuppe_cell_metadata.csv
│   │   └── kuppe_gene_metadata.csv
│   └── references/
├── scripts/
│   └── 01_preprocessing/
│       └── analyze_kuppe.py
├── results/
│   └── tables/
│       ├── kuppe_celltype_summary.csv
│       └── kuppe_sample_metadata.csv
├── docs/
│   └── H5AD_FORMAT_GUIDE.md
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

### Next Steps

- [ ] Extract and merge epicardial/mesothelial cells from all datasets
- [ ] Define proliferation and EMT signatures
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

```python
import scanpy as sc
import anndata as ad

# Load data (use backed mode for large files)
adata = ad.read_h5ad("data.h5ad", backed='r')

# For in-memory processing
adata = ad.read_h5ad("data.h5ad")

# QC metrics
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)

# Filter
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_cells(adata, max_genes=6000)
adata = adata[adata.obs.pct_counts_mt < 15, :]

# Normalize
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Variable genes and PCA
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
sc.tl.pca(adata)

# Batch correction with Harmony
import scanpy.external as sce
sce.pp.harmony_integrate(adata, key='batch')

# Clustering
sc.pp.neighbors(adata, use_rep='X_pca_harmony')
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.5)
```

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

```python
# Score cells for epicardial identity
epicardial_markers = ['WT1', 'TBX18', 'ALDH1A2', 'UPK3B', 'MSLN']

# If using gene symbols
sc.tl.score_genes(adata, epicardial_markers, score_name='epicardial_score')

# Subset epicardial population (by annotation or score)
# For Linna-Kuosmanen data:
epicardial = adata[adata.obs['cell_type'] == 'mesothelial cell']

# Or by score threshold:
threshold = adata.obs['epicardial_score'].quantile(0.95)
epicardial = adata[adata.obs['epicardial_score'] > threshold]

# Subcluster to identify states
sc.pp.neighbors(epicardial, n_neighbors=15)
sc.tl.leiden(epicardial, resolution=0.3)
```

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

```python
proliferation_genes = ['MKI67', 'TOP2A', 'PCNA', 'CDK1', 'CCNB1',
                       'CCNB2', 'CCNA2', 'MCM2', 'MCM6', 'AURKA']

sc.tl.score_genes(adata, proliferation_genes, score_name='proliferation_score')

# Or use cell cycle scoring
sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)
```

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

```python
emt_up_genes = ['SNAI1', 'SNAI2', 'TWIST1', 'ZEB1', 'ZEB2',
                'VIM', 'CDH2', 'FN1', 'ACTA2']
emt_down_genes = ['CDH1', 'CLDN1', 'TJP1', 'OCLN']

sc.tl.score_genes(adata, emt_up_genes, score_name='emt_up_score')
sc.tl.score_genes(adata, emt_down_genes, score_name='emt_down_score')

# Combined EMT score
adata.obs['emt_score'] = adata.obs['emt_up_score'] - adata.obs['emt_down_score']
```

### 4.3 Temporal Labeling

Label cells by:
- Timepoint post-MI (Day 0, 1, 3, 7, 14, etc.)
- Spatial zone (ischemic, border, remote)
- Activation state (quiescent vs activated)

**Target transition**: Quiescent (sham/day 0) → Activated (day 3-7)

---

## 5. Phase 3: Ligand-Receptor Analysis

### 5.1 Interaction Database Resources

Curated resources to merge:
- [CellPhoneDB](https://www.cellphonedb.org/)
- [CellTalkDB](http://tcm.zju.edu.cn/celltalkdb/)
- [NicheNet](https://github.com/saeyslab/nichenetr) ligand-target matrix
- NATMI database
- Ramilowski et al. ligand-receptor pairs

**Filter criteria:**
- Receptors expressed in epicardial cells (mean > 0.1, pct > 10%)
- Ligands expressed in sender populations
- Documented signaling activity

### 5.2 LIANA Analysis (Python, Recommended)

LIANA combines multiple methods (CellPhoneDB, NATMI, Connectome, etc.) and provides consensus rankings.

```python
import liana as li

# Prepare data
adata.raw = adata  # store raw counts
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)

# Run LIANA with multiple methods
li.mt.rank_aggregate(
    adata,
    groupby='cell_type',
    resource_name='consensus',  # or 'cellphonedb', 'celltalkdb'
    expr_prop=0.1,
    verbose=True
)

# Results stored in adata.uns['liana_res']
liana_results = adata.uns['liana_res']

# Filter for epicardial/mesothelial as receiver
epi_interactions = liana_results[liana_results['target'] == 'mesothelial cell']

# Top ligand-receptor pairs
top_pairs = epi_interactions.nsmallest(50, 'magnitude_rank')

# Check if FGF10-FGFR2 is in top pairs (positive control)
fgf10_pairs = epi_interactions[epi_interactions['ligand_complex'].str.contains('FGF10')]
print(fgf10_pairs)
```

### 5.3 NicheNet Analysis (R Alternative)

NicheNet predicts which ligands best explain transcriptional changes in receiver cells.

```r
library(nichenetr)

# Download required matrices
ligand_target_matrix <- readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
lr_network <- readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))

# Define sender and receiver
sender_celltypes <- c("Cardiomyocyte", "Macrophage", "Fibroblast", "Endothelial")
receiver_celltype <- "Epicardial"

# Get expressed genes
expressed_genes_receiver <- get_expressed_genes(receiver_celltype, seurat_obj, pct = 0.10)

# Define gene set of interest (activation signature)
geneset_oi <- c("MKI67", "TOP2A", "VIM", "SNAI1", "SNAI2", "TWIST1",
                "POSTN", "COL1A1", "FN1", "ACTA2")

# Predict ligand activities
ligand_activities <- predict_ligand_activities(
  geneset = geneset_oi,
  background_expressed_genes = expressed_genes_receiver,
  ligand_target_matrix = ligand_target_matrix,
  potential_ligands = potential_ligands
)

# Rank ligands - FGF10 should appear in top if working correctly
top_ligands <- ligand_activities %>%
  arrange(-pearson) %>%
  head(20)
```

### 5.4 Map Ligands to Receptors

```python
# From LIANA results, extract receptor partners for top ligands
top_ligands = ['FGF10', 'PDGFA', 'TGFB1', 'WNT5A', 'HGF']

receptor_mapping = liana_results[
    liana_results['ligand_complex'].isin(top_ligands)
][['ligand_complex', 'receptor_complex', 'magnitude_rank']].drop_duplicates()

# Priority receptors for FGF10:
# FGFR1, FGFR2 (specifically FGFR2b isoform), FGFR3
```

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

```bash
pip install transformers datasets torch
pip install geneformer  # or clone from GitHub
```

**Step 1: Prepare Data**

```python
import scanpy as sc
from geneformer import TranscriptomeTokenizer

# Load epicardial cells from processed dataset
adata = sc.read_h5ad("epicardial_cells.h5ad")

# Add cell state labels
# 0 = quiescent (sham/control), 1 = activated (post-MI day 3-7)
adata.obs['activation_state'] = adata.obs['condition'].map({
    'control': 0, 'sham': 0,
    'MI_day3': 1, 'MI_day7': 1,
    'ischemic': 1, 'border': 1
})

# Ensure gene names are Ensembl IDs (required by Geneformer)
# CellxGene data already uses Ensembl IDs as var_names

# Tokenize transcriptomes
tokenizer = TranscriptomeTokenizer(
    custom_attr_name_dict={"cell_type": "cell_type"},
    nproc=4
)
tokenized_data = tokenizer.tokenize_data(
    adata,
    output_directory="./tokenized_data/",
    output_prefix="epicardial"
)
```

**Step 2: Fine-tune for Cell State Classification**

```python
from transformers import BertForSequenceClassification, Trainer, TrainingArguments
from datasets import load_from_disk
from sklearn.metrics import accuracy_score, f1_score
import numpy as np

# Load tokenized data
dataset = load_from_disk("./tokenized_data/epicardial.dataset")
dataset = dataset.train_test_split(test_size=0.2, seed=42)

# Load pretrained Geneformer
model = BertForSequenceClassification.from_pretrained(
    "ctheodoris/Geneformer",
    num_labels=2,  # quiescent vs activated
    output_attentions=False,
    output_hidden_states=False
)

# Training arguments
training_args = TrainingArguments(
    output_dir="./geneformer_epicardial/",
    num_train_epochs=5,
    per_device_train_batch_size=8,
    per_device_eval_batch_size=8,
    warmup_steps=500,
    weight_decay=0.01,
    logging_dir="./logs/",
    logging_steps=100,
    evaluation_strategy="epoch",
    save_strategy="epoch",
    load_best_model_at_end=True,
    metric_for_best_model="accuracy"
)

# Define metrics
def compute_metrics(pred):
    labels = pred.label_ids
    preds = np.argmax(pred.predictions, axis=1)
    return {
        "accuracy": accuracy_score(labels, preds),
        "f1": f1_score(labels, preds)
    }

# Train
trainer = Trainer(
    model=model,
    args=training_args,
    train_dataset=dataset["train"],
    eval_dataset=dataset["test"],
    compute_metrics=compute_metrics
)
trainer.train()
```

**Step 3: In Silico Perturbation**

```python
from geneformer import InSilicoPerturber

# Initialize perturber with fine-tuned model
perturber = InSilicoPerturber(
    perturb_type="delete",  # or "overexpress"
    model_type="CellClassifier",
    num_classes=2,
    emb_mode="cell",
    filter_data={"cell_type": ["epicardial"]}
)

# Define genes to perturb (candidate receptors)
candidate_receptors = [
    "FGFR2", "FGFR1", "PDGFRA", "PDGFRB",
    "TGFBR1", "TGFBR2", "MET", "IGF1R",
    "FZD1", "FZD2", "FZD7", "BMPR1A", "BMPR2"
]

# Run perturbation analysis
perturber.perturb_data(
    model_directory="./geneformer_epicardial/",
    input_data_file="./tokenized_data/epicardial.dataset",
    output_directory="./perturbation_results/",
    genes_to_perturb=candidate_receptors
)

# Receptors whose deletion PREVENTS activation (shifts cells toward quiescent)
# are potential therapeutic targets for ACTIVATING epicardium
```

**Step 4: Rank Receptors by Perturbation Effect**

```python
import pandas as pd
import matplotlib.pyplot as plt

# Load perturbation results
results = pd.read_csv("./perturbation_results/perturbation_stats.csv")

# Calculate effect size: shift toward quiescent state upon deletion
# Negative values = receptor is required for activation
results['activation_effect'] = results['prob_class1_perturbed'] - results['prob_class1_original']

# Rank by effect (most negative = most important for activation)
results_sorted = results.sort_values('activation_effect')

# Visualize
plt.figure(figsize=(10, 8))
plt.barh(results_sorted['gene'], results_sorted['activation_effect'])
plt.xlabel('Change in Activation Probability (Deletion Effect)')
plt.ylabel('Receptor Gene')
plt.title('In Silico Perturbation: Receptor Importance for Epicardial Activation')
plt.axvline(x=0, color='k', linestyle='--')
plt.tight_layout()
plt.savefig('receptor_perturbation_effects.png', dpi=300)

# Top candidates: receptors with most negative activation_effect
top_receptors = results_sorted.head(10)
print("Top receptor candidates (deletion inhibits activation):")
print(top_receptors[['gene', 'activation_effect']])
```

**Step 5: Validate Positive Control**

```python
# FGFR2 should show strong negative effect (deletion prevents activation)
fgfr2_result = results[results['gene'] == 'FGFR2']
print(f"FGFR2 perturbation effect: {fgfr2_result['activation_effect'].values[0]}")

# If FGFR2 is not in top 10, investigate:
# - Check FGFR2 expression in epicardial cells
# - Verify correct gene ID mapping
# - Consider FGFR1 as alternative (may compensate)
```

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

### Epicardial Biology

12. Cao J, Poss KD (2018) "The epicardium as a hub for heart regeneration" *Nat Rev Cardiol* 15:631-647

13. Wang S, et al. (2015) "Epicardial regeneration is guided by cardiac outflow tract and Hedgehog signalling" *Nature* 522:226-230

---

## Summary Notes

- **FGF10/FGFR2b** serves as positive control throughout pipeline
- Cardiomyocytes are **SENDER** cells; epicardial cells are **RECEIVER** cells
- Epicardial markers: WT1, TBX18, ALDH1A2, UPK3B, MSLN (**TCF21 excluded** - low expression)
- Mesothelial cells in Linna-Kuosmanen = epicardial cells (~49,000 available!)
- Data format: `.h5ad` (AnnData/HDF5), use `backed='r'` for large files
- GPU with **16GB+ VRAM** required for Geneformer fine-tuning
- Start with **LIANA/NicheNet** analysis before Geneformer for initial candidate identification
