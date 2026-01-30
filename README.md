# AI-Driven Epicardial Target Discovery in Myocardial Infarction

> **Purpose**: Identify potential therapeutic targets on epicardial cells that promote proliferation and EMT for cardioprotection after myocardial infarction

## Table of Contents

- [Overview and Rationale](#overview-and-rationale)
- [Phase 1: Data Acquisition and Preprocessing](#phase-1-data-acquisition-and-preprocessing)
- [Phase 2: Define Target Phenotypes](#phase-2-define-target-phenotypes)
- [Phase 3: Ligand-Receptor Analysis](#phase-3-ligand-receptor-analysis)
- [Phase 4: Machine Learning Model](#phase-4-machine-learning-model)
- [Phase 5: Model Interpretation and Prioritization](#phase-5-model-interpretation-and-prioritization)
- [Phase 6: Validation Strategy](#phase-6-validation-strategy)
- [Tools and Resources](#tools-and-resources)
- [Key References](#key-references)

---

## Overview and Rationale

### Goal

Identify receptors on epicardial cells whose activation leads to:
- Epicardial cell proliferation
- Epithelial-to-mesenchymal transition (EMT)
- Cardioprotective paracrine secretion

### Biological Context

After MI, epicardial cells can be activated to support cardiac repair. FGF10 delivered via pericardial cavity activates epicardium through FGFR2b. Epicardium-derived cells (EPDCs) contribute to repair through:
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

## Phase 1: Data Acquisition and Preprocessing

### 1.1 Recommended Datasets

#### Human MI Data (highest priority)
- **Kuppe et al. 2022** (Nature) - Spatial multi-omic map of human MI
- Contains: snRNA-seq, snATAC-seq, spatial transcriptomics
- 23 patients, 250,000+ cells, multiple timepoints and zones
- [Paper](https://www.nature.com/articles/s41586-022-05060-x)

#### Mouse MI Data
- **Farbehi et al. 2019** (eLife) - Non-cardiomyocyte populations
- GEO: GSE130699
- Contains: >30,000 cells including epicardial populations
- Timepoints: Day 3 and 7 post-MI
- [Paper](https://elifesciences.org/articles/43882)

#### Additional datasets
- GSE136088 - Mouse MI scRNA-seq
- GSE185289 - Pig cardiac regeneration model

### 1.2 Standard scRNA-seq Processing

**Quality Control:**
- Filter cells: 200 < nFeature_RNA < 6000
- Mitochondrial content: < 10-20% (tissue-dependent)
- Doublet removal: Scrublet or DoubletFinder

**Normalization options:**
- SCTransform (Seurat) - recommended for integration
- scran pooling method
- Standard log-normalization

**Batch correction/Integration:**
- Harmony (fast, effective)
- scVI (deep learning-based)
- Seurat RPCA integration

<details>
<summary><b>Code Example (Seurat)</b></summary>

```r
library(Seurat)
library(harmony)

# Load and QC
seurat_obj <- CreateSeuratObject(counts = raw_counts)
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
seurat_obj <- subset(seurat_obj,
                     subset = nFeature_RNA > 200 &
                              nFeature_RNA < 6000 &
                              percent.mt < 15)

# Normalize and integrate
seurat_obj <- SCTransform(seurat_obj)
seurat_obj <- RunPCA(seurat_obj)
seurat_obj <- RunHarmony(seurat_obj, group.by.vars = "batch")
seurat_obj <- RunUMAP(seurat_obj, reduction = "harmony", dims = 1:30)
seurat_obj <- FindNeighbors(seurat_obj, reduction = "harmony", dims = 1:30)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
```

</details>

### 1.3 Epicardial Cell Identification

**Canonical epicardial markers:**
- WT1 (Wilms tumor 1) - master regulator
- TBX18 (T-box 18)
- TCF21 (Transcription factor 21)
- ALDH1A2 (Retinaldehyde dehydrogenase)
- UPK3B (Uroplakin 3B)
- MSLN (Mesothelin)

**Activated/EMT epicardial markers:**
- POSTN (Periostin)
- COL1A1, COL3A1 (Collagens)
- ACTA2 (alpha-SMA)
- VIM (Vimentin)

<details>
<summary><b>Code Example</b></summary>

```r
# Identify epicardial clusters
epicardial_markers <- c("WT1", "TBX18", "TCF21", "ALDH1A2", "UPK3B")

# Score cells
seurat_obj <- AddModuleScore(seurat_obj,
                             features = list(epicardial_markers),
                             name = "Epicardial_Score")

# Subset epicardial population
epicardial_cells <- subset(seurat_obj,
                           subset = Epicardial_Score1 > threshold)

# Subcluster to identify states
epicardial_cells <- FindNeighbors(epicardial_cells, dims = 1:20)
epicardial_cells <- FindClusters(epicardial_cells, resolution = 0.3)
```

</details>

---

## Phase 2: Define Target Phenotypes

### 2.1 Proliferation Signature

**Core proliferation genes:**
- MKI67 (Ki-67)
- TOP2A (Topoisomerase II alpha)
- PCNA (Proliferating cell nuclear antigen)
- CDK1 (Cyclin-dependent kinase 1)
- CCNB1, CCNB2 (Cyclin B1/B2)
- CCNA2 (Cyclin A2)
- MCM2, MCM6 (Minichromosome maintenance)
- AURKA, AURKB (Aurora kinases)

<details>
<summary><b>Code Example</b></summary>

```r
proliferation_genes <- c("MKI67", "TOP2A", "PCNA", "CDK1", "CCNB1",
                         "CCNB2", "CCNA2", "MCM2", "MCM6", "AURKA")

seurat_obj <- AddModuleScore(seurat_obj,
                             features = list(proliferation_genes),
                             name = "Proliferation")

# Or use CellCycleScoring
seurat_obj <- CellCycleScoring(seurat_obj,
                                s.features = cc.genes$s.genes,
                                g2m.features = cc.genes$g2m.genes)
```

</details>

### 2.2 EMT Signature

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

<details>
<summary><b>Code Example</b></summary>

```r
emt_up_genes <- c("SNAI1", "SNAI2", "TWIST1", "ZEB1", "ZEB2",
                  "VIM", "CDH2", "FN1")
emt_down_genes <- c("CDH1", "CLDN1", "TJP1", "OCLN")

seurat_obj <- AddModuleScore(seurat_obj,
                             features = list(emt_up_genes),
                             name = "EMT_up")
seurat_obj <- AddModuleScore(seurat_obj,
                             features = list(emt_down_genes),
                             name = "EMT_down")

# Combined EMT score
seurat_obj$EMT_score <- seurat_obj$EMT_up1 - seurat_obj$EMT_down1
```

</details>

### 2.3 Temporal Labeling

Label cells by:
- Timepoint post-MI (Day 0, 1, 3, 7, 14, etc.)
- Spatial zone (ischemic, border, remote)
- Activation state (quiescent vs activated)

**Target transition**: Quiescent (sham/day 0) → Activated (day 3-7)

---

## Phase 3: Ligand-Receptor Analysis

### 3.1 Interaction Database Resources

Curated resources to merge:
- [CellPhoneDB](https://www.cellphonedb.org/)
- [CellTalkDB](http://tcm.zju.edu.cn/celltalkdb/)
- NicheNet ligand-target matrix
- NATMI database
- Ramilowski et al. ligand-receptor pairs

**Filter criteria:**
- Receptors expressed in epicardial cells (mean > 0.1, pct > 10%)
- Ligands expressed in sender populations
- Documented signaling activity

### 3.2 NicheNet Analysis (Recommended Primary Method)

NicheNet predicts which ligands best explain transcriptional changes in receiver cells based on prior knowledge of signaling networks.

<details>
<summary><b>Installation</b></summary>

```r
# Install NicheNet
devtools::install_github("saeyslab/nichenetr")
library(nichenetr)

# Download required matrices
ligand_target_matrix <- readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
lr_network <- readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
weighted_networks <- readRDS(url("https://zenodo.org/record/3260758/files/weighted_networks.rds"))
```

</details>

<details>
<summary><b>Analysis Workflow</b></summary>

```r
# Define sender and receiver
sender_celltypes <- c("Cardiomyocyte", "Macrophage", "Fibroblast", "Endothelial")
receiver_celltype <- "Epicardial"

# Get expressed genes
expressed_genes_receiver <- get_expressed_genes(receiver_celltype,
                                                 seurat_obj,
                                                 pct = 0.10)
expressed_genes_sender <- unique(unlist(lapply(sender_celltypes, function(ct) {
  get_expressed_genes(ct, seurat_obj, pct = 0.10)
})))

# Define gene set of interest (your activation signature)
geneset_oi <- c("MKI67", "TOP2A", "VIM", "SNAI1", "SNAI2", "TWIST1",
                "POSTN", "COL1A1", "FN1", "ACTA2")

# Background genes
background_genes <- expressed_genes_receiver

# Get potential ligands
potential_ligands <- lr_network %>%
  filter(from %in% expressed_genes_sender & to %in% expressed_genes_receiver) %>%
  pull(from) %>% unique()

# Predict ligand activities
ligand_activities <- predict_ligand_activities(
  geneset = geneset_oi,
  background_expressed_genes = background_genes,
  ligand_target_matrix = ligand_target_matrix,
  potential_ligands = potential_ligands
)

# Rank ligands
ligand_activities <- ligand_activities %>%
  arrange(-pearson) %>%
  mutate(rank = row_number())

# Top ligands
top_ligands <- ligand_activities %>%
  top_n(20, pearson) %>%
  pull(test_ligand)

# FGF10 should appear in top ligands if model is working correctly
print(ligand_activities %>% filter(test_ligand == "FGF10"))
```

</details>

### 3.3 Map Ligands to Receptors

```r
# Get receptors for top ligands
top_receptors <- lr_network %>%
  filter(from %in% top_ligands & to %in% expressed_genes_receiver) %>%
  select(from, to) %>%
  rename(ligand = from, receptor = to)

# Priority receptors for FGF10:
# FGFR1, FGFR2 (specifically FGFR2b isoform), FGFR3
```

### 3.4 CellChat Analysis (Complementary Method)

<details>
<summary><b>Code Example</b></summary>

```r
library(CellChat)

# Create CellChat object
cellchat <- createCellChat(object = seurat_obj, group.by = "cell_type")

# Set database
CellChatDB <- CellChatDB.human  # or CellChatDB.mouse
cellchat@DB <- CellChatDB

# Preprocessing
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# Inference
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

# Extract interactions targeting epicardium
df.net <- subsetCommunication(cellchat, targets.use = "Epicardial")

# Visualize
netVisual_chord_gene(cellchat,
                     sources.use = sender_celltypes,
                     targets.use = "Epicardial",
                     lab.cex = 0.5)
```

</details>

---

## Phase 4: Machine Learning Model

### 4.1 Feature Engineering

For each candidate ligand-receptor pair, construct features:

| Category | Features |
|----------|----------|
| **Expression** | Receptor expression level, ligand expression, fold change post-MI vs sham |
| **NicheNet/CellChat** | Ligand activity score, communication probability, regulatory potential |
| **Correlation** | Correlation with proliferation score, EMT score, temporal activation |
| **Pathway** | Pathway membership, downstream target enrichment, cardiac regeneration involvement |
| **Druggability** | Existing agonists/antagonists, protein class, structural information |

### 4.2 Label Construction

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
- Receptors for pro-apoptotic or anti-proliferative signals

### 4.3 Model Options

<details>
<summary><b>Option A: Gradient Boosting (XGBoost)</b></summary>

```python
import xgboost as xgb
from sklearn.model_selection import LeaveOneOut
import shap

# Prepare features and labels
X = feature_matrix  # ligand-receptor pairs x features
y = labels  # 1 for positive, 0 for negative

# Leave-one-out cross-validation (for small positive set)
loo = LeaveOneOut()
predictions = []

for train_idx, test_idx in loo.split(X):
    X_train, X_test = X[train_idx], X[test_idx]
    y_train, y_test = y[train_idx], y[test_idx]

    model = xgb.XGBClassifier(
        n_estimators=100,
        max_depth=4,
        learning_rate=0.1,
        scale_pos_weight=len(y[y==0])/len(y[y==1])  # handle imbalance
    )
    model.fit(X_train, y_train)
    predictions.append(model.predict_proba(X_test)[0, 1])

# SHAP for interpretability
explainer = shap.TreeExplainer(model)
shap_values = explainer.shap_values(X)
shap.summary_plot(shap_values, X, feature_names=feature_names)
```

</details>

<details>
<summary><b>Option B: Graph Neural Network</b></summary>

```python
import torch
from torch_geometric.nn import GCNConv, GATConv

class LigandReceptorGNN(torch.nn.Module):
    def __init__(self, num_node_features, hidden_channels):
        super().__init__()
        self.conv1 = GATConv(num_node_features, hidden_channels)
        self.conv2 = GATConv(hidden_channels, hidden_channels)
        self.lin = torch.nn.Linear(hidden_channels * 2, 1)

    def forward(self, x, edge_index, ligand_idx, receptor_idx):
        x = self.conv1(x, edge_index).relu()
        x = self.conv2(x, edge_index)

        # Concatenate ligand and receptor embeddings
        ligand_emb = x[ligand_idx]
        receptor_emb = x[receptor_idx]
        pair_emb = torch.cat([ligand_emb, receptor_emb], dim=-1)

        return torch.sigmoid(self.lin(pair_emb))
```

</details>

<details>
<summary><b>Option C: Geneformer Fine-tuning</b></summary>

```python
# Requires GPU and huggingface transformers
from transformers import BertForSequenceClassification
from geneformer import TranscriptomeTokenizer

# Tokenize epicardial cells
tokenizer = TranscriptomeTokenizer()
tokenized_data = tokenizer.tokenize_data(epicardial_adata)

# Fine-tune for cell state classification
model = BertForSequenceClassification.from_pretrained(
    "ctheodoris/Geneformer",
    num_labels=2  # quiescent vs activated
)

# In silico perturbation
# Delete/activate receptor genes and measure embedding shift
```

</details>

### 4.4 Training and Validation

**Cross-validation strategy:**
- Leave-one-out CV given small positive set
- Hold out FGF10/FGFR2b as primary validation case
- If model ranks FGF10 highly when blinded, proceed with confidence

**Metrics:**
- AUROC, AUPRC (for ranking)
- Precision at top-k
- Recall of known positives

---

## Phase 5: Model Interpretation and Prioritization

### 5.1 Feature Importance

Extract and visualize:
- SHAP values for gradient boosting models
- Attention weights for transformer models
- Node importance for GNN models

**Key questions:**
- Is the model learning biologically meaningful patterns?
- Which features drive predictions for top candidates?

### 5.2 Ranked Output

Generate final ranked list:
1. Rank all receptor candidates by model score
2. Annotate with known biology, druggability status, expression specificity, pathway membership
3. Flag novel predictions (not in training positives)

### 5.3 Pathway Enrichment

Check if top candidates converge on expected pathways:
- FGF signaling
- PDGF signaling
- TGF-beta superfamily
- Wnt signaling
- Receptor tyrosine kinases
- Hippo/YAP pathway

*If top hits are scattered randomly, revisit model design.*

---

## Phase 6: Validation Strategy

### 6.1 In Silico Validation

- **Independent dataset validation**: Test on datasets not used in training
- **Geneformer in silico perturbation**: Virtual knockout/overexpression
- **Literature validation**: Cross-reference with published studies, GWAS

### 6.2 In Vitro Validation

**Cell models:**
- Primary epicardial cells from mouse hearts
- Human iPSC-derived epicardial cells

**Functional assays:**
- Proliferation: EdU incorporation, Ki67 staining
- EMT markers: Immunofluorescence for VIM, SNAI1, CDH1
- Migration/invasion: Transwell assays
- Ligand stimulation with RNA-seq or qPCR readout

### 6.3 In Vivo Validation

- **Pericardial delivery model**: Deliver candidate ligands via pericardial injection in mouse MI model
- **Functional outcomes**: Echocardiography, infarct size, fibrosis quantification
- **Genetic models**: Epicardial-specific receptor knockout, inducible overexpression

---

## Tools and Resources

### Software Tools

| Category | Tools |
|----------|-------|
| **scRNA-seq** | [Seurat](https://satijalab.org/seurat/), [Scanpy](https://scanpy.readthedocs.io/), [scVI](https://scvi-tools.org/) |
| **Cell-cell communication** | [NicheNet](https://github.com/saeyslab/nichenetr), [CellChat](https://github.com/sqjin/CellChat), [CellPhoneDB](https://www.cellphonedb.org/), [LIANA](https://github.com/saezlab/liana) |
| **Deep learning** | [Geneformer](https://huggingface.co/ctheodoris/Geneformer), [TxGNN](https://github.com/mims-harvard/TxGNN), [PyTorch Geometric](https://pytorch-geometric.readthedocs.io/) |
| **Protein design** | [RFdiffusion](https://github.com/RosettaCommons/RFdiffusion), [BindCraft](https://github.com/martinpacesa/BindCraft), [AlphaFold](https://github.com/deepmind/alphafold) |

### Databases

| Category | Resources |
|----------|-----------|
| **Ligand-receptor** | [CellTalkDB](http://tcm.zju.edu.cn/celltalkdb/), [CellPhoneDB](https://www.cellphonedb.org/), [NATMI](https://github.com/forrest-lab/NATMI) |
| **Cardiac atlases** | [Heart Cell Atlas](https://www.heartcellatlas.org/), [Kuppe et al. portal](https://www.singlecell.rwth-aachen.de/) |
| **Drug-target** | [DrugBank](https://go.drugbank.com/), [ChEMBL](https://www.ebi.ac.uk/chembl/), [BindingDB](https://www.bindingdb.org/) |

---

## Key References

### Spatial/Single-Cell MI Atlases

1. Kuppe C, et al. (2022) "Spatial multi-omic map of human myocardial infarction" *Nature* 608:766-777. [DOI](https://doi.org/10.1038/s41586-022-05060-x)

2. Farbehi N, et al. (2019) "Single-cell expression profiling reveals dynamic flux of cardiac stromal, vascular and immune cells in health and injury" *eLife* 8:e43882. [DOI](https://doi.org/10.7554/eLife.43882)

3. Molenaar B, et al. (2021) "Single-cell transcriptomics following ischemic injury identifies a role for B2M in cardiac repair" *Commun Biol* 4:146. [DOI](https://doi.org/10.1038/s42003-020-01636-3)

4. Yamada S, et al. (2022) "Spatiotemporal transcriptome analysis reveals critical roles for mechano-sensing genes at the border zone in remodeling after MI" *Nat Cardiovasc Res* 1:1039-1055. [DOI](https://doi.org/10.1038/s44161-022-00140-7)

### AI/ML for Target Discovery

5. Theodoris CV, et al. (2023) "Transfer learning enables predictions in network biology" *Nature* 618:616-624. [DOI](https://doi.org/10.1038/s41586-023-06139-9) — *Geneformer*

6. Huang K, et al. (2024) "TxGNN: A foundation model for clinician-centered drug repurposing" *Nat Med*. [DOI](https://doi.org/10.1038/s41591-024-03233-x)

7. Singh R, et al. (2023) "Contrastive learning in protein language space predicts interactions between drugs and protein targets" *PNAS* 120:e2220778120. [DOI](https://doi.org/10.1073/pnas.2220778120) — *ConPLex*

8. Pun FW, et al. (2024) "Drug target prediction through deep learning functional representation of gene signatures" *Nat Commun* 15:1853. [DOI](https://doi.org/10.1038/s41467-024-46089-y) — *FRoGS*

### Protein Design for Therapeutics

9. Watson JL, et al. (2023) "De novo design of protein structure and function with RFdiffusion" *Nature* 620:1089-1100. [DOI](https://doi.org/10.1038/s41586-023-06415-8)

10. Pacesa M, et al. (2025) "BindCraft: One-shot design of functional protein binders" *Nature*. [DOI](https://doi.org/10.1038/s41586-025-09429-6)

### Cell-Cell Communication Methods

11. Browaeys R, et al. (2020) "NicheNet: modeling intercellular communication by linking ligands to target genes" *Nat Methods* 17:159-162. [DOI](https://doi.org/10.1038/s41592-019-0667-5)

12. Jin S, et al. (2021) "Inference and analysis of cell-cell communication using CellChat" *Nat Commun* 12:1088. [DOI](https://doi.org/10.1038/s41467-021-21246-9)

### Epicardial Biology

13. Cao J, Poss KD (2018) "The epicardium as a hub for heart regeneration" *Nat Rev Cardiol* 15:631-647. [DOI](https://doi.org/10.1038/s41569-018-0046-4)

14. Wang S, et al. (2015) "Epicardial regeneration is guided by cardiac outflow tract and Hedgehog signalling" *Nature* 522:226-230. [DOI](https://doi.org/10.1038/nature14325)

---

## Checklist

- [ ] Download Kuppe et al. or Farbehi et al. dataset
- [ ] Process and identify epicardial populations
- [ ] Define proliferation and EMT signatures
- [ ] Run NicheNet with cardiomyocytes as sender, epicardium as receiver
- [ ] Verify FGF10 ranks in top ligands (positive control)
- [ ] Map top ligands to epicardial receptors
- [ ] Build ML model with feature engineering
- [ ] Validate with leave-one-out CV
- [ ] Prioritize novel receptor candidates
- [ ] Design experimental validation strategy

---

## Notes

- **FGF10/FGFR2b** should serve as positive control throughout pipeline
- Cardiomyocytes are **SENDER** cells; epicardial cells are **RECEIVER** cells
- Start with NicheNet analysis as the most accessible entry point
- GPU resources required for Geneformer fine-tuning
- Experimental validation essential before clinical translation
