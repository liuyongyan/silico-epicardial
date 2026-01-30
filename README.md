# AI-Driven Epicardial Target Discovery in Myocardial Infarction

> **Purpose**: Identify potential therapeutic targets on epicardial cells that promote proliferation and EMT for cardioprotection after myocardial infarction

**Version**: 2.3 (January 2026)

## Table of Contents

- [Overview and Rationale](#overview-and-rationale)
- [Environment Setup](#environment-setup)
- [Data Sources](#data-sources)
- [Current Progress](#current-progress)
- [Project Structure](#project-structure)
- [Pipeline Overview](#pipeline-overview)
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

## Environment Setup

### Why Python Instead of R?

This project uses **Python** instead of R/Seurat for several reasons:
- Data is in `.h5ad` format (AnnData), native to Python ecosystem
- No format conversion needed (R requires `SeuratDisk` to read h5ad)
- Better integration with deep learning tools (Geneformer, PyTorch)
- Unified workflow from preprocessing to model training

### Python Dependencies

**Core packages:**

```bash
# Create virtual environment (recommended)
python -m venv venv
source venv/bin/activate  # Linux/Mac
# or: venv\Scripts\activate  # Windows

# Install core dependencies
pip install anndata>=0.10.0
pip install scanpy>=1.10.0
pip install pandas numpy scipy
pip install matplotlib seaborn

# Cell-cell communication
pip install liana
pip install cellphonedb

# Deep learning (for Geneformer)
pip install torch transformers datasets
pip install geneformer  # or clone from GitHub
```

**Full requirements:**

| Package | Version | Purpose |
|---------|---------|---------|
| `anndata` | >=0.10.0 | h5ad file I/O, AnnData structure |
| `scanpy` | >=1.10.0 | Single-cell analysis, clustering, visualization |
| `pandas` | >=2.0.0 | Data manipulation |
| `numpy` | >=1.24.0 | Numerical computing |
| `scipy` | >=1.10.0 | Sparse matrices, statistics |
| `matplotlib` | >=3.7.0 | Plotting |
| `seaborn` | >=0.12.0 | Statistical visualization |
| `liana` | >=1.0.0 | Ligand-receptor analysis |
| `torch` | >=2.0.0 | Deep learning framework |
| `transformers` | >=4.30.0 | Geneformer model loading |

### Hardware Requirements

| Component | Minimum | Recommended |
|-----------|---------|-------------|
| RAM | 16 GB | 32 GB |
| Disk | 10 GB | 20 GB |
| GPU | - | 16GB VRAM (for Geneformer) |

---

## Data Sources

All data downloaded from [CellxGene](https://cellxgene.cziscience.com/gene-expression) in `.h5ad` format.

### Kuppe et al. 2022 (Nature)

| Property | Value |
|----------|-------|
| **File** | `data/raw/kuppe/54d24dbe-d39a-4844-bb21-07b5f4e173ad.h5ad` |
| **Size** | 1.8 GB |
| **Tissue** | Heart left ventricle |
| **Cells** | 191,795 |
| **Genes** | 28,380 |
| **Conditions** | MI (various zones) vs normal |
| **Epicardial** | ~523 adipocytes (epicardial fat); true mesothelial cells likely in Fibroblast population |

**Key findings from exploration:**
- Disease states: myocardial infarction (168,044 cells), normal (23,751 cells)
- Patient groups: myogenic, ischemic, fibrotic
- Epicardial markers present: WT1, TBX18, ALDH1A2, UPK3B (MSLN not found in this dataset)

### Linna-Kuosmanen et al. 2024 (Cell Reports Medicine)

Two cohorts from Finnish cardiac surgery patients (right atrium tissue):

| Dataset | File | Cells | Mesothelial (Epicardial) |
|---------|------|-------|--------------------------|
| **PERIHEART** | `77df1b04-4905-4dac-b9aa-6243f10201ae.h5ad` | 392,819 | 28,851 (7.3%) |
| **CAREBANK** | `8e1afde5-6a0a-4ed7-ba0b-1bcdc335ddfe.h5ad` | 221,756 | 20,493 (9.2%) |

**Key findings:**
- **Rich epicardial cell annotations**: "mesothelial cell" = epicardial cells (~49,000 total)
- All 6 epicardial markers present including MSLN
- 35,477 genes per dataset
- Multiple cell types: cardiomyocyte, fibroblast, endothelial, macrophage, etc.

> **Important**: Linna-Kuosmanen data provides the best source of epicardial cells for this project due to explicit mesothelial cell annotations.

---

## Current Progress

### Completed

- [x] Download all datasets from CellxGene
  - Kuppe 2022 MI heart data
  - Linna-Kuosmanen PERIHEART cohort
  - Linna-Kuosmanen CAREBANK cohort
- [x] Explore data structure (cells, genes, metadata)
- [x] Identify epicardial cell populations in each dataset
- [x] Set up Python-based analysis workflow (scanpy/anndata)
- [x] Document h5ad format and CellxGene schema

### Next Steps

- [ ] Extract and merge epicardial cells from all datasets
- [ ] Define proliferation and EMT gene signatures
- [ ] Run ligand-receptor analysis (LIANA/CellPhoneDB)
- [ ] Fine-tune Geneformer for epicardial activation
- [ ] In silico perturbation to rank receptor candidates
- [ ] Validate FGF10/FGFR2b as positive control

---

## Project Structure

```
epicardial-target-discovery/
├── README.md
├── data/
│   ├── raw/
│   │   ├── kuppe/
│   │   │   └── 54d24dbe-*.h5ad          # Kuppe MI data (1.8GB)
│   │   └── linna_kuosmanen/
│   │       ├── 77df1b04-*.h5ad          # PERIHEART (2.6GB)
│   │       └── 8e1afde5-*.h5ad          # CAREBANK (1.5GB)
│   ├── processed/
│   │   ├── kuppe_cell_metadata.csv      # 191,795 cells
│   │   └── kuppe_gene_metadata.csv      # 28,380 genes
│   └── references/
├── scripts/
│   └── 01_preprocessing/
│       └── analyze_kuppe.py             # Data exploration script
├── results/
│   └── tables/
│       ├── kuppe_celltype_summary.csv
│       └── kuppe_sample_metadata.csv
├── docs/
│   └── H5AD_FORMAT_GUIDE.md             # h5ad format reference
└── papers/
```

---

## Pipeline Overview

### Phase 1: Data Preprocessing
- Load h5ad files with `anndata` (backed mode for large files)
- Extract epicardial/mesothelial cells using annotations
- Identify cell states by marker expression

### Phase 2: Define Target Phenotypes
- **Proliferation**: MKI67, TOP2A, PCNA, CDK1, CCNB1/2, MCM2/6, AURKA/B
- **EMT (up)**: SNAI1/2, TWIST1, ZEB1/2, VIM, CDH2, FN1
- **EMT (down)**: CDH1, CLDN1, TJP1, OCLN

### Phase 3: Ligand-Receptor Analysis
- Use LIANA or CellPhoneDB to identify signaling
- Focus on sender→epicardium interactions
- Verify FGF10→FGFR2 appears in top hits

### Phase 4: Geneformer Fine-tuning
- Tokenize epicardial transcriptomes
- Fine-tune for quiescent vs activated classification
- Run in silico perturbation on candidate receptors

### Phase 5: Target Prioritization
- Rank receptors by perturbation effect
- Validate FGFR2 as positive control
- Identify novel therapeutic targets

---

## Epicardial Markers

| Marker | Ensembl ID | Function |
|--------|------------|----------|
| WT1 | ENSG00000184937 | Master regulator |
| TBX18 | ENSG00000112837 | Transcription factor |
| ALDH1A2 | ENSG00000128918 | Retinoic acid synthesis |
| UPK3B | ENSG00000243566 | Mesothelial marker |
| MSLN | ENSG00000102854 | Mesothelin |

> **Note**: TCF21 excluded due to low expression in adult human epicardium

---

## Tools and Resources

### Software Stack

| Category | Package | Documentation |
|----------|---------|---------------|
| **Data I/O** | [anndata](https://anndata.readthedocs.io/) | h5ad format, AnnData structure |
| **Single-cell analysis** | [scanpy](https://scanpy.readthedocs.io/) | Preprocessing, clustering, DE analysis |
| **Cell communication** | [LIANA](https://liana-py.readthedocs.io/) | Multi-method ligand-receptor analysis |
| **Cell communication** | [CellPhoneDB](https://www.cellphonedb.org/) | Statistical framework for interactions |
| **Deep learning** | [Geneformer](https://huggingface.co/ctheodoris/Geneformer) | Foundation model for single-cell |
| **ML framework** | [PyTorch](https://pytorch.org/docs/) | Deep learning backend |
| **Transformers** | [HuggingFace](https://huggingface.co/docs/transformers/) | Model loading and fine-tuning |

### Databases

| Category | Resource | URL |
|----------|----------|-----|
| **Single-cell data** | CellxGene | https://cellxgene.cziscience.com/ |
| **Ligand-receptor** | CellTalkDB | http://tcm.zju.edu.cn/celltalkdb/ |
| **Ligand-receptor** | CellPhoneDB | https://www.cellphonedb.org/ |
| **Cardiac atlas** | Heart Cell Atlas | https://www.heartcellatlas.org/ |
| **Gene ontology** | Cell Ontology | https://obofoundry.org/ontology/cl.html |

### Key Tutorials

- [Scanpy Tutorial](https://scanpy-tutorials.readthedocs.io/en/latest/pbmc3k.html) - Basic single-cell analysis
- [AnnData Guide](https://anndata.readthedocs.io/en/latest/tutorials/notebooks/getting-started.html) - Working with h5ad files
- [LIANA Tutorial](https://liana-py.readthedocs.io/en/latest/notebooks/basic_usage.html) - Cell-cell communication
- [Geneformer Tutorial](https://huggingface.co/ctheodoris/Geneformer) - Fine-tuning and perturbation

---

## Key References

### Data Sources

1. Kuppe C, et al. (2022) "Spatial multi-omic map of human myocardial infarction" *Nature* 608:766-777. [DOI](https://doi.org/10.1038/s41586-022-05060-x)

2. Linna-Kuosmanen S, et al. (2024) "Single-nucleus RNA sequencing of human heart" *Cell Reports Medicine*

### Methods

3. Wolf FA, et al. (2018) "SCANPY: large-scale single-cell gene expression data analysis" *Genome Biology* 19:15. [DOI](https://doi.org/10.1186/s13059-017-1382-0)

4. Dimitrov D, et al. (2022) "Comparison of methods and resources for cell-cell communication inference from single-cell RNA-Seq data" *Nat Commun* 13:3224. [DOI](https://doi.org/10.1038/s41467-022-30755-0) - LIANA

5. Theodoris CV, et al. (2023) "Transfer learning enables predictions in network biology" *Nature* 618:616-624. [DOI](https://doi.org/10.1038/s41586-023-06139-9) - Geneformer

### Biology

6. Cao J, Poss KD (2018) "The epicardium as a hub for heart regeneration" *Nat Rev Cardiol* 15:631-647. [DOI](https://doi.org/10.1038/s41569-018-0046-4)

---

## Notes

- **FGF10/FGFR2b** serves as positive control throughout pipeline
- Cardiomyocytes are **SENDER** cells; epicardial cells are **RECEIVER** cells
- Mesothelial cells in Linna-Kuosmanen = epicardial cells
- Data format: `.h5ad` (AnnData/HDF5), use `backed='r'` mode for large files
- Python environment: recommend using virtual environment or conda
