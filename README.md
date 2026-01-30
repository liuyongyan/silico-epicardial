# AI-Driven Epicardial Target Discovery in Myocardial Infarction

> **Purpose**: Identify potential therapeutic targets on epicardial cells that promote proliferation and EMT for cardioprotection after myocardial infarction

**Version**: 2.2 (January 2026)

## Table of Contents

- [Overview and Rationale](#overview-and-rationale)
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

### Software Stack (Python)

| Category | Tools |
|----------|-------|
| **Data handling** | anndata, scanpy |
| **Cell communication** | LIANA, CellPhoneDB |
| **Deep learning** | Geneformer, PyTorch |
| **Visualization** | matplotlib, seaborn |

### Requirements

- Python 3.8+
- GPU with 16GB+ VRAM (for Geneformer)
- ~10GB disk space for raw data

---

## Key References

1. Kuppe C, et al. (2022) "Spatial multi-omic map of human myocardial infarction" *Nature* 608:766-777
2. Linna-Kuosmanen S, et al. (2024) "Single-nucleus RNA sequencing of human heart" *Cell Reports Medicine*
3. Theodoris CV, et al. (2023) "Transfer learning enables predictions in network biology" *Nature* 618:616-624 (Geneformer)
4. Cao J, Poss KD (2018) "The epicardium as a hub for heart regeneration" *Nat Rev Cardiol* 15:631-647

---

## Notes

- **FGF10/FGFR2b** serves as positive control throughout pipeline
- Cardiomyocytes are **SENDER** cells; epicardial cells are **RECEIVER** cells
- Mesothelial cells in Linna-Kuosmanen = epicardial cells
- Data format: `.h5ad` (AnnData/HDF5), use `backed='r'` mode for large files
