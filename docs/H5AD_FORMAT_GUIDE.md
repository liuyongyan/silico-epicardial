# h5ad File Format Guide

## What is h5ad?

**h5ad** = HDF5 + AnnData, the most commonly used data format in single-cell omics.

- **HDF5** (Hierarchical Data Format 5): An efficient binary file format supporting compression and hierarchical storage
- **AnnData**: The core data structure in Python's `scanpy` ecosystem, designed for single-cell data

## AnnData Data Structure

```
AnnData Object
│
├── X              # Expression matrix (cells × genes)
│                  # Can be dense numpy array or sparse matrix
│
├── obs            # Cell metadata (DataFrame, n_obs × n_obs_keys)
│   ├── cell_type  # Cell type annotations
│   ├── sample     # Sample origin
│   ├── disease    # Disease state
│   └── ...        # Other cell-level annotations
│
├── var            # Gene metadata (DataFrame, n_vars × n_var_keys)
│   ├── gene_name  # Gene symbol
│   ├── ensembl_id # Ensembl ID
│   └── ...        # Other gene-level annotations
│
├── uns            # Unstructured data (dict)
│   ├── neighbors  # KNN graph parameters
│   ├── umap       # UMAP parameters
│   └── ...        # Analysis parameters, color schemes, etc.
│
├── obsm           # Cell embedding matrices (dict of arrays)
│   ├── X_pca      # PCA coordinates
│   ├── X_umap     # UMAP coordinates
│   └── X_tsne     # t-SNE coordinates
│
├── varm           # Gene embedding matrices (dict of arrays)
│   └── PCs        # Gene loadings on PCs
│
├── layers         # Other expression matrices (dict)
│   ├── counts     # Raw counts
│   ├── normalized # Normalized data
│   └── ...
│
└── obsp / varp    # Cell/gene pairwise relationship matrices
    └── connectivities  # KNN connectivity matrix
```

## How to Read h5ad Files

### Python (Recommended)

```python
import anndata as ad

# Method 1: Load fully into memory (small datasets)
adata = ad.read_h5ad("data.h5ad")

# Method 2: Backed mode (large datasets, memory efficient)
adata = ad.read_h5ad("data.h5ad", backed='r')  # read-only
adata = ad.read_h5ad("data.h5ad", backed='r+') # read-write

# View basic info
print(f"Cells: {adata.n_obs}")
print(f"Genes: {adata.n_vars}")

# View cell metadata
print(adata.obs.columns)  # all column names
print(adata.obs['cell_type'].value_counts())  # cell type distribution

# View gene metadata
print(adata.var.columns)
print(adata.var_names[:10])  # first 10 gene IDs

# Get expression matrix
X = adata.X  # sparse matrix
X_dense = adata.X.toarray()  # convert to dense (caution: memory!)

# Get expression for specific gene
gene_idx = list(adata.var_names).index('ENSG00000184937')  # WT1
wt1_expr = adata.X[:, gene_idx]
```

### R (via SeuratDisk)

```r
library(SeuratDisk)
library(Seurat)

# First convert to h5seurat format
Convert("data.h5ad", dest = "h5seurat", overwrite = TRUE)

# Then load as Seurat object
seurat_obj <- LoadH5Seurat("data.h5seurat")

# Or use anndata package to read directly
library(anndata)
adata <- read_h5ad("data.h5ad")
```

## CellxGene Standardized Metadata

h5ad files downloaded from CellxGene follow a standardized schema with these columns:

### obs (Cell Metadata)

| Column | Description | Example |
|--------|-------------|---------|
| `cell_type` | Cell type (human readable) | "fibroblast" |
| `cell_type_ontology_term_id` | Cell Ontology ID | "CL:0000057" |
| `disease` | Disease state | "myocardial infarction" |
| `disease_ontology_term_id` | MONDO ID | "MONDO:0005068" |
| `tissue` | Tissue origin | "heart left ventricle" |
| `tissue_ontology_term_id` | UBERON ID | "UBERON:0002084" |
| `assay` | Sequencing technology | "10x 3' v3" |
| `sex` | Sex | "male" / "female" |
| `development_stage` | Development stage | "adult" |
| `self_reported_ethnicity` | Ethnicity | "European" |
| `suspension_type` | Suspension type | "nucleus" / "cell" |

### var (Gene Metadata)

| Column | Description | Example |
|--------|-------------|---------|
| `feature_name` | Gene symbol | "WT1" |
| `feature_type` | Feature type | "gene" |
| `feature_biotype` | Biotype | "protein_coding" |
| `feature_reference` | Reference genome | "NCBITaxon:9606" |

**Note**: var_names (index) use Ensembl IDs, e.g., "ENSG00000184937"

## Kuppe 2022 Data Structure Example

```
Kuppe MI Heart Data
├── n_obs: 191,795 cells
├── n_vars: 28,380 genes
│
├── obs columns (30):
│   ├── sample              # Sample ID
│   ├── donor_id            # Donor ID
│   ├── patient_group       # myogenic/ischemic/fibrotic
│   ├── disease             # myocardial infarction / normal
│   ├── cell_type           # CellxGene standardized cell type
│   ├── cell_type_original  # Original paper annotations
│   └── ...
│
├── var columns (6):
│   ├── feature_name        # Gene symbol
│   ├── feature_biotype     # protein_coding, etc.
│   └── ...
│
└── var_names: Ensembl IDs (ENSG...)
```

## Memory Management Tips

| Data Size | Recommended Method |
|-----------|-------------------|
| < 500 MB | Direct load: `ad.read_h5ad()` |
| 500 MB - 2 GB | Backed mode: `ad.read_h5ad(backed='r')` |
| > 2 GB | Chunk processing or use Dask |

## Common Operations

```python
# Subset cells
fibroblasts = adata[adata.obs['cell_type'] == 'fibroblast']

# Subset genes
markers = adata[:, ['ENSG00000184937', 'ENSG00000112837']]

# Filter by condition
mi_cells = adata[adata.obs['disease'] == 'myocardial infarction']

# Save processed data
adata.write('processed.h5ad')
```

## Reference Resources

- [AnnData Documentation](https://anndata.readthedocs.io/)
- [Scanpy Tutorials](https://scanpy.readthedocs.io/)
- [CellxGene Schema](https://github.com/chanzuckerberg/single-cell-curation)
