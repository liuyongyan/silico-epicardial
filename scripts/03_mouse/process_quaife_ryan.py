#!/usr/bin/env python3
"""
Process Quaife-Ryan et al. 2021 mouse epicardial dataset (E-MTAB-10035)
Goal: Identify epicardial cells, classify activation states, analyze FGF10/FGFR2
"""

import scanpy as sc
import pandas as pd
import numpy as np
from scipy.io import mmread
from scipy.sparse import csr_matrix
import matplotlib.pyplot as plt
import os
import mygene

# Paths
DATA_DIR = "/Users/yliu/Desktop/Columbia - Biostatistics/Cheng Lab/Silico Epicardial/data/mouse/quaife_ryan_2021"
OUTPUT_DIR = "/Users/yliu/Desktop/Columbia - Biostatistics/Cheng Lab/Silico Epicardial/data/processed"
RESULTS_DIR = "/Users/yliu/Desktop/Columbia - Biostatistics/Cheng Lab/Silico Epicardial/results/mouse"

os.makedirs(RESULTS_DIR, exist_ok=True)

print("=" * 60)
print("Loading Quaife-Ryan 2021 Mouse Epicardial Dataset")
print("=" * 60)

# Load matrix
print("\n1. Loading expression matrix (this may take a few minutes)...")
mtx_path = f"{DATA_DIR}/E-MTAB-10035.aggregated_filtered_normalised_counts.mtx"
matrix = mmread(mtx_path).T.tocsr()  # Transpose: cells x genes
print(f"   Matrix shape: {matrix.shape}")

# Load gene names (rows)
print("\n2. Loading gene names...")
genes = pd.read_csv(f"{DATA_DIR}/E-MTAB-10035.aggregated_filtered_normalised_counts.mtx_rows",
                    header=None, sep='\t', names=['gene_id', 'gene_id2'])
print(f"   Genes: {len(genes)}")
print(f"   Example gene IDs: {genes['gene_id'].head().tolist()}")

# Load cell barcodes (cols)
print("\n3. Loading cell barcodes...")
cells = pd.read_csv(f"{DATA_DIR}/E-MTAB-10035.aggregated_filtered_normalised_counts.mtx_cols",
                    header=None, names=['cell_id'])
print(f"   Cells: {len(cells)}")

# Convert Ensembl IDs to gene symbols
print("\n4. Converting Ensembl IDs to gene symbols...")
ensembl_ids = genes['gene_id'].str.replace(r'\.\d+$', '', regex=True).tolist()  # Remove version numbers

# Check for cached mapping
cache_file = f"{DATA_DIR}/ensembl_to_symbol_cache.csv"
if os.path.exists(cache_file):
    print("   Loading cached gene mapping...")
    gene_mapping = pd.read_csv(cache_file)
else:
    print("   Querying mygene for gene symbols (this may take a minute)...")
    mg = mygene.MyGeneInfo()
    result = mg.querymany(ensembl_ids, scopes='ensembl.gene', fields='symbol', species='mouse', verbose=False)

    gene_mapping = pd.DataFrame([
        {'ensembl_id': r.get('query'), 'symbol': r.get('symbol', r.get('query'))}
        for r in result
    ])
    gene_mapping.to_csv(cache_file, index=False)
    print(f"   Cached gene mapping to: {cache_file}")

# Create mapping dict
ensembl_to_symbol = dict(zip(gene_mapping['ensembl_id'], gene_mapping['symbol']))

# Map gene names
gene_symbols = [ensembl_to_symbol.get(eid.split('.')[0], eid) for eid in genes['gene_id']]
print(f"   Mapped {sum(1 for s in gene_symbols if not s.startswith('ENSMUSG'))}/{len(gene_symbols)} genes to symbols")

# Create AnnData
print("\n5. Creating AnnData object...")
adata = sc.AnnData(X=matrix)
adata.obs_names = cells['cell_id'].values
adata.var_names = gene_symbols
adata.var['ensembl_id'] = genes['gene_id'].values
adata.var_names_make_unique()
print(f"   AnnData: {adata.shape[0]} cells x {adata.shape[1]} genes")

# Load clustering info
print("\n6. Loading clustering annotations...")
clustering = pd.read_csv(f"{DATA_DIR}/clustering.tsv", sep='\t')
print(f"   Clustering shape: {clustering.shape}")

# Find the selected K row
selected_row = clustering[clustering['sel.K'] == True]
if len(selected_row) == 0:
    selected_row = clustering[clustering['K'] == 12]

selected_k = selected_row['K'].values[0]
print(f"   Using K={selected_k} clustering")

# Extract cluster assignments
cell_barcodes = clustering.columns[2:]
cluster_assignments = selected_row.iloc[0, 2:].values
cluster_map = dict(zip(cell_barcodes, cluster_assignments.astype(int)))
adata.obs['cluster'] = adata.obs_names.map(cluster_map)
print(f"   Clusters: {adata.obs['cluster'].nunique()} unique")

# Load experiment design for condition info
print("\n7. Loading experiment design...")
exp_design = pd.read_csv(f"{DATA_DIR}/experiment_design.tsv", sep='\t')
print(f"   Experiment design columns of interest:")

# Check for condition/disease columns
condition_cols = [c for c in exp_design.columns if any(x in c.lower() for x in ['disease', 'condition', 'stimulus', 'infar'])]
print(f"   Condition columns: {condition_cols}")

# Add experiment metadata to adata
exp_design.index = exp_design['Assay']
for col in ['Sample Characteristic[disease]', 'Sample Characteristic[stimulus]', 'Sample Characteristic[individual]']:
    if col in exp_design.columns:
        adata.obs[col.split('[')[1].rstrip(']')] = adata.obs_names.map(
            lambda x: exp_design.loc[x, col] if x in exp_design.index else None
        )
        print(f"   Added: {col.split('[')[1].rstrip(']')}")

# Basic QC
print("\n8. Computing QC metrics...")
sc.pp.calculate_qc_metrics(adata, inplace=True)
print(f"   Mean genes per cell: {adata.obs['n_genes_by_counts'].mean():.0f}")
print(f"   Mean counts per cell: {adata.obs['total_counts'].mean():.0f}")

# Check for epicardial markers
print("\n9. Checking key markers...")
epicardial_markers = ['Wt1', 'Upk3b', 'Msln', 'Tcf21', 'Tbx18']
found_epi = [g for g in epicardial_markers if g in adata.var_names]
print(f"   Epicardial markers: {found_epi}")

# FGF genes
fgf_genes = ['Fgf10', 'Fgfr1', 'Fgfr2', 'Fgfr3', 'Fgf2', 'Fgf7', 'Fgf1']
found_fgf = [g for g in fgf_genes if g in adata.var_names]
print(f"   FGF genes: {found_fgf}")

# Activation markers
activation_markers = ['Vim', 'Fn1', 'Col1a1', 'Col3a1', 'Postn', 'Acta2', 'Snai1', 'Snai2', 'Twist1']
found_act = [g for g in activation_markers if g in adata.var_names]
print(f"   Activation markers: {found_act}")

# Quiescent markers
quiescent_markers = ['Cdh1', 'Cldn1', 'Cldn3', 'Tjp1', 'Ocln', 'Krt19']
found_qui = [g for g in quiescent_markers if g in adata.var_names]
print(f"   Quiescent markers: {found_qui}")

# Save intermediate result
print("\n10. Saving AnnData with gene symbols...")
adata.write(f"{OUTPUT_DIR}/mouse_quaife_ryan_raw.h5ad")
print(f"    Saved to: {OUTPUT_DIR}/mouse_quaife_ryan_raw.h5ad")

print("\n" + "=" * 60)
print("Initial loading complete!")
print("=" * 60)
print(f"\nDataset summary:")
print(f"  - Total cells: {adata.shape[0]:,}")
print(f"  - Total genes: {adata.shape[1]:,}")
print(f"  - Clusters: {adata.obs['cluster'].nunique()}")
print(f"  - Key markers found: epicardial={len(found_epi)}, FGF={len(found_fgf)}, activation={len(found_act)}")
