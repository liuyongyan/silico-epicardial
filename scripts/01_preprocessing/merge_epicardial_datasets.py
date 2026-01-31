#!/usr/bin/env python3
"""
Merge All Epicardial Datasets

Combines log-normalized epicardial cells from:
- PERIHEART (Linna-Kuosmanen)
- CAREBANK (Linna-Kuosmanen)
- Kuppe MI

Adds batch labels for downstream batch correction.
"""

import anndata as ad
import numpy as np
import pandas as pd
from pathlib import Path

PROJECT_DIR = Path(__file__).parent.parent.parent
PROCESSED_DIR = PROJECT_DIR / "data/processed"

def main():
    print("="*60)
    print("Merging Epicardial Datasets")
    print("="*60)

    datasets = []

    # Load PERIHEART
    print("\nLoading PERIHEART...")
    adata_ph = ad.read_h5ad(PROCESSED_DIR / "epicardial_periheart_log.h5ad")
    adata_ph.obs['dataset'] = 'PERIHEART'
    adata_ph.obs['source'] = 'Linna-Kuosmanen'
    print(f"  Cells: {adata_ph.n_obs:,}, Genes: {adata_ph.n_vars:,}")
    datasets.append(adata_ph)

    # Load CAREBANK
    print("\nLoading CAREBANK...")
    adata_cb = ad.read_h5ad(PROCESSED_DIR / "epicardial_carebank_log.h5ad")
    adata_cb.obs['dataset'] = 'CAREBANK'
    adata_cb.obs['source'] = 'Linna-Kuosmanen'
    print(f"  Cells: {adata_cb.n_obs:,}, Genes: {adata_cb.n_vars:,}")
    datasets.append(adata_cb)

    # Load Kuppe
    print("\nLoading Kuppe...")
    adata_kuppe = ad.read_h5ad(PROCESSED_DIR / "epicardial_kuppe.h5ad")
    adata_kuppe.obs['dataset'] = 'Kuppe_MI'
    adata_kuppe.obs['source'] = 'Kuppe'
    print(f"  Cells: {adata_kuppe.n_obs:,}, Genes: {adata_kuppe.n_vars:,}")
    datasets.append(adata_kuppe)

    # Find common genes
    print("\n" + "="*60)
    print("Finding common genes...")
    common_genes = set(datasets[0].var_names)
    for adata in datasets[1:]:
        common_genes = common_genes.intersection(set(adata.var_names))
    common_genes = sorted(list(common_genes))
    print(f"  Common genes: {len(common_genes):,}")

    # Subset to common genes
    print("\nSubsetting to common genes...")
    for i, adata in enumerate(datasets):
        datasets[i] = adata[:, common_genes].copy()

    # Concatenate
    print("\nConcatenating datasets...")
    adata_merged = ad.concat(datasets, join='outer', label='batch', keys=['PERIHEART', 'CAREBANK', 'Kuppe_MI'])

    # Clean up obs index
    adata_merged.obs_names_make_unique()

    print(f"\nMerged dataset:")
    print(f"  Total cells: {adata_merged.n_obs:,}")
    print(f"  Total genes: {adata_merged.n_vars:,}")

    # Summary by dataset
    print("\nCells per dataset:")
    for ds in ['PERIHEART', 'CAREBANK', 'Kuppe_MI']:
        n = (adata_merged.obs['dataset'] == ds).sum()
        print(f"  {ds}: {n:,}")

    # Save
    output_path = PROCESSED_DIR / "epicardial_all_merged.h5ad"
    print(f"\nSaving to {output_path}...")
    adata_merged.write_h5ad(output_path)

    print("\n" + "="*60)
    print("Done!")
    print("="*60)

if __name__ == "__main__":
    main()
