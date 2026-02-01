#!/usr/bin/env python3
"""
Extract Epicardial Cells from Kuppe Dataset

Uses marker expression score to identify epicardial cells since
this dataset doesn't have explicit mesothelial annotations.
"""

import anndata as ad
import numpy as np
import pandas as pd
from pathlib import Path
from scipy import sparse
import gc

PROJECT_DIR = Path(__file__).parent.parent.parent
RAW_DIR = PROJECT_DIR / "data/raw"
PROCESSED_DIR = PROJECT_DIR / "data/processed"

# Epicardial markers (Ensembl IDs)
MARKERS = {
    'WT1': 'ENSG00000184937',
    'TBX18': 'ENSG00000112837',
    'ALDH1A2': 'ENSG00000128918',
    'UPK3B': 'ENSG00000243566',
}

# Disease conditions to keep (MI + Normal only)
KEEP_DISEASES = ['myocardial infarction', 'normal']


def main(filter_disease=True):
    print("\n#" + "="*58 + "#")
    print("# Kuppe Epicardial Cell Extraction")
    print("#" + "="*58 + "#")

    kuppe_file = RAW_DIR / "kuppe" / "54d24dbe-d39a-4844-bb21-07b5f4e173ad.h5ad"

    # Load in backed mode
    print(f"\nLoading {kuppe_file.name}...")
    adata = ad.read_h5ad(kuppe_file, backed='r')
    print(f"Total cells: {adata.n_obs:,}")
    print(f"Total genes: {adata.n_vars:,}")

    # Show disease distribution
    print("\nDisease distribution:")
    print(adata.obs['disease'].value_counts().to_string())

    # Create disease mask
    if filter_disease:
        print(f"\nFiltering to diseases: {KEEP_DISEASES}")
        disease_mask = adata.obs['disease'].isin(KEEP_DISEASES).values
        print(f"Cells with MI/Normal: {disease_mask.sum():,} / {adata.n_obs:,}")
    else:
        disease_mask = np.ones(adata.n_obs, dtype=bool)

    # Find marker gene indices
    print("\nChecking markers:")
    marker_idx = []
    for name, ensembl in MARKERS.items():
        if ensembl in adata.var_names:
            idx = list(adata.var_names).index(ensembl)
            marker_idx.append(idx)
            print(f"  + {name}: present (idx={idx})")
        else:
            print(f"  - {name}: NOT FOUND")

    # Calculate epicardial score in chunks
    print(f"\nCalculating epicardial score using {len(marker_idx)} markers...")
    n_cells = adata.n_obs
    scores = np.zeros(n_cells, dtype=np.float32)
    chunk_size = 20000

    for i in range(0, n_cells, chunk_size):
        end = min(i + chunk_size, n_cells)
        X_chunk = adata.X[i:end, marker_idx]
        if sparse.issparse(X_chunk):
            X_chunk = X_chunk.toarray()
        scores[i:end] = X_chunk.mean(axis=1)

        if end % 50000 == 0 or end == n_cells:
            print(f"  {end:,} / {n_cells:,}")

    # Score stats
    print(f"\nScore distribution:")
    print(f"  Min: {scores.min():.4f}")
    print(f"  Max: {scores.max():.4f}")
    print(f"  Mean: {scores.mean():.4f}")
    print(f"  Std: {scores.std():.4f}")

    # Use 95th percentile as threshold
    threshold = np.percentile(scores, 95)
    print(f"  95th percentile: {threshold:.4f}")

    # Get high-scoring cells (apply disease filter)
    high_mask = (scores > threshold) & disease_mask
    n_high = high_mask.sum()
    print(f"\nHigh-scoring cells (>95th, disease-filtered): {n_high:,}")

    # Also get adipocytes (apply disease filter)
    adipo_mask = (adata.obs['cell_type_original'] == 'Adipocyte').values & disease_mask
    n_adipo = adipo_mask.sum()
    print(f"Adipocytes (disease-filtered): {n_adipo:,}")

    # Combine
    combined_mask = high_mask | adipo_mask
    n_combined = combined_mask.sum()
    print(f"Combined: {n_combined:,}")

    # Show disease breakdown
    if filter_disease:
        print("\nDisease breakdown of selected cells:")
        selected_obs = adata.obs.iloc[combined_mask]
        for disease in KEEP_DISEASES:
            count = (selected_obs['disease'] == disease).sum()
            print(f"  {disease}: {count:,}")

    # Cell type distribution
    print("\nCell type distribution of selected cells:")
    selected_ct = adata.obs.iloc[combined_mask]['cell_type_original'].value_counts()
    for ct, count in selected_ct.head(10).items():
        print(f"  {ct}: {count:,}")

    # Get indices
    keep_idx = np.where(combined_mask)[0]

    # Extract data
    print(f"\nExtracting {n_combined:,} cells...")
    obs_subset = adata.obs.iloc[keep_idx].copy()
    var_df = adata.var.copy()

    # Read X in batches
    batch_size = 1000
    X_parts = []
    for i in range(0, n_combined, batch_size):
        end = min(i + batch_size, n_combined)
        batch_idx = keep_idx[i:end]
        X_batch = adata.X[batch_idx, :]
        if sparse.issparse(X_batch):
            X_parts.append(X_batch.copy())
        else:
            X_parts.append(sparse.csr_matrix(X_batch))

        if end % 5000 == 0 or end == n_combined:
            print(f"  {end:,} / {n_combined:,}")

    print("Stacking...")
    X_combined = sparse.vstack(X_parts)
    del X_parts
    gc.collect()

    # Create AnnData
    print("Creating AnnData...")
    adata_epi = ad.AnnData(
        X=X_combined,
        obs=obs_subset.reset_index(drop=True),
        var=var_df
    )

    # Add epicardial score
    adata_epi.obs['epicardial_score'] = scores[combined_mask]
    adata_epi.obs['source_dataset'] = 'Kuppe_MI'
    adata_epi.obs['source_study'] = 'Kuppe_2022'

    # Save
    output_file = PROCESSED_DIR / "epicardial_kuppe.h5ad"
    print(f"\nSaving to {output_file}...")
    adata_epi.write(output_file)
    print(f"Saved: {output_file.stat().st_size / 1e6:.1f} MB")

    print("\n" + "="*60)
    print(f"Extracted {adata_epi.n_obs:,} epicardial cells from Kuppe")
    print("="*60)


if __name__ == "__main__":
    main()
