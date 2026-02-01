#!/usr/bin/env python3
"""
Extract Epicardial Cells - Simple Memory-Efficient Version

Uses h5py directly to extract only mesothelial cells.
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
PROCESSED_DIR.mkdir(parents=True, exist_ok=True)

# Disease conditions to keep (MI + Normal only)
KEEP_DISEASES = ['myocardial infarction', 'normal']


def extract_mesothelial_simple(filepath, output_name, filter_disease=True):
    """Extract mesothelial cells using simple row-by-row approach."""
    print(f"\n{'='*60}")
    print(f"Processing: {filepath.name}")
    print(f"{'='*60}")

    # Load in backed mode
    adata = ad.read_h5ad(filepath, backed='r')
    print(f"Total cells: {adata.n_obs:,}")
    print(f"Total genes: {adata.n_vars:,}")

    # Filter by disease condition if requested
    if filter_disease:
        print(f"\nFiltering to diseases: {KEEP_DISEASES}")
        print("Disease distribution before filtering:")
        print(adata.obs['disease'].value_counts().to_string())

        disease_mask = adata.obs['disease'].isin(KEEP_DISEASES)
        print(f"\nCells with MI/Normal: {disease_mask.sum():,} / {adata.n_obs:,}")
    else:
        disease_mask = pd.Series([True] * adata.n_obs, index=adata.obs_names)

    # Find mesothelial cells (combined with disease filter)
    meso_mask = (adata.obs['cell_type'] == 'mesothelial cell') & disease_mask.values
    meso_indices = np.where(meso_mask)[0]
    n_meso = len(meso_indices)
    print(f"\nMesothelial cells (after disease filter): {n_meso:,}")

    if filter_disease and n_meso > 0:
        print("Disease breakdown of selected cells:")
        selected_diseases = adata.obs.iloc[meso_indices]['disease'].value_counts()
        for disease, count in selected_diseases.items():
            print(f"  {disease}: {count:,}")

    if n_meso == 0:
        return None

    # Extract obs and var
    obs_meso = adata.obs.iloc[meso_indices].copy()
    var_df = adata.var.copy()

    # Step 1: Check if normalization is needed BEFORE extraction
    # Sample a few cells to detect data scale
    print("\nChecking data scale...")
    sample_idx = meso_indices[:min(100, n_meso)]
    X_sample = adata.X[sample_idx, :]
    if sparse.issparse(X_sample):
        sample_max = X_sample.data.max() if len(X_sample.data) > 0 else 0
    else:
        sample_max = X_sample.max()

    needs_log1p = sample_max > 50  # Linear scale (CPM) needs transformation
    if needs_log1p:
        print(f"  Data is in CPM format (max={sample_max:.1f}), will apply log1p during extraction")
    else:
        print(f"  Data is already log-normalized (max={sample_max:.2f})")

    # Step 2: Extract and normalize in batches
    print(f"\nExtracting expression data...")
    batch_size = 1000
    X_parts = []

    for i in range(0, n_meso, batch_size):
        batch_end = min(i + batch_size, n_meso)
        batch_idx = meso_indices[i:batch_end]

        # Read batch
        X_batch = adata.X[batch_idx, :]
        if not sparse.issparse(X_batch):
            X_batch = sparse.csr_matrix(X_batch)
        else:
            X_batch = X_batch.copy()

        # Apply log1p normalization during extraction (before assembly)
        if needs_log1p:
            X_batch.data = np.log1p(X_batch.data)

        X_parts.append(X_batch)

        if batch_end % 5000 == 0 or batch_end == n_meso:
            print(f"  {batch_end:,} / {n_meso:,}")

    # Stack
    print("Stacking...")
    X_combined = sparse.vstack(X_parts)
    del X_parts
    gc.collect()

    # Create AnnData (data is already normalized)
    print("Creating AnnData...")
    adata_meso = ad.AnnData(
        X=X_combined,
        obs=obs_meso.reset_index(drop=True),
        var=var_df
    )

    # Verify final data scale
    if sparse.issparse(adata_meso.X):
        final_max = adata_meso.X.data.max() if len(adata_meso.X.data) > 0 else 0
    else:
        final_max = adata_meso.X.max()
    print(f"Final data max value: {final_max:.2f} (log1p scale)")

    # Save
    output_file = PROCESSED_DIR / output_name
    print(f"Saving to {output_file}...")
    adata_meso.write(output_file)
    print(f"Saved: {output_file.stat().st_size / 1e6:.1f} MB")

    del adata
    gc.collect()

    return adata_meso


def main():
    print("\n#" + "="*58 + "#")
    print("# Epicardial Cell Extraction (Simple Version)")
    print("#" + "="*58 + "#")

    linna_dir = RAW_DIR / "linna_kuosmanen"

    # Process PERIHEART
    periheart = linna_dir / "77df1b04-4905-4dac-b9aa-6243f10201ae.h5ad"
    adata1 = extract_mesothelial_simple(periheart, "epicardial_periheart.h5ad")

    gc.collect()

    # Process CAREBANK
    carebank = linna_dir / "8e1afde5-6a0a-4ed7-ba0b-1bcdc335ddfe.h5ad"
    adata2 = extract_mesothelial_simple(carebank, "epicardial_carebank.h5ad")

    # Summary
    print("\n" + "="*60)
    print("SUMMARY")
    print("="*60)

    if adata1 is not None:
        print(f"PERIHEART: {adata1.n_obs:,} mesothelial cells")
    if adata2 is not None:
        print(f"CAREBANK: {adata2.n_obs:,} mesothelial cells")

    total = (adata1.n_obs if adata1 else 0) + (adata2.n_obs if adata2 else 0)
    print(f"\nTotal: {total:,} epicardial cells")
    print("\nOutput files:")
    print("  - data/processed/epicardial_periheart.h5ad")
    print("  - data/processed/epicardial_carebank.h5ad")


if __name__ == "__main__":
    main()
