#!/usr/bin/env python3
"""
Prepare Communication Datasets for LIANA/NicheNet Analysis

Creates two datasets:
1. communication_kuppe.h5ad - Kuppe epicardial + Kuppe sender cells
2. communication_merged.h5ad - All epicardial + all sender cells (batch-corrected)

Sender cell types:
- Cardiomyocytes
- Fibroblasts
- Endothelial cells
- Macrophages/Myeloid cells

Usage:
    python prepare_communication_data.py              # Full data
    python prepare_communication_data.py --sample 0.01  # 1% sample for development
    python prepare_communication_data.py --force      # Force regenerate all
    python prepare_communication_data.py --only merged  # Only generate merged dataset
"""

import argparse
import anndata as ad
import scanpy as sc
import numpy as np
import pandas as pd
from pathlib import Path
from scipy import sparse

PROJECT_DIR = Path(__file__).parent.parent.parent
RAW_DIR = PROJECT_DIR / "data/raw"
PROCESSED_DIR = PROJECT_DIR / "data/processed"

# Global settings (set via command line)
SAMPLE_FRACTION = None
FORCE_REGENERATE = False

# Cell type mapping for sender cells
KUPPE_SENDER_TYPES = {
    'Cardiomyocyte': ['cardiac muscle myoblast'],
    'Fibroblast': ['fibroblast of cardiac tissue'],
    'Endothelial': ['cardiac endothelial cell'],
    'Myeloid': [],  # Use cell_type_original for this
}

KUPPE_SENDER_ORIGINAL = {
    'Myeloid': ['Myeloid'],
}

LINNA_SENDER_TYPES = {
    'Cardiomyocyte': ['cardiac muscle cell'],
    'Fibroblast': ['fibroblast'],
    'Endothelial': ['endocardial cell', 'cardiac blood vessel endothelial cell'],
    'Macrophage': ['macrophage'],
}


def check_output_exists(name):
    """Check if output file already exists."""
    suffix = "_sample" if SAMPLE_FRACTION else ""
    output_path = PROCESSED_DIR / f"{name}{suffix}.h5ad"
    return output_path.exists()


def sample_adata(adata, fraction, random_state=42):
    """Randomly sample a fraction of cells from AnnData object."""
    if fraction is None or fraction >= 1.0:
        return adata
    n_sample = max(int(adata.n_obs * fraction), 100)  # At least 100 cells
    n_sample = min(n_sample, adata.n_obs)
    np.random.seed(random_state)
    indices = np.random.choice(adata.n_obs, n_sample, replace=False)
    return adata[indices].copy()


def extract_kuppe_sender_cells():
    """Extract sender cells from Kuppe dataset."""
    print("\n" + "="*60)
    print("Extracting Kuppe sender cells")
    print("="*60)

    adata = ad.read_h5ad(RAW_DIR / "kuppe/54d24dbe-d39a-4844-bb21-07b5f4e173ad.h5ad")
    print(f"Loaded: {adata.n_obs:,} cells")

    # Build mask for sender cells
    mask = np.zeros(adata.n_obs, dtype=bool)

    # By cell_type
    for sender_type, cell_types in KUPPE_SENDER_TYPES.items():
        for ct in cell_types:
            ct_mask = adata.obs['cell_type'] == ct
            mask |= ct_mask
            print(f"  {sender_type} ({ct}): {ct_mask.sum():,} cells")

    # By cell_type_original (for Myeloid)
    for sender_type, cell_types in KUPPE_SENDER_ORIGINAL.items():
        for ct in cell_types:
            ct_mask = adata.obs['cell_type_original'] == ct
            mask |= ct_mask
            print(f"  {sender_type} (original: {ct}): {ct_mask.sum():,} cells")

    # Subset
    sender = adata[mask].copy()
    print(f"\nTotal sender cells: {sender.n_obs:,}")

    # Add standardized cell type label
    sender.obs['sender_type'] = 'unknown'
    for ct in KUPPE_SENDER_TYPES['Cardiomyocyte']:
        sender.obs.loc[sender.obs['cell_type'] == ct, 'sender_type'] = 'Cardiomyocyte'
    for ct in KUPPE_SENDER_TYPES['Fibroblast']:
        sender.obs.loc[sender.obs['cell_type'] == ct, 'sender_type'] = 'Fibroblast'
    for ct in KUPPE_SENDER_TYPES['Endothelial']:
        sender.obs.loc[sender.obs['cell_type'] == ct, 'sender_type'] = 'Endothelial'
    sender.obs.loc[sender.obs['cell_type_original'] == 'Myeloid', 'sender_type'] = 'Myeloid'

    sender.obs['dataset'] = 'Kuppe'
    sender.obs['cell_role'] = 'sender'

    # Apply sampling if requested
    if SAMPLE_FRACTION is not None:
        sender = sample_adata(sender, SAMPLE_FRACTION)
        print(f"\n  Sampled to {sender.n_obs:,} cells ({SAMPLE_FRACTION*100:.1f}%)")

    print("\nSender type distribution:")
    print(sender.obs['sender_type'].value_counts())

    del adata
    return sender


def extract_linna_sender_cells(filepath, dataset_name):
    """Extract sender cells from a Linna-Kuosmanen dataset."""
    print(f"\n{'='*60}")
    print(f"Extracting {dataset_name} sender cells")
    print("="*60)

    adata = ad.read_h5ad(filepath)
    print(f"Loaded: {adata.n_obs:,} cells")

    # Apply log1p if needed (Linna-Kuosmanen is in CPM format)
    sample = adata.X[:100, :100]
    if sparse.issparse(sample):
        sample = sample.toarray()
    if sample.max() > 50:  # Likely not log-transformed
        print("  Applying log1p transformation...")
        if sparse.issparse(adata.X):
            adata.X = adata.X.log1p()
        else:
            adata.X = np.log1p(adata.X)

    # Build mask for sender cells
    mask = np.zeros(adata.n_obs, dtype=bool)

    for sender_type, cell_types in LINNA_SENDER_TYPES.items():
        for ct in cell_types:
            ct_mask = adata.obs['cell_type'] == ct
            mask |= ct_mask
            print(f"  {sender_type} ({ct}): {ct_mask.sum():,} cells")

    # Subset
    sender = adata[mask].copy()
    print(f"\nTotal sender cells: {sender.n_obs:,}")

    # Add standardized cell type label
    sender.obs['sender_type'] = 'unknown'
    for ct in LINNA_SENDER_TYPES['Cardiomyocyte']:
        sender.obs.loc[sender.obs['cell_type'] == ct, 'sender_type'] = 'Cardiomyocyte'
    for ct in LINNA_SENDER_TYPES['Fibroblast']:
        sender.obs.loc[sender.obs['cell_type'] == ct, 'sender_type'] = 'Fibroblast'
    for ct in LINNA_SENDER_TYPES['Endothelial']:
        sender.obs.loc[sender.obs['cell_type'] == ct, 'sender_type'] = 'Endothelial'
    for ct in LINNA_SENDER_TYPES['Macrophage']:
        sender.obs.loc[sender.obs['cell_type'] == ct, 'sender_type'] = 'Macrophage'

    sender.obs['dataset'] = dataset_name
    sender.obs['cell_role'] = 'sender'

    # Apply sampling if requested
    if SAMPLE_FRACTION is not None:
        sender = sample_adata(sender, SAMPLE_FRACTION)
        print(f"\n  Sampled to {sender.n_obs:,} cells ({SAMPLE_FRACTION*100:.1f}%)")

    print("\nSender type distribution:")
    print(sender.obs['sender_type'].value_counts())

    del adata
    return sender


def create_kuppe_communication_dataset():
    """Create Kuppe-only communication dataset."""
    suffix = "_sample" if SAMPLE_FRACTION else ""
    output_path = PROCESSED_DIR / f"communication_kuppe{suffix}.h5ad"

    # Check if already exists
    if output_path.exists() and not FORCE_REGENERATE:
        print("\n" + "#"*60)
        print(f"# SKIP: communication_kuppe{suffix}.h5ad already exists")
        print("#"*60)
        return

    print("\n" + "#"*60)
    print(f"# Creating communication_kuppe{suffix}.h5ad")
    print("#"*60)

    # Load Kuppe epicardial cells
    print("\nLoading Kuppe epicardial cells...")
    epi = ad.read_h5ad(PROCESSED_DIR / "epicardial_kuppe.h5ad")
    epi.obs['sender_type'] = 'Epicardial'
    epi.obs['cell_role'] = 'receiver'
    epi.obs['dataset'] = 'Kuppe'

    # Apply sampling if requested
    if SAMPLE_FRACTION is not None:
        epi = sample_adata(epi, SAMPLE_FRACTION)
        print(f"  Epicardial: {epi.n_obs:,} cells (sampled {SAMPLE_FRACTION*100:.1f}%)")
    else:
        print(f"  Epicardial: {epi.n_obs:,} cells")

    # Extract Kuppe sender cells
    sender = extract_kuppe_sender_cells()

    # Find common genes
    common_genes = list(set(epi.var_names) & set(sender.var_names))
    print(f"\nCommon genes: {len(common_genes):,}")

    # Subset to common genes
    epi = epi[:, common_genes].copy()
    sender = sender[:, common_genes].copy()

    # Concatenate
    print("\nConcatenating...")
    combined = ad.concat([epi, sender], join='outer')
    combined.obs_names_make_unique()

    print(f"\nCombined dataset: {combined.n_obs:,} cells, {combined.n_vars:,} genes")
    print("\nCell role distribution:")
    print(combined.obs['cell_role'].value_counts())
    print("\nSender type distribution:")
    print(combined.obs['sender_type'].value_counts())

    # Save
    print(f"\nSaving to {output_path}...")
    combined.write_h5ad(output_path)

    del epi, sender, combined
    print("Done!")


def create_merged_communication_dataset():
    """Create merged communication dataset with all cells."""
    suffix = "_sample" if SAMPLE_FRACTION else ""
    output_path = PROCESSED_DIR / f"communication_merged{suffix}.h5ad"

    # Check if already exists
    if output_path.exists() and not FORCE_REGENERATE:
        print("\n" + "#"*60)
        print(f"# SKIP: communication_merged{suffix}.h5ad already exists")
        print("#"*60)
        return

    print("\n" + "#"*60)
    print(f"# Creating communication_merged{suffix}.h5ad")
    print("#"*60)

    datasets = []

    # Load all epicardial cells
    print("\nLoading epicardial cells...")
    epi = ad.read_h5ad(PROCESSED_DIR / "epicardial_with_states.h5ad")
    epi.obs['sender_type'] = 'Epicardial'
    epi.obs['cell_role'] = 'receiver'

    # Apply sampling if requested
    if SAMPLE_FRACTION is not None:
        epi = sample_adata(epi, SAMPLE_FRACTION)
        print(f"  Epicardial: {epi.n_obs:,} cells (sampled {SAMPLE_FRACTION*100:.1f}%)")
    else:
        print(f"  Epicardial: {epi.n_obs:,} cells")
    datasets.append(epi)

    # Extract Kuppe sender cells
    kuppe_sender = extract_kuppe_sender_cells()
    datasets.append(kuppe_sender)

    # Extract Linna-Kuosmanen sender cells
    periheart_sender = extract_linna_sender_cells(
        RAW_DIR / "linna_kuosmanen/77df1b04-4905-4dac-b9aa-6243f10201ae.h5ad",
        "PERIHEART"
    )
    datasets.append(periheart_sender)

    carebank_sender = extract_linna_sender_cells(
        RAW_DIR / "linna_kuosmanen/8e1afde5-6a0a-4ed7-ba0b-1bcdc335ddfe.h5ad",
        "CAREBANK"
    )
    datasets.append(carebank_sender)

    # Find common genes
    print("\n" + "="*60)
    print("Finding common genes...")
    common_genes = set(datasets[0].var_names)
    for adata in datasets[1:]:
        common_genes = common_genes & set(adata.var_names)
    common_genes = sorted(list(common_genes))
    print(f"Common genes: {len(common_genes):,}")

    # Subset to common genes
    for i in range(len(datasets)):
        datasets[i] = datasets[i][:, common_genes].copy()

    # Concatenate
    print("\nConcatenating all datasets...")
    combined = ad.concat(datasets, join='outer')
    combined.obs_names_make_unique()

    print(f"\nCombined dataset: {combined.n_obs:,} cells, {combined.n_vars:,} genes")
    print("\nCell role distribution:")
    print(combined.obs['cell_role'].value_counts())
    print("\nDataset distribution:")
    print(combined.obs['dataset'].value_counts())
    print("\nSender type distribution:")
    print(combined.obs['sender_type'].value_counts())

    # Handle any infinity or NaN values (optimized for sparse matrices)
    print("\n" + "="*60)
    print("Cleaning data...")
    if sparse.issparse(combined.X):
        # Only process non-zero elements (much faster for sparse data)
        combined.X.data = np.nan_to_num(combined.X.data, nan=0.0, posinf=0.0, neginf=0.0)
        # Clip values to prevent expm1 overflow (max safe value ~88 for float32)
        combined.X.data = np.clip(combined.X.data, 0, 50)
        combined.X.eliminate_zeros()
    else:
        combined.X = np.nan_to_num(combined.X, nan=0.0, posinf=0.0, neginf=0.0)
        combined.X = np.clip(combined.X, 0, 50)
    print("  Removed NaN/Inf values and clipped extreme values")

    # Batch correction with Harmony
    print("\n" + "="*60)
    print("Batch correction with Harmony...")

    # Need to compute PCA first
    try:
        # Check if harmony dependencies are available
        import scanpy.external as sce

        # Use seurat_v3 flavor to avoid expm1 overflow issues
        sc.pp.highly_variable_genes(combined, n_top_genes=2000, batch_key='dataset', flavor='seurat_v3')
        sc.tl.pca(combined, n_comps=50)

        # Harmony integration
        sce.pp.harmony_integrate(combined, key='dataset')
        print("Harmony integration complete!")
    except ImportError as e:
        print(f"WARNING: Harmony dependency missing: {e}")
        print("Install with: pip install scikit-misc harmonypy")
        print("Proceeding without batch correction...")
    except Exception as e:
        print(f"WARNING: Batch correction failed: {e}")
        print("Proceeding without batch correction...")

    # Save
    print(f"\nSaving to {output_path}...")
    combined.write_h5ad(output_path)

    print("\n" + "="*60)
    print("Done!")
    print("="*60)


def main():
    global SAMPLE_FRACTION, FORCE_REGENERATE

    parser = argparse.ArgumentParser(description='Prepare communication datasets for LIANA/NicheNet')
    parser.add_argument('--sample', type=float, default=None,
                        help='Sample fraction (e.g., 0.01 for 1%%). Default: use all data')
    parser.add_argument('--force', action='store_true',
                        help='Force regenerate even if output files exist')
    parser.add_argument('--only', choices=['kuppe', 'merged'],
                        help='Only generate specified dataset')
    args = parser.parse_args()

    SAMPLE_FRACTION = args.sample
    FORCE_REGENERATE = args.force

    print("\n" + "#"*60)
    print("# Prepare Communication Datasets")
    if SAMPLE_FRACTION:
        print(f"# MODE: Sampling {SAMPLE_FRACTION*100:.1f}% of data")
    else:
        print("# MODE: Full data")
    if FORCE_REGENERATE:
        print("# FORCE: Will regenerate existing files")
    print("#"*60)

    # Create datasets based on --only flag
    if args.only == 'kuppe':
        create_kuppe_communication_dataset()
    elif args.only == 'merged':
        create_merged_communication_dataset()
    else:
        # Create both
        create_kuppe_communication_dataset()
        create_merged_communication_dataset()


if __name__ == "__main__":
    main()
