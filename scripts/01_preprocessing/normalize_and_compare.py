#!/usr/bin/env python3
"""
Normalize Linna-Kuosmanen Data and Compare Epicardial Scores

This script:
1. Applies log1p transformation to Linna-Kuosmanen datasets (PERIHEART, CAREBANK)
2. Saves the normalized versions
3. Compares epicardial scores across all three datasets using consistent normalization
"""

import anndata as ad
import numpy as np
from scipy import sparse
from pathlib import Path

PROJECT_DIR = Path(__file__).parent.parent.parent
PROCESSED_DIR = PROJECT_DIR / "data/processed"

# Epicardial markers
MARKERS = {
    'WT1': 'ENSG00000184937',
    'TBX18': 'ENSG00000112837',
    'ALDH1A2': 'ENSG00000128918',
    'UPK3B': 'ENSG00000243566',
}


def apply_log1p_and_save(input_path, output_path, dataset_name):
    """Apply log1p transformation and save."""
    print(f"\n{'='*60}")
    print(f"Processing: {dataset_name}")
    print(f"{'='*60}")

    print(f"Loading {input_path}...")
    adata = ad.read_h5ad(input_path)
    print(f"  Shape: {adata.shape}")

    # Check current range
    if sparse.issparse(adata.X):
        sample = adata.X[:100, :].toarray()
    else:
        sample = adata.X[:100, :]
    print(f"  Before log1p - range: {sample.min():.3f} - {sample.max():.3f}")

    # Apply log1p
    print("  Applying log1p transformation...")
    if sparse.issparse(adata.X):
        adata.X = adata.X.log1p()
    else:
        adata.X = np.log1p(adata.X)

    # Check new range
    if sparse.issparse(adata.X):
        sample = adata.X[:100, :].toarray()
    else:
        sample = adata.X[:100, :]
    print(f"  After log1p - range: {sample.min():.3f} - {sample.max():.3f}")

    # Add normalization info to uns
    adata.uns['normalization'] = 'log1p'
    adata.uns['X_approximate_distribution'] = 'normal'

    # Save
    print(f"  Saving to {output_path}...")
    adata.write_h5ad(output_path)
    print(f"  Done!")

    return adata


def calculate_epicardial_score(adata, dataset_name):
    """Calculate epicardial score using the 4 canonical markers."""
    # Find marker indices
    marker_indices = []
    found_markers = []
    for name, ensembl in MARKERS.items():
        if ensembl in adata.var_names:
            idx = list(adata.var_names).index(ensembl)
            marker_indices.append(idx)
            found_markers.append(name)

    if len(marker_indices) == 0:
        return None, None

    # Extract marker expression
    X_markers = adata.X[:, marker_indices]
    if sparse.issparse(X_markers):
        X_markers = X_markers.toarray()

    # Calculate mean expression (epicardial score)
    scores = X_markers.mean(axis=1)

    # Individual marker stats
    marker_stats = {}
    for i, name in enumerate(found_markers):
        marker_expr = X_markers[:, i]
        marker_stats[name] = {
            'mean': float(marker_expr.mean()),
            'max': float(marker_expr.max()),
            'pct_nonzero': float(100 * (marker_expr > 0).sum() / len(marker_expr))
        }

    return scores, marker_stats


def main():
    print("\n" + "#"*60)
    print("# Normalize Linna-Kuosmanen and Compare Epicardial Scores")
    print("#"*60)

    results = {}

    # Step 1: Apply log1p to PERIHEART
    periheart_in = PROCESSED_DIR / "epicardial_periheart.h5ad"
    periheart_out = PROCESSED_DIR / "epicardial_periheart_log.h5ad"
    adata_periheart = apply_log1p_and_save(periheart_in, periheart_out, "PERIHEART")

    scores, marker_stats = calculate_epicardial_score(adata_periheart, "PERIHEART")
    results["PERIHEART"] = {
        'n_cells': len(scores),
        'scores': scores,
        'marker_stats': marker_stats,
    }
    del adata_periheart

    # Step 2: Apply log1p to CAREBANK
    carebank_in = PROCESSED_DIR / "epicardial_carebank.h5ad"
    carebank_out = PROCESSED_DIR / "epicardial_carebank_log.h5ad"
    adata_carebank = apply_log1p_and_save(carebank_in, carebank_out, "CAREBANK")

    scores, marker_stats = calculate_epicardial_score(adata_carebank, "CAREBANK")
    results["CAREBANK"] = {
        'n_cells': len(scores),
        'scores': scores,
        'marker_stats': marker_stats,
    }
    del adata_carebank

    # Step 3: Load Kuppe (already log-normalized)
    print(f"\n{'='*60}")
    print("Loading: Kuppe_MI (already log-normalized)")
    print(f"{'='*60}")
    kuppe_path = PROCESSED_DIR / "epicardial_kuppe.h5ad"
    adata_kuppe = ad.read_h5ad(kuppe_path)
    print(f"  Shape: {adata_kuppe.shape}")

    scores, marker_stats = calculate_epicardial_score(adata_kuppe, "Kuppe_MI")
    results["Kuppe_MI"] = {
        'n_cells': len(scores),
        'scores': scores,
        'marker_stats': marker_stats,
    }
    del adata_kuppe

    # Step 4: Comparison
    print("\n" + "="*60)
    print("COMPARISON AFTER LOG1P NORMALIZATION")
    print("="*60)

    print("\n{:<12} {:>10} {:>10} {:>10} {:>10} {:>10}".format(
        "Dataset", "N cells", "Mean", "Median", "Std", "95th %"
    ))
    print("-"*60)

    for name, data in results.items():
        scores = data['scores']
        print("{:<12} {:>10,} {:>10.4f} {:>10.4f} {:>10.4f} {:>10.4f}".format(
            name,
            data['n_cells'],
            float(np.mean(scores)),
            float(np.median(scores)),
            float(np.std(scores)),
            float(np.percentile(scores, 95))
        ))

    # Individual marker comparison
    print("\n" + "="*60)
    print("INDIVIDUAL MARKER EXPRESSION (Mean)")
    print("="*60)
    print("\n{:<12} {:>10} {:>10} {:>10} {:>10}".format(
        "Dataset", "WT1", "TBX18", "ALDH1A2", "UPK3B"
    ))
    print("-"*60)

    for name, data in results.items():
        ms = data['marker_stats']
        print("{:<12} {:>10.4f} {:>10.4f} {:>10.4f} {:>10.4f}".format(
            name,
            ms.get('WT1', {}).get('mean', 0),
            ms.get('TBX18', {}).get('mean', 0),
            ms.get('ALDH1A2', {}).get('mean', 0),
            ms.get('UPK3B', {}).get('mean', 0)
        ))

    # Non-zero percentage
    print("\n{:<12} {:>10} {:>10} {:>10} {:>10}".format(
        "Dataset", "WT1 %", "TBX18 %", "ALDH1A2 %", "UPK3B %"
    ))
    print("-"*60)

    for name, data in results.items():
        ms = data['marker_stats']
        print("{:<12} {:>10.1f} {:>10.1f} {:>10.1f} {:>10.1f}".format(
            name,
            ms.get('WT1', {}).get('pct_nonzero', 0),
            ms.get('TBX18', {}).get('pct_nonzero', 0),
            ms.get('ALDH1A2', {}).get('pct_nonzero', 0),
            ms.get('UPK3B', {}).get('pct_nonzero', 0)
        ))

    # Score ratio check
    print("\n" + "="*60)
    print("INTERPRETATION")
    print("="*60)

    means = {name: float(np.mean(data['scores'])) for name, data in results.items()}
    max_mean = max(means.values())
    min_mean = min(means.values())
    ratio = max_mean / min_mean if min_mean > 0 else float('inf')

    print(f"\nScore range: {min_mean:.4f} - {max_mean:.4f}")
    print(f"Max/Min ratio: {ratio:.2f}x")

    if ratio < 2:
        print("\n✓ After log1p normalization, scores are now COMPARABLE!")
        print("  The datasets can be safely combined for downstream analysis.")
    else:
        print(f"\n⚠ Still significant difference (ratio={ratio:.2f}x)")
        print("  May need additional batch correction.")

    print("\n" + "="*60)
    print("Output files:")
    print(f"  - {PROCESSED_DIR}/epicardial_periheart_log.h5ad")
    print(f"  - {PROCESSED_DIR}/epicardial_carebank_log.h5ad")
    print(f"  - {PROCESSED_DIR}/epicardial_kuppe.h5ad (unchanged)")
    print("="*60 + "\n")


if __name__ == "__main__":
    main()
