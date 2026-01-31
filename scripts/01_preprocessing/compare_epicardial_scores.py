#!/usr/bin/env python3
"""
Compare Epicardial Scores Across All Three Datasets

Uses the same 4 markers (WT1, TBX18, ALDH1A2, UPK3B) to calculate
epicardial scores and compare distributions across:
- PERIHEART (Linna-Kuosmanen)
- CAREBANK (Linna-Kuosmanen)
- Kuppe MI
"""

import anndata as ad
import numpy as np
import pandas as pd
from pathlib import Path
from scipy import sparse

PROJECT_DIR = Path(__file__).parent.parent.parent
PROCESSED_DIR = PROJECT_DIR / "data/processed"

# Same markers used in Kuppe extraction
MARKERS = {
    'WT1': 'ENSG00000184937',
    'TBX18': 'ENSG00000112837',
    'ALDH1A2': 'ENSG00000128918',
    'UPK3B': 'ENSG00000243566',
}


def calculate_epicardial_score(adata, dataset_name):
    """Calculate epicardial score using the 4 canonical markers."""
    print(f"\n{'='*60}")
    print(f"Dataset: {dataset_name}")
    print(f"{'='*60}")
    print(f"Cells: {adata.n_obs:,}")
    print(f"Genes: {adata.n_vars:,}")

    # Find marker indices
    print("\nMarker availability:")
    marker_indices = []
    found_markers = []
    for name, ensembl in MARKERS.items():
        if ensembl in adata.var_names:
            idx = list(adata.var_names).index(ensembl)
            marker_indices.append(idx)
            found_markers.append(name)
            print(f"  + {name}: found (idx={idx})")
        else:
            print(f"  - {name}: NOT FOUND")

    if len(marker_indices) == 0:
        print("  ERROR: No markers found!")
        return None

    print(f"\nUsing {len(marker_indices)}/{len(MARKERS)} markers: {', '.join(found_markers)}")

    # Extract marker expression
    X_markers = adata.X[:, marker_indices]
    if sparse.issparse(X_markers):
        X_markers = X_markers.toarray()

    # Calculate mean expression (epicardial score)
    scores = X_markers.mean(axis=1)

    # Statistics
    print(f"\nScore distribution:")
    print(f"  Min:    {scores.min():.4f}")
    print(f"  25th:   {np.percentile(scores, 25):.4f}")
    print(f"  Median: {np.median(scores):.4f}")
    print(f"  Mean:   {scores.mean():.4f}")
    print(f"  75th:   {np.percentile(scores, 75):.4f}")
    print(f"  95th:   {np.percentile(scores, 95):.4f}")
    print(f"  Max:    {scores.max():.4f}")
    print(f"  Std:    {scores.std():.4f}")

    # Individual marker expression
    print(f"\nIndividual marker mean expression:")
    for i, name in enumerate(found_markers):
        marker_expr = X_markers[:, i]
        print(f"  {name}: mean={marker_expr.mean():.4f}, max={marker_expr.max():.4f}, "
              f"nonzero={100*(marker_expr > 0).sum()/len(marker_expr):.1f}%")

    return scores


def main():
    print("\n#" + "="*58 + "#")
    print("# Epicardial Score Comparison Across Datasets")
    print("#" + "="*58 + "#")

    results = {}

    # Load and score each dataset
    datasets = [
        ("PERIHEART", PROCESSED_DIR / "epicardial_periheart.h5ad"),
        ("CAREBANK", PROCESSED_DIR / "epicardial_carebank.h5ad"),
        ("Kuppe_MI", PROCESSED_DIR / "epicardial_kuppe.h5ad"),
    ]

    for name, filepath in datasets:
        if not filepath.exists():
            print(f"\nSkipping {name}: file not found")
            continue

        adata = ad.read_h5ad(filepath)
        scores = calculate_epicardial_score(adata, name)

        if scores is not None:
            results[name] = {
                'n_cells': len(scores),
                'scores': scores,
                'mean': scores.mean(),
                'median': np.median(scores),
                'std': scores.std(),
                'p95': np.percentile(scores, 95),
            }

        del adata

    # Comparison summary
    print("\n" + "="*60)
    print("COMPARISON SUMMARY")
    print("="*60)

    print("\n{:<15} {:>10} {:>10} {:>10} {:>10} {:>10}".format(
        "Dataset", "N cells", "Mean", "Median", "Std", "95th %"
    ))
    print("-"*60)

    for name, data in results.items():
        print("{:<15} {:>10,} {:>10.4f} {:>10.4f} {:>10.4f} {:>10.4f}".format(
            name, data['n_cells'], data['mean'], data['median'],
            data['std'], data['p95']
        ))

    # Check if scores are comparable
    print("\n" + "="*60)
    print("INTERPRETATION")
    print("="*60)

    if len(results) >= 2:
        means = [data['mean'] for data in results.values()]
        max_mean = max(means)
        min_mean = min(means)

        if max_mean > 0:
            ratio = max_mean / min_mean if min_mean > 0 else float('inf')
            print(f"\nMean score ratio (max/min): {ratio:.2f}x")

            if ratio < 2:
                print("  -> Scores are comparable across datasets")
            elif ratio < 5:
                print("  -> Moderate difference in score levels")
                print("     May need normalization before combining")
            else:
                print("  -> Large difference in score levels")
                print("     Recommend careful normalization")

        # Check if Kuppe has pre-computed scores
        if "Kuppe_MI" in results:
            kuppe_file = PROCESSED_DIR / "epicardial_kuppe.h5ad"
            adata_kuppe = ad.read_h5ad(kuppe_file)
            if 'epicardial_score' in adata_kuppe.obs.columns:
                stored_scores = adata_kuppe.obs['epicardial_score'].values
                print(f"\nKuppe stored epicardial_score vs recalculated:")
                print(f"  Stored mean:  {stored_scores.mean():.4f}")
                print(f"  Recalc mean:  {results['Kuppe_MI']['mean']:.4f}")
                print(f"  Correlation:  {np.corrcoef(stored_scores, results['Kuppe_MI']['scores'])[0,1]:.4f}")

    print("\n" + "="*60)
    print("Done!")
    print("="*60 + "\n")


if __name__ == "__main__":
    main()
