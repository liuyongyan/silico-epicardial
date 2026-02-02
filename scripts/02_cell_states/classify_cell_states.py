#!/usr/bin/env python3
"""
Phase 2: Classify Epicardial Cell States

Two-method approach with cross-validation:
1. Condition-based: Using disease/spatial zone annotations
2. Score-based: Using EMT and proliferation molecular signatures

Improved methodology:
- sc.tl.score_genes() with control gene sets (background correction)
- GMM-based threshold selection instead of arbitrary percentile

Output: Consensus classification (quiescent vs activated)
"""

import anndata as ad
import scanpy as sc
import numpy as np
import pandas as pd
from pathlib import Path
from scipy import sparse
from sklearn.mixture import GaussianMixture
import warnings
warnings.filterwarnings('ignore')

PROJECT_DIR = Path(__file__).parent.parent.parent
PROCESSED_DIR = PROJECT_DIR / "data/processed"

# =============================================================================
# Gene signatures (Ensembl IDs)
# =============================================================================

# Proliferation markers
PROLIFERATION_GENES = {
    'MKI67': 'ENSG00000148773',
    'TOP2A': 'ENSG00000131747',
    'PCNA': 'ENSG00000132646',
    'CDK1': 'ENSG00000170461',
    'CCNB1': 'ENSG00000134057',
    'CCNB2': 'ENSG00000157456',
    'CCNA2': 'ENSG00000145386',
    'MCM2': 'ENSG00000073111',
    'MCM6': 'ENSG00000076003',
    'AURKA': 'ENSG00000087586',
}

# EMT upregulated markers
EMT_UP_GENES = {
    'SNAI1': 'ENSG00000124216',
    'SNAI2': 'ENSG00000019549',
    'TWIST1': 'ENSG00000122691',
    'TWIST2': 'ENSG00000233608',
    'ZEB1': 'ENSG00000148516',
    'ZEB2': 'ENSG00000169554',
    'VIM': 'ENSG00000026025',
    'CDH2': 'ENSG00000170558',  # N-cadherin
    'FN1': 'ENSG00000115414',   # Fibronectin
    'ACTA2': 'ENSG00000107796', # alpha-SMA
    'POSTN': 'ENSG00000133110', # Periostin
    'COL1A1': 'ENSG00000108821',
    'COL3A1': 'ENSG00000168542',
}

# EMT downregulated markers (epithelial)
EMT_DOWN_GENES = {
    'CDH1': 'ENSG00000039068',  # E-cadherin
    'CLDN1': 'ENSG00000163347',
    'CLDN3': 'ENSG00000165215',
    'TJP1': 'ENSG00000104067',  # ZO-1
    'OCLN': 'ENSG00000197822',
}


def get_available_genes(adata, gene_dict):
    """Filter gene dictionary to only include genes present in dataset."""
    available = {}
    for name, ensembl in gene_dict.items():
        if ensembl in adata.var_names:
            available[name] = ensembl
    return available


def find_threshold_gmm(scores, n_components=2):
    """
    Find threshold using Gaussian Mixture Model.

    Fits a 2-component GMM and finds the intersection point between the two distributions.
    This is more principled than using an arbitrary percentile.
    """
    scores_array = scores.values.reshape(-1, 1)

    gmm = GaussianMixture(n_components=n_components, random_state=42, n_init=10)
    gmm.fit(scores_array)

    # Get component means and identify low/high
    means = gmm.means_.flatten()
    low_idx = np.argmin(means)
    high_idx = np.argmax(means)

    # Find threshold at intersection (where posterior probabilities are equal)
    # Use the midpoint between means as approximation
    threshold = (means[low_idx] + means[high_idx]) / 2

    # Get proportion in each component
    labels = gmm.predict(scores_array)
    n_low = (labels == low_idx).sum()
    n_high = (labels == high_idx).sum()

    return threshold, means, (n_low, n_high), gmm


def assign_condition_labels(adata):
    """
    Method 1: Assign cell state labels based on condition annotations.

    Quiescent: normal/healthy tissue
    Activated: MI tissue
    """
    print("\n" + "="*60)
    print("Method 1: Condition-based Classification")
    print("="*60)

    # Simple classification based on disease annotation
    # (Linna-Kuosmanen data only - PERIHEART and CAREBANK)
    disease_map = {'normal': 'quiescent', 'myocardial infarction': 'activated'}
    adata.obs['condition_state'] = adata.obs['disease'].astype(str).map(disease_map).fillna('unknown')

    # Binary classification (same as condition_state for our filtered data)
    adata.obs['condition_binary'] = adata.obs['condition_state']

    # Summary
    print("\nCondition state distribution:")
    print(adata.obs['condition_state'].value_counts())
    print("\nBinary classification:")
    print(adata.obs['condition_binary'].value_counts())


def calculate_molecular_scores(adata):
    """
    Method 2: Calculate EMT and proliferation scores using sc.tl.score_genes().

    Improved methodology:
    - Uses sc.tl.score_genes() which compares gene set expression to a control
      gene set with similar expression levels (background correction)
    - Provides more robust scoring that accounts for technical variation
    """
    print("\n" + "="*60)
    print("Method 2: Molecular Signature Scores (sc.tl.score_genes)")
    print("="*60)

    # Proliferation score
    print("\nProliferation signature:")
    available_prolif = get_available_genes(adata, PROLIFERATION_GENES)
    prolif_genes = list(available_prolif.values())
    print(f"  Available genes: {len(prolif_genes)}/{len(PROLIFERATION_GENES)}")
    print(f"  Genes: {list(available_prolif.keys())}")

    sc.tl.score_genes(adata, gene_list=prolif_genes, score_name='proliferation_score',
                      ctrl_size=50, n_bins=25, use_raw=False)
    print(f"  Score: mean={adata.obs['proliferation_score'].mean():.4f}, "
          f"std={adata.obs['proliferation_score'].std():.4f}")

    # EMT up score
    print("\nEMT upregulated signature:")
    available_emt_up = get_available_genes(adata, EMT_UP_GENES)
    emt_up_genes = list(available_emt_up.values())
    print(f"  Available genes: {len(emt_up_genes)}/{len(EMT_UP_GENES)}")
    print(f"  Genes: {list(available_emt_up.keys())}")

    sc.tl.score_genes(adata, gene_list=emt_up_genes, score_name='emt_up_score',
                      ctrl_size=50, n_bins=25, use_raw=False)
    print(f"  Score: mean={adata.obs['emt_up_score'].mean():.4f}, "
          f"std={adata.obs['emt_up_score'].std():.4f}")

    # EMT down score
    print("\nEMT downregulated (epithelial) signature:")
    available_emt_down = get_available_genes(adata, EMT_DOWN_GENES)
    emt_down_genes = list(available_emt_down.values())
    print(f"  Available genes: {len(emt_down_genes)}/{len(EMT_DOWN_GENES)}")
    print(f"  Genes: {list(available_emt_down.keys())}")

    sc.tl.score_genes(adata, gene_list=emt_down_genes, score_name='emt_down_score',
                      ctrl_size=50, n_bins=25, use_raw=False)
    print(f"  Score: mean={adata.obs['emt_down_score'].mean():.4f}, "
          f"std={adata.obs['emt_down_score'].std():.4f}")

    # Combined EMT score (up - down)
    adata.obs['emt_score'] = adata.obs['emt_up_score'] - adata.obs['emt_down_score']
    print(f"\nCombined EMT score: mean={adata.obs['emt_score'].mean():.4f}, "
          f"std={adata.obs['emt_score'].std():.4f}")

    # Combined activation score (average of z-normalized scores)
    prolif_z = (adata.obs['proliferation_score'] - adata.obs['proliferation_score'].mean()) / adata.obs['proliferation_score'].std()
    emt_z = (adata.obs['emt_score'] - adata.obs['emt_score'].mean()) / adata.obs['emt_score'].std()
    adata.obs['activation_score'] = (prolif_z + emt_z) / 2

    print(f"\nCombined activation score (z-normalized): "
          f"mean={adata.obs['activation_score'].mean():.4f}, "
          f"std={adata.obs['activation_score'].std():.4f}")


def classify_by_scores(adata):
    """
    Classify cells as activated based on molecular scores using GMM threshold.

    Uses Gaussian Mixture Model to find a data-driven threshold instead of
    arbitrary percentile cutoff.
    """
    print("\n" + "="*60)
    print("Score-based Classification (GMM threshold)")
    print("="*60)

    scores = adata.obs['activation_score']
    threshold, means, counts, gmm = find_threshold_gmm(scores)

    print(f"\nGMM analysis:")
    print(f"  Low component:  mean={means[np.argmin(means)]:.4f}, n={counts[0]:,}")
    print(f"  High component: mean={means[np.argmax(means)]:.4f}, n={counts[1]:,}")
    print(f"  Threshold (intersection): {threshold:.4f}")

    # For comparison, show what percentile this corresponds to
    percentile = (scores < threshold).mean() * 100
    print(f"  Equivalent percentile: {percentile:.1f}%")

    adata.obs['score_binary'] = adata.obs['activation_score'].apply(
        lambda x: 'activated' if x > threshold else 'quiescent'
    )

    print("\nScore-based classification:")
    print(adata.obs['score_binary'].value_counts())

    return threshold


def cross_validate(adata):
    """
    Cross-validate the two methods and create consensus labels.
    """
    print("\n" + "="*60)
    print("Cross-Validation")
    print("="*60)

    # Filter to cells with valid classifications from both methods
    valid_mask = (adata.obs['condition_binary'].isin(['quiescent', 'activated'])) & \
                 (adata.obs['score_binary'].isin(['quiescent', 'activated']))

    valid_obs = adata.obs[valid_mask]

    print(f"\nCells with valid labels from both methods: {valid_mask.sum():,}")

    # Concordance matrix
    concordance = pd.crosstab(
        valid_obs['condition_binary'],
        valid_obs['score_binary'],
        margins=True
    )
    print("\nConcordance matrix:")
    print(concordance)

    # Calculate agreement
    agree = ((valid_obs['condition_binary'] == 'quiescent') &
             (valid_obs['score_binary'] == 'quiescent')).sum() + \
            ((valid_obs['condition_binary'] == 'activated') &
             (valid_obs['score_binary'] == 'activated')).sum()

    agreement_rate = agree / len(valid_obs) * 100
    print(f"\nAgreement rate: {agreement_rate:.1f}%")

    # Compare scores between condition groups
    print("\n" + "-"*60)
    print("Score comparison by condition:")
    for state in ['quiescent', 'activated']:
        mask = adata.obs['condition_binary'] == state
        if mask.sum() > 0:
            prolif = adata.obs.loc[mask, 'proliferation_score'].mean()
            emt = adata.obs.loc[mask, 'emt_score'].mean()
            act = adata.obs.loc[mask, 'activation_score'].mean()
            print(f"  {state:12}: prolif={prolif:.4f}, emt={emt:.4f}, activation={act:.4f} (n={mask.sum():,})")

    # Create consensus labels using Option A naming:
    # - activated: MI + high score (true responders)
    # - bystander: MI + low score (in injury zone but not responding)
    # - primed: Normal + high score (constitutively active)
    # - quiescent: Normal + low score (resting)
    print("\n" + "-"*60)
    print("Creating cell state labels (Option A naming)...")

    adata.obs['cell_state'] = 'unknown'

    # MI + high score = activated (true responders)
    activated_mask = (adata.obs['condition_binary'] == 'activated') & \
                     (adata.obs['score_binary'] == 'activated')
    adata.obs.loc[activated_mask, 'cell_state'] = 'activated'

    # MI + low score = bystander (in injury zone but not responding)
    bystander_mask = (adata.obs['condition_binary'] == 'activated') & \
                     (adata.obs['score_binary'] == 'quiescent')
    adata.obs.loc[bystander_mask, 'cell_state'] = 'bystander'

    # Normal + high score = primed (constitutively active)
    primed_mask = (adata.obs['condition_binary'] == 'quiescent') & \
                  (adata.obs['score_binary'] == 'activated')
    adata.obs.loc[primed_mask, 'cell_state'] = 'primed'

    # Normal + low score = quiescent (resting)
    quiescent_mask = (adata.obs['condition_binary'] == 'quiescent') & \
                     (adata.obs['score_binary'] == 'quiescent')
    adata.obs.loc[quiescent_mask, 'cell_state'] = 'quiescent'

    print("\nCell state distribution:")
    print(adata.obs['cell_state'].value_counts())

    # Also keep legacy columns for compatibility
    adata.obs['consensus_state'] = adata.obs['cell_state']
    adata.obs['final_state'] = adata.obs['cell_state']

    print("\nCell state summary:")
    print("  activated  = MI + high molecular score (true responders)")
    print("  bystander  = MI + low molecular score (non-responding)")
    print("  primed     = Normal + high molecular score (constitutively active)")
    print("  quiescent  = Normal + low molecular score (resting)")


def main():
    print("\n" + "#"*60)
    print("# Phase 2: Epicardial Cell State Classification")
    print("#"*60)

    # Load and merge individual epicardial files (Linna-Kuosmanen only)
    # NOTE: Kuppe excluded - investigation showed no true epicardial cells in dataset
    # (see README.md Section 3.1 for details)
    print("\nLoading individual epicardial datasets (Linna-Kuosmanen only)...")
    datasets = []

    # PERIHEART
    periheart = ad.read_h5ad(PROCESSED_DIR / "epicardial_periheart.h5ad")
    periheart.obs['dataset'] = 'PERIHEART'
    print(f"  PERIHEART: {periheart.n_obs:,} cells")
    datasets.append(periheart)

    # CAREBANK
    carebank = ad.read_h5ad(PROCESSED_DIR / "epicardial_carebank.h5ad")
    carebank.obs['dataset'] = 'CAREBANK'
    print(f"  CAREBANK: {carebank.n_obs:,} cells")
    datasets.append(carebank)

    # Find common genes and merge
    print("\nMerging datasets...")
    common_genes = set(datasets[0].var_names)
    for d in datasets[1:]:
        common_genes &= set(d.var_names)
    common_genes = sorted(list(common_genes))
    print(f"  Common genes: {len(common_genes):,}")

    # Save var metadata before subsetting (from first dataset)
    var_metadata = datasets[0].var.loc[common_genes].copy()

    for i in range(len(datasets)):
        datasets[i] = datasets[i][:, common_genes].copy()

    adata = ad.concat(datasets, join='outer')
    adata.obs_names_make_unique()

    # Restore var metadata (concat loses it)
    adata.var = var_metadata

    print(f"  Total: {adata.n_obs:,} cells, {adata.n_vars:,} genes")

    # Method 1: Condition-based
    assign_condition_labels(adata)

    # Method 2: Score-based (with improved sc.tl.score_genes + GMM threshold)
    calculate_molecular_scores(adata)
    classify_by_scores(adata)

    # Cross-validate
    cross_validate(adata)

    # Per-dataset summary
    print("\n" + "="*60)
    print("Per-Dataset Summary")
    print("="*60)

    for ds in ['PERIHEART', 'CAREBANK']:
        mask = adata.obs['dataset'] == ds
        subset = adata.obs[mask]
        print(f"\n[{ds}]")
        print(f"  Cell states: {dict(subset['cell_state'].value_counts())}")
        print(f"  Mean activation score: {subset['activation_score'].mean():.4f}")

    # Save
    output_path = PROCESSED_DIR / "epicardial_with_states.h5ad"
    print(f"\n{'='*60}")
    print(f"Saving to {output_path}...")
    adata.write_h5ad(output_path)

    print("\n" + "="*60)
    print("Done!")
    print("="*60 + "\n")

    return adata


if __name__ == "__main__":
    main()
