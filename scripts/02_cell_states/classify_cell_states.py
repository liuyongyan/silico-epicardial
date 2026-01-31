#!/usr/bin/env python3
"""
Phase 2: Classify Epicardial Cell States

Two-method approach with cross-validation:
1. Condition-based: Using disease/spatial zone annotations
2. Score-based: Using EMT and proliferation molecular signatures

Output: Consensus classification (quiescent vs activated)
"""

import anndata as ad
import scanpy as sc
import numpy as np
import pandas as pd
from pathlib import Path
from scipy import sparse

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


def calculate_gene_score(adata, gene_dict, score_name):
    """Calculate mean expression score for a gene set."""
    available = get_available_genes(adata, gene_dict)

    if len(available) == 0:
        print(f"  WARNING: No genes found for {score_name}")
        adata.obs[score_name] = 0
        return 0

    gene_ids = list(available.values())
    gene_indices = [list(adata.var_names).index(g) for g in gene_ids]

    X = adata.X[:, gene_indices]
    if sparse.issparse(X):
        X = X.toarray()

    scores = X.mean(axis=1)
    adata.obs[score_name] = scores

    print(f"  {score_name}: {len(available)}/{len(gene_dict)} genes, "
          f"mean={scores.mean():.4f}, std={scores.std():.4f}")

    return len(available)


def assign_condition_labels(adata):
    """
    Method 1: Assign cell state labels based on condition annotations.

    Quiescent: normal/healthy tissue
    Activated: MI, ischemia, diseased tissue
    """
    print("\n" + "="*60)
    print("Method 1: Condition-based Classification")
    print("="*60)

    # Initialize
    adata.obs['condition_state'] = 'unknown'

    # Classification rules
    for idx in adata.obs.index:
        disease = adata.obs.loc[idx, 'disease']
        dataset = adata.obs.loc[idx, 'dataset']

        # Get spatial zone for Kuppe
        major_labl = adata.obs.loc[idx, 'major_labl'] if 'major_labl' in adata.obs.columns else None

        # Classify based on disease
        if disease == 'normal':
            if dataset == 'Kuppe_MI' and major_labl == 'CTRL':
                adata.obs.loc[idx, 'condition_state'] = 'quiescent'
            elif dataset in ['PERIHEART', 'CAREBANK']:
                adata.obs.loc[idx, 'condition_state'] = 'quiescent'
            else:
                adata.obs.loc[idx, 'condition_state'] = 'quiescent'
        elif disease in ['myocardial infarction', 'myocardial ischemia']:
            # For Kuppe, further classify by spatial zone
            if major_labl in ['IZ', 'BZ']:  # Ischemic/Border zone = high activation
                adata.obs.loc[idx, 'condition_state'] = 'activated_high'
            elif major_labl == 'FZ':  # Fibrotic zone = medium
                adata.obs.loc[idx, 'condition_state'] = 'activated_medium'
            elif major_labl == 'RZ':  # Remote zone = low activation
                adata.obs.loc[idx, 'condition_state'] = 'activated_low'
            else:
                adata.obs.loc[idx, 'condition_state'] = 'activated'
        elif disease in ['heart valve disorder', 'coronary artery disorder', 'heart failure']:
            adata.obs.loc[idx, 'condition_state'] = 'diseased'
        else:
            adata.obs.loc[idx, 'condition_state'] = 'unknown'

    # Simplify to binary for main analysis
    adata.obs['condition_binary'] = adata.obs['condition_state'].apply(
        lambda x: 'quiescent' if x == 'quiescent' else
                  ('activated' if 'activated' in x else 'other')
    )

    # Summary
    print("\nCondition state distribution:")
    print(adata.obs['condition_state'].value_counts())
    print("\nBinary classification:")
    print(adata.obs['condition_binary'].value_counts())


def calculate_molecular_scores(adata):
    """
    Method 2: Calculate EMT and proliferation scores.
    """
    print("\n" + "="*60)
    print("Method 2: Molecular Signature Scores")
    print("="*60)

    # Proliferation score
    print("\nProliferation signature:")
    calculate_gene_score(adata, PROLIFERATION_GENES, 'proliferation_score')

    # EMT up score
    print("\nEMT upregulated signature:")
    calculate_gene_score(adata, EMT_UP_GENES, 'emt_up_score')

    # EMT down score
    print("\nEMT downregulated (epithelial) signature:")
    calculate_gene_score(adata, EMT_DOWN_GENES, 'emt_down_score')

    # Combined EMT score
    adata.obs['emt_score'] = adata.obs['emt_up_score'] - adata.obs['emt_down_score']
    print(f"\nCombined EMT score: mean={adata.obs['emt_score'].mean():.4f}, "
          f"std={adata.obs['emt_score'].std():.4f}")

    # Combined activation score (average of normalized scores)
    prolif_z = (adata.obs['proliferation_score'] - adata.obs['proliferation_score'].mean()) / adata.obs['proliferation_score'].std()
    emt_z = (adata.obs['emt_score'] - adata.obs['emt_score'].mean()) / adata.obs['emt_score'].std()
    adata.obs['activation_score'] = (prolif_z + emt_z) / 2

    print(f"\nCombined activation score (z-normalized): "
          f"mean={adata.obs['activation_score'].mean():.4f}, "
          f"std={adata.obs['activation_score'].std():.4f}")


def classify_by_scores(adata, threshold_percentile=75):
    """
    Classify cells as activated based on molecular scores.
    """
    print("\n" + "="*60)
    print("Score-based Classification")
    print("="*60)

    threshold = np.percentile(adata.obs['activation_score'], threshold_percentile)
    print(f"Activation threshold ({threshold_percentile}th percentile): {threshold:.4f}")

    adata.obs['score_binary'] = adata.obs['activation_score'].apply(
        lambda x: 'activated' if x > threshold else 'quiescent'
    )

    print("\nScore-based classification:")
    print(adata.obs['score_binary'].value_counts())


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

    # Create consensus labels
    print("\n" + "-"*60)
    print("Creating consensus labels...")

    # Consensus: both methods agree
    adata.obs['consensus_state'] = 'ambiguous'

    # Both say quiescent
    quiescent_mask = (adata.obs['condition_binary'] == 'quiescent') & \
                     (adata.obs['score_binary'] == 'quiescent')
    adata.obs.loc[quiescent_mask, 'consensus_state'] = 'quiescent'

    # Both say activated
    activated_mask = (adata.obs['condition_binary'] == 'activated') & \
                     (adata.obs['score_binary'] == 'activated')
    adata.obs.loc[activated_mask, 'consensus_state'] = 'activated'

    # Condition says activated but score is low (possible early activation or resistant)
    early_mask = (adata.obs['condition_binary'] == 'activated') & \
                 (adata.obs['score_binary'] == 'quiescent')
    adata.obs.loc[early_mask, 'consensus_state'] = 'early_activated'

    # Score high but condition normal (possibly pre-disease or heterogeneous)
    predisease_mask = (adata.obs['condition_binary'] == 'quiescent') & \
                      (adata.obs['score_binary'] == 'activated')
    adata.obs.loc[predisease_mask, 'consensus_state'] = 'pre_activated'

    print("\nConsensus state distribution:")
    print(adata.obs['consensus_state'].value_counts())

    # Final binary for downstream analysis
    adata.obs['final_state'] = adata.obs['consensus_state'].apply(
        lambda x: 'activated' if x in ['activated', 'early_activated'] else
                  ('quiescent' if x == 'quiescent' else 'other')
    )

    print("\nFinal binary state:")
    print(adata.obs['final_state'].value_counts())


def main():
    print("\n" + "#"*60)
    print("# Phase 2: Epicardial Cell State Classification")
    print("#"*60)

    # Load data
    print("\nLoading merged epicardial dataset...")
    adata = ad.read_h5ad(PROCESSED_DIR / "epicardial_all_merged.h5ad")
    print(f"Loaded: {adata.n_obs:,} cells, {adata.n_vars:,} genes")

    # Method 1: Condition-based
    assign_condition_labels(adata)

    # Method 2: Score-based
    calculate_molecular_scores(adata)
    classify_by_scores(adata, threshold_percentile=75)

    # Cross-validate
    cross_validate(adata)

    # Per-dataset summary
    print("\n" + "="*60)
    print("Per-Dataset Summary")
    print("="*60)

    for ds in ['PERIHEART', 'CAREBANK', 'Kuppe_MI']:
        mask = adata.obs['dataset'] == ds
        subset = adata.obs[mask]
        print(f"\n[{ds}]")
        print(f"  Final state: {dict(subset['final_state'].value_counts())}")
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
