#!/usr/bin/env python3
"""
Improved Cell State Classification (PERIHEART Only)

Key improvements over original:
1. Uses sc.tl.score_genes() with control gene sets (background correction)
2. Uses Gaussian Mixture Model (GMM) for threshold selection instead of arbitrary percentile
3. PERIHEART dataset only (excludes CAREBANK for cleaner analysis)

Output: New classification files (does NOT modify existing files)
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
RESULTS_DIR = PROJECT_DIR / "results"
DEG_DIR = RESULTS_DIR / "deg"

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


def calculate_scores_improved(adata):
    """
    Calculate EMT and proliferation scores using sc.tl.score_genes().

    This method:
    - Compares gene set expression to a control gene set with similar expression levels
    - Provides background correction to account for technical variation
    - Is robust to differences in highly vs lowly expressed genes
    """
    print("\n" + "="*60)
    print("Improved Molecular Signature Scoring (sc.tl.score_genes)")
    print("="*60)

    # Proliferation score
    print("\nProliferation signature:")
    available_prolif = get_available_genes(adata, PROLIFERATION_GENES)
    prolif_genes = list(available_prolif.values())
    print(f"  Available genes: {len(prolif_genes)}/{len(PROLIFERATION_GENES)}")
    print(f"  Genes: {list(available_prolif.keys())}")

    sc.tl.score_genes(adata, gene_list=prolif_genes, score_name='proliferation_score_improved',
                      ctrl_size=50, n_bins=25, use_raw=False)
    print(f"  Score: mean={adata.obs['proliferation_score_improved'].mean():.4f}, "
          f"std={adata.obs['proliferation_score_improved'].std():.4f}")

    # EMT up score
    print("\nEMT upregulated signature:")
    available_emt_up = get_available_genes(adata, EMT_UP_GENES)
    emt_up_genes = list(available_emt_up.values())
    print(f"  Available genes: {len(emt_up_genes)}/{len(EMT_UP_GENES)}")
    print(f"  Genes: {list(available_emt_up.keys())}")

    sc.tl.score_genes(adata, gene_list=emt_up_genes, score_name='emt_up_score_improved',
                      ctrl_size=50, n_bins=25, use_raw=False)
    print(f"  Score: mean={adata.obs['emt_up_score_improved'].mean():.4f}, "
          f"std={adata.obs['emt_up_score_improved'].std():.4f}")

    # EMT down score
    print("\nEMT downregulated (epithelial) signature:")
    available_emt_down = get_available_genes(adata, EMT_DOWN_GENES)
    emt_down_genes = list(available_emt_down.values())
    print(f"  Available genes: {len(emt_down_genes)}/{len(EMT_DOWN_GENES)}")
    print(f"  Genes: {list(available_emt_down.keys())}")

    sc.tl.score_genes(adata, gene_list=emt_down_genes, score_name='emt_down_score_improved',
                      ctrl_size=50, n_bins=25, use_raw=False)
    print(f"  Score: mean={adata.obs['emt_down_score_improved'].mean():.4f}, "
          f"std={adata.obs['emt_down_score_improved'].std():.4f}")

    # Combined EMT score (up - down)
    adata.obs['emt_score_improved'] = adata.obs['emt_up_score_improved'] - adata.obs['emt_down_score_improved']
    print(f"\nCombined EMT score: mean={adata.obs['emt_score_improved'].mean():.4f}, "
          f"std={adata.obs['emt_score_improved'].std():.4f}")

    # Combined activation score (average of z-normalized scores)
    prolif_z = (adata.obs['proliferation_score_improved'] - adata.obs['proliferation_score_improved'].mean()) / adata.obs['proliferation_score_improved'].std()
    emt_z = (adata.obs['emt_score_improved'] - adata.obs['emt_score_improved'].mean()) / adata.obs['emt_score_improved'].std()
    adata.obs['activation_score_improved'] = (prolif_z + emt_z) / 2

    print(f"\nCombined activation score: "
          f"mean={adata.obs['activation_score_improved'].mean():.4f}, "
          f"std={adata.obs['activation_score_improved'].std():.4f}")


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


def classify_with_gmm(adata):
    """
    Classify cells using GMM-based threshold instead of arbitrary percentile.
    """
    print("\n" + "="*60)
    print("GMM-based Classification")
    print("="*60)

    scores = adata.obs['activation_score_improved']
    threshold, means, counts, gmm = find_threshold_gmm(scores)

    print(f"\nGMM analysis:")
    print(f"  Low component:  mean={means[np.argmin(means)]:.4f}, n={counts[0]:,}")
    print(f"  High component: mean={means[np.argmax(means)]:.4f}, n={counts[1]:,}")
    print(f"  Threshold (intersection): {threshold:.4f}")

    # For comparison, show what percentile this corresponds to
    percentile = (scores < threshold).mean() * 100
    print(f"  Equivalent percentile: {percentile:.1f}%")

    # Classify
    adata.obs['score_binary_improved'] = adata.obs['activation_score_improved'].apply(
        lambda x: 'activated' if x > threshold else 'quiescent'
    )

    print("\nScore-based classification:")
    print(adata.obs['score_binary_improved'].value_counts())

    return threshold


def assign_condition_labels(adata):
    """
    Assign condition labels based on disease annotation.
    """
    print("\n" + "="*60)
    print("Condition-based Classification")
    print("="*60)

    disease_map = {'normal': 'quiescent', 'myocardial infarction': 'activated'}
    adata.obs['condition_binary'] = adata.obs['disease'].astype(str).map(disease_map).fillna('unknown')

    print("\nCondition distribution:")
    print(adata.obs['condition_binary'].value_counts())


def create_cell_states(adata):
    """
    Create 4-state classification combining condition and molecular score.
    """
    print("\n" + "="*60)
    print("Creating Cell State Labels")
    print("="*60)

    adata.obs['cell_state_improved'] = 'unknown'

    # MI + high score = activated (true responders)
    activated_mask = (adata.obs['condition_binary'] == 'activated') & \
                     (adata.obs['score_binary_improved'] == 'activated')
    adata.obs.loc[activated_mask, 'cell_state_improved'] = 'activated'

    # MI + low score = bystander
    bystander_mask = (adata.obs['condition_binary'] == 'activated') & \
                     (adata.obs['score_binary_improved'] == 'quiescent')
    adata.obs.loc[bystander_mask, 'cell_state_improved'] = 'bystander'

    # Normal + high score = primed
    primed_mask = (adata.obs['condition_binary'] == 'quiescent') & \
                  (adata.obs['score_binary_improved'] == 'activated')
    adata.obs.loc[primed_mask, 'cell_state_improved'] = 'primed'

    # Normal + low score = quiescent
    quiescent_mask = (adata.obs['condition_binary'] == 'quiescent') & \
                     (adata.obs['score_binary_improved'] == 'quiescent')
    adata.obs.loc[quiescent_mask, 'cell_state_improved'] = 'quiescent'

    print("\nCell state distribution:")
    print(adata.obs['cell_state_improved'].value_counts())

    # Cross-tabulation
    print("\nCondition vs Score classification:")
    print(pd.crosstab(adata.obs['condition_binary'], adata.obs['score_binary_improved'], margins=True))


def run_deg_analysis(adata):
    """
    Run differential expression analysis: Normal vs MI cells.
    """
    print("\n" + "="*60)
    print("Differential Expression Analysis (Normal vs MI)")
    print("="*60)

    # Load receptor list
    receptor_file = RESULTS_DIR / "geneformer/receptors_for_perturbation.csv"
    receptors = pd.read_csv(receptor_file)
    receptor_ids = set(receptors['ensembl_id'].tolist())
    print(f"Loaded {len(receptor_ids)} receptors of interest")

    # Filter to cells with valid conditions
    mask = adata.obs['condition_binary'].isin(['quiescent', 'activated'])
    adata_filtered = adata[mask].copy()
    print(f"Cells for comparison: {adata_filtered.n_obs:,}")
    print(f"  Normal (quiescent): {(adata_filtered.obs['condition_binary'] == 'quiescent').sum():,}")
    print(f"  MI (activated): {(adata_filtered.obs['condition_binary'] == 'activated').sum():,}")

    # Run Wilcoxon rank-sum test
    print("\nRunning Wilcoxon rank-sum test...")
    sc.tl.rank_genes_groups(
        adata_filtered,
        groupby='condition_binary',
        groups=['activated'],  # Compare activated (MI) vs quiescent (Normal)
        reference='quiescent',
        method='wilcoxon',
        pts=True  # Include percent expressing
    )

    # Extract results
    result = adata_filtered.uns['rank_genes_groups']
    names = result['names']['activated']
    scores = result['scores']['activated']
    pvals = result['pvals']['activated']
    pvals_adj = result['pvals_adj']['activated']
    logfoldchanges = result['logfoldchanges']['activated']

    # pts keys vary by scanpy version
    if 'pts' in result:
        pts = result['pts']['activated']
        pts_rest = result['pts_rest']['activated'] if 'pts_rest' in result else None
    else:
        pts = None
        pts_rest = None

    # Build results dataframe
    deg_df = pd.DataFrame({
        'ensembl_id': names,
        'score': scores,
        'logFC': logfoldchanges,
        'pval': pvals,
        'padj': pvals_adj,
    })

    if pts is not None:
        deg_df['pct_MI'] = pts
    if pts_rest is not None:
        deg_df['pct_Normal'] = pts_rest

    # Filter to receptors only
    receptor_deg = deg_df[deg_df['ensembl_id'].isin(receptor_ids)].copy()

    # Add gene symbols
    symbol_map = dict(zip(receptors['ensembl_id'], receptors['gene_symbol']))
    receptor_deg['gene_symbol'] = receptor_deg['ensembl_id'].map(symbol_map)

    # Sort by adjusted p-value
    receptor_deg = receptor_deg.sort_values('padj')
    receptor_deg['rank'] = range(1, len(receptor_deg) + 1)

    # Reorder columns (only include columns that exist)
    base_cols = ['rank', 'gene_symbol', 'ensembl_id', 'logFC', 'padj', 'pval']
    optional_cols = ['pct_MI', 'pct_Normal', 'score']
    cols = base_cols + [c for c in optional_cols if c in receptor_deg.columns]
    receptor_deg = receptor_deg[cols]

    print(f"\nReceptor DEG results: {len(receptor_deg)} receptors")

    # Show top 20
    print("\nTop 20 receptors (sorted by p-value):")
    pd.set_option('display.max_columns', None)
    pd.set_option('display.width', None)
    print(receptor_deg.head(20).to_string(index=False))

    # Show FGFR2 specifically
    fgfr2 = receptor_deg[receptor_deg['gene_symbol'] == 'FGFR2']
    if len(fgfr2) > 0:
        print(f"\n*** FGFR2 ranking: {fgfr2['rank'].values[0]}/{len(receptor_deg)} ***")
        print(f"    logFC: {fgfr2['logFC'].values[0]:.4f}")
        print(f"    padj: {fgfr2['padj'].values[0]:.2e}")
        if 'pct_Normal' in fgfr2.columns and 'pct_MI' in fgfr2.columns:
            print(f"    % expressing: {fgfr2['pct_Normal'].values[0]*100:.1f}% (Normal) -> {fgfr2['pct_MI'].values[0]*100:.1f}% (MI)")

    return receptor_deg, adata_filtered


def main():
    print("\n" + "#"*60)
    print("# Improved Cell State Classification (PERIHEART Only)")
    print("#"*60)

    # Load PERIHEART data only
    print("\nLoading PERIHEART epicardial data...")
    adata = ad.read_h5ad(PROCESSED_DIR / "epicardial_periheart.h5ad")
    adata.obs['dataset'] = 'PERIHEART'
    print(f"  Cells: {adata.n_obs:,}")
    print(f"  Genes: {adata.n_vars:,}")

    # Show disease distribution
    print("\nDisease distribution:")
    print(adata.obs['disease'].value_counts())

    # Condition-based labels
    assign_condition_labels(adata)

    # Calculate improved molecular scores
    calculate_scores_improved(adata)

    # GMM-based classification
    threshold = classify_with_gmm(adata)

    # Create 4-state labels
    create_cell_states(adata)

    # Compare activation scores by condition
    print("\n" + "="*60)
    print("Score Comparison by Condition")
    print("="*60)
    for state in ['quiescent', 'activated']:
        mask = adata.obs['condition_binary'] == state
        if mask.sum() > 0:
            prolif = adata.obs.loc[mask, 'proliferation_score_improved'].mean()
            emt = adata.obs.loc[mask, 'emt_score_improved'].mean()
            act = adata.obs.loc[mask, 'activation_score_improved'].mean()
            label = 'Normal' if state == 'quiescent' else 'MI'
            print(f"  {label:8}: prolif={prolif:.4f}, emt={emt:.4f}, activation={act:.4f} (n={mask.sum():,})")

    # Run DEG analysis
    receptor_deg, adata_deg = run_deg_analysis(adata)

    # Save results to NEW files (don't modify existing)
    DEG_DIR.mkdir(parents=True, exist_ok=True)

    # Save DEG results
    output_deg = DEG_DIR / "receptor_normal_vs_mi_PERIHEART_improved.csv"
    receptor_deg.to_csv(output_deg, index=False)
    print(f"\nSaved receptor DEG results to: {output_deg}")

    # Save annotated data
    output_h5ad = PROCESSED_DIR / "epicardial_with_states_PERIHEART_improved.h5ad"
    adata.write_h5ad(output_h5ad)
    print(f"Saved annotated data to: {output_h5ad}")

    print("\n" + "="*60)
    print("Done!")
    print("="*60 + "\n")

    return adata, receptor_deg


if __name__ == "__main__":
    adata, receptor_deg = main()
