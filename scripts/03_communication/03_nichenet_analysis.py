#!/usr/bin/env python3
"""
NicheNet Analysis: Ligand Activity Prediction

Uses NicheNet to predict which ligands from sender cells
best explain gene expression changes in receiver (epicardial) cells.

NicheNet requires pre-computed networks from Zenodo:
- ligand_target_matrix.rds
- lr_network.rds
- weighted_networks.rds

Usage:
    python 03_nichenet_analysis.py              # Full data
    python 03_nichenet_analysis.py --sample     # Use sample data
"""

import argparse
import anndata as ad
import numpy as np
import pandas as pd
from pathlib import Path
from scipy import sparse

PROJECT_DIR = Path(__file__).parent.parent.parent
PROCESSED_DIR = PROJECT_DIR / "data/processed"
RESULTS_DIR = PROJECT_DIR / "results/nichenet"
DEG_DIR = PROJECT_DIR / "results/deg"
NICHENET_DATA_DIR = PROJECT_DIR / "data/nichenet"


def load_nichenet_networks():
    """Load NicheNet pre-computed networks."""
    print("\nLoading NicheNet networks...")

    # Check if networks exist
    ligand_target_file = NICHENET_DATA_DIR / "ligand_target_matrix.csv"
    lr_network_file = NICHENET_DATA_DIR / "lr_network.csv"

    if not ligand_target_file.exists() or not lr_network_file.exists():
        print("ERROR: NicheNet networks not found.")
        print(f"Expected files in: {NICHENET_DATA_DIR}")
        print("\nTo download NicheNet data:")
        print("  1. Go to https://zenodo.org/record/7074291")
        print("  2. Download ligand_target_matrix.rds and lr_network.rds")
        print("  3. Convert to CSV using R:")
        print("     ligand_target_matrix <- readRDS('ligand_target_matrix.rds')")
        print("     write.csv(as.matrix(ligand_target_matrix), 'ligand_target_matrix.csv')")
        print("     lr_network <- readRDS('lr_network.rds')")
        print("     write.csv(lr_network, 'lr_network.csv', row.names=FALSE)")
        return None, None

    print(f"  Loading ligand-target matrix from {ligand_target_file}...")
    ligand_target = pd.read_csv(ligand_target_file, index_col=0)
    print(f"    Shape: {ligand_target.shape}")

    print(f"  Loading L-R network from {lr_network_file}...")
    lr_network = pd.read_csv(lr_network_file)
    print(f"    Interactions: {len(lr_network):,}")

    return ligand_target, lr_network


def get_expressed_genes(adata, cell_type, min_pct=0.1):
    """Get genes expressed in a specific cell type."""
    mask = adata.obs['sender_type'] == cell_type
    subset = adata[mask]

    if sparse.issparse(subset.X):
        X = subset.X.toarray()
    else:
        X = subset.X

    # Calculate percentage of cells expressing each gene
    pct_expressed = (X > 0).sum(axis=0) / X.shape[0]
    expressed_genes = adata.var_names[pct_expressed >= min_pct].tolist()

    return expressed_genes


def run_nichenet_analysis(use_sample=False):
    """Run NicheNet ligand activity prediction."""

    # Create results directory
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    NICHENET_DATA_DIR.mkdir(parents=True, exist_ok=True)

    suffix = "_sample" if use_sample else ""

    print("="*60)
    print("NicheNet Analysis: Ligand Activity Prediction")
    print("="*60)

    # Load NicheNet networks
    ligand_target, lr_network = load_nichenet_networks()
    if ligand_target is None:
        return None

    # Load DEG results
    deg_file = DEG_DIR / f"upregulated_genes{suffix}.txt"
    if not deg_file.exists():
        # Try without suffix
        deg_file = DEG_DIR / "upregulated_genes.txt"
    if not deg_file.exists():
        deg_file = DEG_DIR / "upregulated_genes_sample.txt"

    if not deg_file.exists():
        print(f"ERROR: DEG results not found. Run 01_deg_analysis.py first.")
        return None

    print(f"\nLoading DEG results from {deg_file}...")
    with open(deg_file, 'r') as f:
        geneset = [line.strip() for line in f if line.strip()]
    print(f"  Upregulated genes: {len(geneset)}")

    # Load communication dataset
    input_file = PROCESSED_DIR / f"communication_merged{'_sample' if use_sample else ''}.h5ad"
    if not input_file.exists():
        input_file = PROCESSED_DIR / "communication_kuppe.h5ad"

    print(f"\nLoading {input_file}...")
    adata = ad.read_h5ad(input_file)
    print(f"  Total cells: {adata.n_obs:,}")

    # Get expressed genes in receivers (epicardial)
    print("\nIdentifying expressed genes...")
    receiver_genes = get_expressed_genes(adata, 'Epicardial', min_pct=0.1)
    print(f"  Receiver (Epicardial) expressed genes: {len(receiver_genes)}")

    # Get potential ligands from sender cell types
    sender_types = ['Cardiomyocyte', 'Fibroblast', 'Endothelial', 'Myeloid', 'Macrophage']
    sender_genes = set()
    for sender in sender_types:
        if sender in adata.obs['sender_type'].values:
            genes = get_expressed_genes(adata, sender, min_pct=0.1)
            sender_genes.update(genes)
            print(f"  {sender} expressed genes: {len(genes)}")
    sender_genes = list(sender_genes)
    print(f"  Total sender genes: {len(sender_genes)}")

    # Filter ligand-receptor network for expressed genes
    print("\nFiltering L-R network for expressed genes...")
    lr_filtered = lr_network[
        (lr_network['from'].isin(sender_genes)) &
        (lr_network['to'].isin(receiver_genes))
    ]
    print(f"  Filtered L-R pairs: {len(lr_filtered)}")

    # Get potential ligands
    potential_ligands = lr_filtered['from'].unique().tolist()
    potential_ligands = [l for l in potential_ligands if l in ligand_target.columns]
    print(f"  Potential ligands (in NicheNet database): {len(potential_ligands)}")

    if len(potential_ligands) == 0:
        print("ERROR: No potential ligands found in NicheNet database")
        return None

    # Filter geneset to genes in ligand-target matrix
    geneset_filtered = [g for g in geneset if g in ligand_target.index]
    print(f"  Geneset genes in NicheNet database: {len(geneset_filtered)}")

    if len(geneset_filtered) == 0:
        print("ERROR: No geneset genes found in NicheNet database")
        return None

    # Calculate ligand activity scores
    print("\nCalculating ligand activity scores...")

    # Ligand activity = mean regulatory potential of ligand for target genes
    ligand_activities = []
    for ligand in potential_ligands:
        # Get regulatory potential for this ligand
        target_scores = ligand_target.loc[geneset_filtered, ligand]
        activity = target_scores.mean()
        ligand_activities.append({
            'ligand': ligand,
            'activity': activity,
            'n_targets': (target_scores > 0).sum(),
            'mean_target_score': target_scores[target_scores > 0].mean() if (target_scores > 0).sum() > 0 else 0
        })

    # Create results DataFrame
    results_df = pd.DataFrame(ligand_activities)
    results_df = results_df.sort_values('activity', ascending=False)
    results_df['rank'] = range(1, len(results_df) + 1)

    print(f"\nNicheNet Results: {len(results_df)} ligands scored")

    # Save results
    output_file = RESULTS_DIR / f"nichenet_ligand_activities{suffix}.csv"
    results_df.to_csv(output_file, index=False)
    print(f"\nSaved results to: {output_file}")

    # Print top ligands
    print("\nTop 30 Ligands by Activity Score:")
    print("-" * 60)
    for _, row in results_df.head(30).iterrows():
        print(f"  {row['rank']:3d}. {row['ligand']:15s} activity={row['activity']:.4f}, "
              f"targets={row['n_targets']}")

    # Check for expected ligands
    print("\nExpected Ligands Check:")
    expected = ['FGF10', 'PDGFA', 'TGFB1', 'WNT5A', 'HGF', 'PDGFB', 'IGF1', 'EGF', 'FGF2', 'BMP4']
    for ligand in expected:
        if ligand in results_df['ligand'].values:
            row = results_df[results_df['ligand'] == ligand].iloc[0]
            print(f"  {ligand}: rank={row['rank']}, activity={row['activity']:.4f}")
        else:
            print(f"  {ligand}: Not found in results")

    # Find receptors for top ligands
    print("\nReceptors for Top 10 Ligands:")
    for _, row in results_df.head(10).iterrows():
        ligand = row['ligand']
        receptors = lr_filtered[lr_filtered['from'] == ligand]['to'].unique()
        print(f"  {ligand}: {', '.join(receptors[:5])}")

    return results_df


def main():
    parser = argparse.ArgumentParser(description='NicheNet ligand activity analysis')
    parser.add_argument('--sample', action='store_true',
                        help='Use sample data for faster development')
    args = parser.parse_args()

    run_nichenet_analysis(use_sample=args.sample)


if __name__ == "__main__":
    main()
