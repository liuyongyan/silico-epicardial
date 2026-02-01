#!/usr/bin/env python3
"""
DEG Analysis: Activated vs Quiescent Epicardial Cells

Uses pyDESeq2 to find differentially expressed genes between
activated and quiescent epicardial cells.

Usage:
    python 01_deg_analysis.py              # Full data
    python 01_deg_analysis.py --sample     # Use sample data
"""

import argparse
import anndata as ad
import numpy as np
import pandas as pd
from pathlib import Path
from scipy import sparse

PROJECT_DIR = Path(__file__).parent.parent.parent
PROCESSED_DIR = PROJECT_DIR / "data/processed"
RESULTS_DIR = PROJECT_DIR / "results/deg"


def run_deg_analysis(use_sample=False):
    """Run DEG analysis comparing activated vs quiescent epicardial cells."""

    # Create results directory
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    # Load data
    suffix = "_sample" if use_sample else ""
    input_file = PROCESSED_DIR / "epicardial_with_states.h5ad"

    print("="*60)
    print("DEG Analysis: Activated vs Quiescent Epicardial Cells")
    print("="*60)

    print(f"\nLoading {input_file}...")
    adata = ad.read_h5ad(input_file)
    print(f"  Total cells: {adata.n_obs:,}")
    print(f"  Total genes: {adata.n_vars:,}")

    # Filter to activated and quiescent cells only
    print("\nFiltering to activated vs quiescent cells...")
    mask = adata.obs['final_state'].isin(['activated', 'quiescent'])
    adata = adata[mask].copy()
    print(f"  Activated: {(adata.obs['final_state'] == 'activated').sum():,}")
    print(f"  Quiescent: {(adata.obs['final_state'] == 'quiescent').sum():,}")

    # Sample if requested (for faster development)
    if use_sample:
        print(f"\nSampling 10% for development...")
        np.random.seed(42)
        n_sample = int(adata.n_obs * 0.1)
        indices = np.random.choice(adata.n_obs, n_sample, replace=False)
        adata = adata[indices].copy()
        print(f"  Sampled cells: {adata.n_obs:,}")
        print(f"  Activated: {(adata.obs['final_state'] == 'activated').sum():,}")
        print(f"  Quiescent: {(adata.obs['final_state'] == 'quiescent').sum():,}")

    # Get raw counts (pyDESeq2 needs counts, not normalized data)
    # Since our data is log1p(CPM), we need to reverse it
    print("\nPreparing count matrix...")
    if sparse.issparse(adata.X):
        X = adata.X.toarray()
    else:
        X = adata.X.copy()

    # Reverse log1p: counts ≈ expm1(X)
    # Note: This gives pseudo-counts, not true counts
    counts = np.expm1(X)
    counts = np.round(counts).astype(int)
    counts = np.clip(counts, 0, None)  # Ensure non-negative

    # Filter low-expressed genes
    print("Filtering low-expressed genes...")
    min_cells = int(adata.n_obs * 0.01)  # At least 1% of cells
    gene_mask = (counts > 0).sum(axis=0) >= min_cells
    counts = counts[:, gene_mask]
    gene_names = adata.var_names[gene_mask]
    print(f"  Genes after filtering: {len(gene_names):,}")

    # Create count DataFrame
    count_df = pd.DataFrame(
        counts,
        index=adata.obs_names,
        columns=gene_names
    )

    # Create metadata DataFrame
    metadata = pd.DataFrame({
        'condition': adata.obs['final_state'].values
    }, index=adata.obs_names)

    print("\nRunning pyDESeq2...")
    try:
        from pydeseq2.dds import DeseqDataSet
        from pydeseq2.ds import DeseqStats

        # Create DESeq2 dataset
        dds = DeseqDataSet(
            counts=count_df,
            metadata=metadata,
            design_factors="condition",
            refit_cooks=True,
            n_cpus=4
        )

        # Run DESeq2
        dds.deseq2()

        # Get results (activated vs quiescent)
        stat_res = DeseqStats(dds, contrast=["condition", "activated", "quiescent"])
        stat_res.summary()

        # Get results DataFrame
        results_df = stat_res.results_df.copy()
        results_df['gene'] = results_df.index
        results_df = results_df.reset_index(drop=True)

        # Sort by adjusted p-value
        results_df = results_df.sort_values('padj')

        print(f"\nDEG Results Summary:")
        print(f"  Total genes tested: {len(results_df):,}")

        # Significant genes (padj < 0.05)
        sig_mask = results_df['padj'] < 0.05
        print(f"  Significant (padj < 0.05): {sig_mask.sum():,}")

        # Upregulated in activated (log2FC > 0)
        up_mask = sig_mask & (results_df['log2FoldChange'] > 0)
        print(f"  Upregulated in activated: {up_mask.sum():,}")

        # Downregulated in activated (log2FC < 0)
        down_mask = sig_mask & (results_df['log2FoldChange'] < 0)
        print(f"  Downregulated in activated: {down_mask.sum():,}")

        # Save full results
        output_file = RESULTS_DIR / f"deg_results{suffix}.csv"
        results_df.to_csv(output_file, index=False)
        print(f"\nSaved full results to: {output_file}")

        # Save upregulated genes (for NicheNet)
        upregulated = results_df[up_mask]['gene'].tolist()
        up_file = RESULTS_DIR / f"upregulated_genes{suffix}.txt"
        with open(up_file, 'w') as f:
            f.write('\n'.join(upregulated))
        print(f"Saved upregulated genes to: {up_file}")

        # Save downregulated genes
        downregulated = results_df[down_mask]['gene'].tolist()
        down_file = RESULTS_DIR / f"downregulated_genes{suffix}.txt"
        with open(down_file, 'w') as f:
            f.write('\n'.join(downregulated))
        print(f"Saved downregulated genes to: {down_file}")

        # Print top upregulated genes
        print("\nTop 20 Upregulated Genes in Activated Cells:")
        top_up = results_df[up_mask].head(20)
        for _, row in top_up.iterrows():
            print(f"  {row['gene']}: log2FC={row['log2FoldChange']:.2f}, padj={row['padj']:.2e}")

        # Check for expected EMT/activation markers
        print("\nChecking expected markers:")
        markers = ['VIM', 'CDH2', 'SNAI1', 'SNAI2', 'TWIST1', 'ZEB1', 'ZEB2',
                   'MKI67', 'PCNA', 'TOP2A', 'FN1', 'COL1A1', 'ACTA2']
        for marker in markers:
            if marker in results_df['gene'].values:
                row = results_df[results_df['gene'] == marker].iloc[0]
                status = "UP" if row['log2FoldChange'] > 0 else "DOWN"
                sig = "*" if row['padj'] < 0.05 else ""
                print(f"  {marker}: log2FC={row['log2FoldChange']:.2f}, padj={row['padj']:.2e} [{status}]{sig}")
            else:
                print(f"  {marker}: not found")

        return results_df

    except ImportError:
        print("ERROR: pyDESeq2 not installed. Install with: pip install pydeseq2")
        return None
    except Exception as e:
        print(f"ERROR: DEG analysis failed: {e}")
        import traceback
        traceback.print_exc()
        return None


def main():
    parser = argparse.ArgumentParser(description='DEG analysis for epicardial cells')
    parser.add_argument('--sample', action='store_true',
                        help='Use 10% sample for faster development')
    args = parser.parse_args()

    run_deg_analysis(use_sample=args.sample)


if __name__ == "__main__":
    main()
