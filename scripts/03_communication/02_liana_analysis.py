#!/usr/bin/env python3
"""
LIANA Analysis: Ligand-Receptor Inference

Uses LIANA to identify ligand-receptor interactions between
sender cells and epicardial (receiver) cells.

Usage:
    python 02_liana_analysis.py              # Full data
    python 02_liana_analysis.py --sample     # Use sample data
    python 02_liana_analysis.py --kuppe      # Use Kuppe-only data
"""

import argparse
import anndata as ad
import numpy as np
import pandas as pd
from pathlib import Path

PROJECT_DIR = Path(__file__).parent.parent.parent
RAW_DIR = PROJECT_DIR / "data/raw"
PROCESSED_DIR = PROJECT_DIR / "data/processed"
RESULTS_DIR = PROJECT_DIR / "results/liana"


def get_ensembl_to_symbol_mapping():
    """Load Ensembl ID to gene symbol mapping from raw data."""
    print("Loading gene symbol mapping from raw data...")

    # Load from Kuppe raw data (has feature_name column)
    raw_file = RAW_DIR / "kuppe/54d24dbe-d39a-4844-bb21-07b5f4e173ad.h5ad"
    adata_raw = ad.read_h5ad(raw_file, backed='r')

    # Create mapping dictionary
    mapping = dict(zip(adata_raw.var_names, adata_raw.var['feature_name']))
    print(f"  Loaded {len(mapping):,} gene mappings")

    return mapping


def convert_to_gene_symbols(adata, mapping):
    """Convert var_names from Ensembl IDs to gene symbols."""
    print("Converting Ensembl IDs to gene symbols...")

    # Get new names
    new_names = []
    converted = 0
    for ensembl_id in adata.var_names:
        if ensembl_id in mapping:
            new_names.append(mapping[ensembl_id])
            converted += 1
        else:
            new_names.append(ensembl_id)  # Keep original if no mapping

    print(f"  Converted {converted:,}/{len(new_names):,} genes")

    # Handle duplicates by keeping first occurrence
    seen = set()
    unique_mask = []
    for name in new_names:
        if name not in seen:
            seen.add(name)
            unique_mask.append(True)
        else:
            unique_mask.append(False)

    # Subset to unique genes
    adata = adata[:, unique_mask].copy()
    adata.var_names = [new_names[i] for i, keep in enumerate(unique_mask) if keep]
    adata.var_names_make_unique()

    print(f"  Final gene count: {adata.n_vars:,}")
    return adata


def run_liana_analysis(use_sample=False, kuppe_only=False):
    """Run LIANA ligand-receptor analysis."""

    # Create results directory
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    # Determine input file
    if kuppe_only:
        suffix = "_kuppe_sample" if use_sample else "_kuppe"
        input_file = PROCESSED_DIR / f"communication_kuppe{'_sample' if use_sample else ''}.h5ad"
    else:
        suffix = "_merged_sample" if use_sample else "_merged"
        input_file = PROCESSED_DIR / f"communication_merged{'_sample' if use_sample else ''}.h5ad"

    print("="*60)
    print("LIANA Analysis: Ligand-Receptor Inference")
    print("="*60)

    print(f"\nLoading {input_file}...")
    adata = ad.read_h5ad(input_file)
    print(f"  Total cells: {adata.n_obs:,}")
    print(f"  Total genes: {adata.n_vars:,}")

    # Check if gene names are Ensembl IDs and convert if needed
    if adata.var_names[0].startswith('ENSG'):
        print("\nDetected Ensembl IDs - converting to gene symbols...")
        mapping = get_ensembl_to_symbol_mapping()
        adata = convert_to_gene_symbols(adata, mapping)

    print("\nCell type distribution:")
    print(adata.obs['sender_type'].value_counts())

    # Prepare for LIANA
    # LIANA expects a 'cell_type' column for groupby
    adata.obs['cell_type'] = adata.obs['sender_type'].copy()

    print("\nRunning LIANA...")
    try:
        import liana as li

        # Check available methods
        print("\nAvailable LIANA methods:")
        print(li.method.show_methods())

        # Run LIANA with rank_aggregate (consensus of multiple methods)
        # This integrates: CellPhoneDB, NATMI, Connectome, logFC, SingleCellSignalR, etc.
        li.mt.rank_aggregate(
            adata,
            groupby='cell_type',
            resource_name='consensus',  # Use consensus L-R database
            expr_prop=0.1,  # Minimum expression proportion
            verbose=True,
            use_raw=False,  # We already have normalized data
            n_perms=100 if use_sample else 1000,  # Fewer permutations for sample
        )

        # Get results
        liana_res = adata.uns['liana_res'].copy()
        print(f"\nLIANA Results: {len(liana_res):,} interactions")

        # Filter for interactions targeting Epicardial cells
        print("\nFiltering for Epicardial as receiver...")
        epi_res = liana_res[liana_res['target'] == 'Epicardial'].copy()
        print(f"  Epicardial-targeting interactions: {len(epi_res):,}")

        # Sort by magnitude_rank (lower is better)
        epi_res = epi_res.sort_values('magnitude_rank')

        # Save full LIANA results
        output_file = RESULTS_DIR / f"liana_all_results{suffix}.csv"
        liana_res.to_csv(output_file, index=False)
        print(f"\nSaved full results to: {output_file}")

        # Save epicardial-specific results
        epi_output = RESULTS_DIR / f"liana_epicardial_results{suffix}.csv"
        epi_res.to_csv(epi_output, index=False)
        print(f"Saved epicardial results to: {epi_output}")

        # Print top L-R pairs targeting epicardial cells
        print("\nTop 30 L-R Pairs Targeting Epicardial Cells:")
        print("-" * 80)
        top_pairs = epi_res.head(30)
        for i, row in top_pairs.iterrows():
            print(f"  {row['source']} → Epicardial: {row['ligand_complex']} → {row['receptor_complex']}")
            print(f"    magnitude_rank: {row['magnitude_rank']:.4f}, specificity_rank: {row.get('specificity_rank', 'N/A')}")

        # Aggregate by ligand (sum of interactions)
        print("\nTop Ligands (by interaction count):")
        ligand_counts = epi_res['ligand_complex'].value_counts().head(20)
        for ligand, count in ligand_counts.items():
            print(f"  {ligand}: {count} interactions")

        # Aggregate by sender cell type
        print("\nInteractions by Sender Cell Type:")
        sender_counts = epi_res['source'].value_counts()
        for sender, count in sender_counts.items():
            print(f"  {sender}: {count} interactions")

        # Check for expected ligands (FGF10, PDGFA, TGFB1, WNT5A, HGF)
        print("\nExpected Ligands Check:")
        expected_ligands = ['FGF10', 'PDGFA', 'TGFB1', 'WNT5A', 'HGF', 'PDGFB', 'IGF1', 'EGF']
        for ligand in expected_ligands:
            ligand_rows = epi_res[epi_res['ligand_complex'].str.contains(ligand, na=False)]
            if len(ligand_rows) > 0:
                best = ligand_rows.iloc[0]
                print(f"  {ligand}: Found! Best rank={best['magnitude_rank']:.4f}, "
                      f"receptor={best['receptor_complex']}, source={best['source']}")
            else:
                print(f"  {ligand}: Not found in top interactions")

        return liana_res, epi_res

    except ImportError:
        print("ERROR: LIANA not installed. Install with: pip install liana")
        return None, None
    except Exception as e:
        print(f"ERROR: LIANA analysis failed: {e}")
        import traceback
        traceback.print_exc()
        return None, None


def main():
    parser = argparse.ArgumentParser(description='LIANA ligand-receptor analysis')
    parser.add_argument('--sample', action='store_true',
                        help='Use sample data for faster development')
    parser.add_argument('--kuppe', action='store_true',
                        help='Use Kuppe-only data instead of merged')
    args = parser.parse_args()

    run_liana_analysis(use_sample=args.sample, kuppe_only=args.kuppe)


if __name__ == "__main__":
    main()
