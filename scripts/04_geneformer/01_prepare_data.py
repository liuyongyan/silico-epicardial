#!/usr/bin/env python3
"""
Step 1: Prepare Data for Geneformer Analysis

Converts epicardial scRNA-seq data to Geneformer's tokenized format.
Geneformer represents each cell as a sequence of genes ranked by expression.

Output:
- Tokenized dataset ready for Geneformer
- List of 84 upregulated receptors with their Ensembl IDs
- Verification that receptors are in Geneformer vocabulary
"""

import os
import pickle
import json
import numpy as np
import pandas as pd
import anndata as ad
from pathlib import Path
from scipy import sparse
from collections import OrderedDict
from datasets import Dataset
from tqdm import tqdm

# Paths
PROJECT_DIR = Path(__file__).parent.parent.parent
DATA_DIR = PROJECT_DIR / "data/processed"
RESULTS_DIR = PROJECT_DIR / "results/geneformer"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# Geneformer parameters
MAX_INPUT_SIZE = 2048  # Maximum number of genes per cell


def load_geneformer_token_dict():
    """
    Load or create Geneformer's gene token dictionary.
    Geneformer uses Ensembl IDs as gene identifiers.
    """
    # Try to load from Geneformer package
    try:
        from geneformer import TranscriptomeTokenizer
        tokenizer = TranscriptomeTokenizer()
        print("Loaded Geneformer tokenizer from package")
        return tokenizer.gene_token_dict, tokenizer
    except Exception as e:
        print(f"Could not load Geneformer tokenizer: {e}")
        print("Will create custom tokenization based on our gene set")
        return None, None


def rank_genes_by_expression(expression_vector, gene_names):
    """
    Rank genes by expression level (descending).
    Returns list of gene IDs ordered by expression.
    """
    # Get non-zero genes
    if sparse.issparse(expression_vector):
        expression_vector = expression_vector.toarray().flatten()

    # Get indices of non-zero genes, sorted by expression (descending)
    nonzero_mask = expression_vector > 0
    nonzero_indices = np.where(nonzero_mask)[0]

    if len(nonzero_indices) == 0:
        return []

    # Sort by expression (descending)
    sorted_indices = nonzero_indices[np.argsort(-expression_vector[nonzero_indices])]

    # Limit to max input size
    sorted_indices = sorted_indices[:MAX_INPUT_SIZE]

    # Return gene names (Ensembl IDs)
    return [gene_names[i] for i in sorted_indices]


def tokenize_cell(ranked_genes, gene_to_token):
    """
    Convert ranked gene list to token IDs.
    """
    tokens = []
    for gene in ranked_genes:
        if gene in gene_to_token:
            tokens.append(gene_to_token[gene])
    return tokens


def prepare_dataset(adata, gene_to_token, cell_states=['activated', 'quiescent']):
    """
    Prepare tokenized dataset for Geneformer.
    """
    print(f"\nPreparing dataset for cell states: {cell_states}")

    # Filter to desired cell states
    mask = adata.obs['cell_state'].isin(cell_states)
    adata_filtered = adata[mask].copy()
    print(f"Filtered to {adata_filtered.n_obs:,} cells")

    # Get gene names (Ensembl IDs)
    gene_names = list(adata_filtered.var_names)

    # Tokenize each cell
    tokenized_data = []
    cell_metadata = []

    print("Tokenizing cells...")
    for i in tqdm(range(adata_filtered.n_obs)):
        # Get expression vector
        expr = adata_filtered.X[i]

        # Rank genes by expression
        ranked_genes = rank_genes_by_expression(expr, gene_names)

        # Convert to tokens
        tokens = tokenize_cell(ranked_genes, gene_to_token)

        if len(tokens) > 0:
            tokenized_data.append({
                'input_ids': tokens,
                'length': len(tokens),
            })
            cell_metadata.append({
                'cell_id': adata_filtered.obs_names[i],
                'cell_state': adata_filtered.obs['cell_state'].iloc[i],
                'dataset': adata_filtered.obs['dataset'].iloc[i],
                'n_genes': len(ranked_genes),
                'n_tokens': len(tokens),
            })

    print(f"Successfully tokenized {len(tokenized_data):,} cells")

    return tokenized_data, cell_metadata


def load_upregulated_receptors():
    """
    Load the 84 upregulated receptors from DEG analysis.
    """
    deg_file = PROJECT_DIR / "results/deg/receptor_deg_LIANA.csv"
    df = pd.read_csv(deg_file)

    # Filter to significant upregulated
    up_receptors = df[(df['pvals_adj'] < 0.05) & (df['logfoldchanges'] > 0)].copy()
    up_receptors = up_receptors.sort_values('pvals_adj')

    print(f"\nLoaded {len(up_receptors)} upregulated receptors")

    # Create mapping
    receptor_info = []
    for _, row in up_receptors.iterrows():
        receptor_info.append({
            'ensembl_id': row['names'],
            'gene_symbol': row['gene'],
            'logFC': row['logfoldchanges'],
            'padj': row['pvals_adj'],
        })

    return receptor_info


def main():
    print("=" * 60)
    print("Geneformer Data Preparation")
    print("=" * 60)

    # Load epicardial data
    print("\nLoading epicardial data...")
    adata = ad.read_h5ad(DATA_DIR / "epicardial_with_states.h5ad")
    print(f"Loaded {adata.n_obs:,} cells, {adata.n_vars:,} genes")

    # Cell state distribution
    print("\nCell state distribution:")
    print(adata.obs['cell_state'].value_counts())

    # Load Geneformer token dictionary
    print("\n" + "-" * 60)
    print("Setting up tokenization...")
    gene_token_dict, tokenizer = load_geneformer_token_dict()

    if gene_token_dict is None:
        # Create our own token dictionary from the genes in our dataset
        print("Creating custom token dictionary from dataset genes...")
        gene_names = list(adata.var_names)
        gene_token_dict = {gene: i + 1 for i, gene in enumerate(gene_names)}  # 0 reserved for padding
        token_gene_dict = {i + 1: gene for i, gene in enumerate(gene_names)}
        print(f"Created dictionary with {len(gene_token_dict):,} genes")
    else:
        token_gene_dict = {v: k for k, v in gene_token_dict.items()}

    # Check overlap between our genes and Geneformer vocabulary
    our_genes = set(adata.var_names)
    vocab_genes = set(gene_token_dict.keys())
    overlap = our_genes & vocab_genes
    print(f"\nGene overlap with tokenizer vocabulary: {len(overlap):,}/{len(our_genes):,} ({100*len(overlap)/len(our_genes):.1f}%)")

    # Load upregulated receptors
    print("\n" + "-" * 60)
    receptors = load_upregulated_receptors()

    # Check which receptors are in vocabulary
    receptors_in_vocab = []
    receptors_missing = []
    for r in receptors:
        if r['ensembl_id'] in gene_token_dict:
            r['token_id'] = gene_token_dict[r['ensembl_id']]
            receptors_in_vocab.append(r)
        else:
            receptors_missing.append(r)

    print(f"\nReceptors in tokenizer vocabulary: {len(receptors_in_vocab)}/84")
    if receptors_missing:
        print(f"Missing receptors: {[r['gene_symbol'] for r in receptors_missing]}")

    # Save receptor info
    receptor_df = pd.DataFrame(receptors_in_vocab)
    receptor_df.to_csv(RESULTS_DIR / "receptors_for_perturbation.csv", index=False)
    print(f"Saved receptor info to {RESULTS_DIR / 'receptors_for_perturbation.csv'}")

    # Prepare tokenized dataset (activated + quiescent only)
    print("\n" + "-" * 60)
    tokenized_data, cell_metadata = prepare_dataset(
        adata,
        gene_token_dict,
        cell_states=['activated', 'quiescent']
    )

    # Create labels (0 = quiescent, 1 = activated)
    labels = [1 if m['cell_state'] == 'activated' else 0 for m in cell_metadata]

    # Add labels to tokenized data
    for i, td in enumerate(tokenized_data):
        td['label'] = labels[i]

    # Save as HuggingFace Dataset
    print("\n" + "-" * 60)
    print("Saving tokenized dataset...")

    # Convert to HuggingFace Dataset format
    dataset_dict = {
        'input_ids': [td['input_ids'] for td in tokenized_data],
        'length': [td['length'] for td in tokenized_data],
        'label': labels,
    }

    dataset = Dataset.from_dict(dataset_dict)
    dataset.save_to_disk(str(RESULTS_DIR / "tokenized_dataset"))
    print(f"Saved tokenized dataset to {RESULTS_DIR / 'tokenized_dataset'}")

    # Save cell metadata
    metadata_df = pd.DataFrame(cell_metadata)
    metadata_df['label'] = labels
    metadata_df.to_csv(RESULTS_DIR / "cell_metadata.csv", index=False)
    print(f"Saved cell metadata to {RESULTS_DIR / 'cell_metadata.csv'}")

    # Save token dictionaries
    with open(RESULTS_DIR / "gene_token_dict.pkl", 'wb') as f:
        pickle.dump(gene_token_dict, f)
    with open(RESULTS_DIR / "token_gene_dict.pkl", 'wb') as f:
        pickle.dump(token_gene_dict, f)
    print(f"Saved token dictionaries")

    # Summary statistics
    print("\n" + "=" * 60)
    print("Summary")
    print("=" * 60)

    activated_count = sum(labels)
    quiescent_count = len(labels) - activated_count

    print(f"Total cells: {len(tokenized_data):,}")
    print(f"  Activated: {activated_count:,}")
    print(f"  Quiescent: {quiescent_count:,}")
    print(f"\nTokens per cell:")
    lengths = [td['length'] for td in tokenized_data]
    print(f"  Mean: {np.mean(lengths):.0f}")
    print(f"  Median: {np.median(lengths):.0f}")
    print(f"  Range: {np.min(lengths)} - {np.max(lengths)}")
    print(f"\nReceptors for perturbation: {len(receptors_in_vocab)}")

    # Show FGFR2 info
    fgfr2_info = [r for r in receptors_in_vocab if r['gene_symbol'] == 'FGFR2']
    if fgfr2_info:
        print(f"\nFGFR2 info:")
        print(f"  Ensembl ID: {fgfr2_info[0]['ensembl_id']}")
        print(f"  Token ID: {fgfr2_info[0]['token_id']}")
        print(f"  logFC: {fgfr2_info[0]['logFC']:.3f}")
        print(f"  padj: {fgfr2_info[0]['padj']:.2e}")

    print("\n" + "=" * 60)
    print("Data preparation complete!")
    print("=" * 60)


if __name__ == "__main__":
    main()
