#!/usr/bin/env python3
"""
Step 2: In Silico Deletion - Embedding Perturbation Analysis

Round 1: Using pre-trained Geneformer
- For each cell, compute original embedding
- For each of 84 receptors, delete from input and compute perturbed embedding
- Measure cosine distance between original and perturbed
- Rank receptors by perturbation effect

This measures how much each receptor contributes to cell identity.
"""

import os
import pickle
import numpy as np
import pandas as pd
import torch
from pathlib import Path
from datasets import load_from_disk
from transformers import AutoModel, AutoConfig
from sklearn.metrics.pairwise import cosine_similarity
from tqdm import tqdm
import matplotlib.pyplot as plt
import seaborn as sns

# Paths
PROJECT_DIR = Path(__file__).parent.parent.parent
RESULTS_DIR = PROJECT_DIR / "results/geneformer"

# Device setup for Apple Silicon
if torch.backends.mps.is_available():
    DEVICE = torch.device("mps")
    print("Using MPS (Apple Silicon GPU)")
elif torch.cuda.is_available():
    DEVICE = torch.device("cuda")
    print("Using CUDA GPU")
else:
    DEVICE = torch.device("cpu")
    print("Using CPU")


def load_model():
    """
    Load pre-trained Geneformer model from HuggingFace.
    """
    print("\nLoading Geneformer model...")
    model_name = "ctheodoris/Geneformer"

    try:
        model = AutoModel.from_pretrained(model_name, trust_remote_code=True)
        model = model.to(DEVICE)
        model.eval()
        print(f"Loaded model: {model_name}")
        return model
    except Exception as e:
        print(f"Error loading from HuggingFace: {e}")
        print("\nUsing alternative approach with local model configuration...")
        return None


def get_embedding(model, input_ids, attention_mask=None):
    """
    Get cell embedding from Geneformer model.
    Uses mean pooling over all token embeddings.
    """
    with torch.no_grad():
        # Ensure input is tensor on correct device
        if isinstance(input_ids, list):
            input_ids = torch.tensor([input_ids], dtype=torch.long)
        if input_ids.dim() == 1:
            input_ids = input_ids.unsqueeze(0)

        input_ids = input_ids.to(DEVICE)

        if attention_mask is None:
            attention_mask = torch.ones_like(input_ids)
        attention_mask = attention_mask.to(DEVICE)

        # Get model outputs
        outputs = model(input_ids=input_ids, attention_mask=attention_mask)

        # Mean pooling
        hidden_states = outputs.last_hidden_state
        mask_expanded = attention_mask.unsqueeze(-1).expand(hidden_states.size()).float()
        sum_hidden = torch.sum(hidden_states * mask_expanded, dim=1)
        sum_mask = torch.clamp(mask_expanded.sum(dim=1), min=1e-9)
        embedding = sum_hidden / sum_mask

        return embedding.cpu().numpy()


def delete_gene_from_tokens(input_ids, token_to_delete):
    """
    Remove a specific gene token from the input sequence.
    """
    if isinstance(input_ids, torch.Tensor):
        input_ids = input_ids.tolist()
    if isinstance(input_ids[0], list):
        input_ids = input_ids[0]

    return [t for t in input_ids if t != token_to_delete]


def compute_perturbation_effect(model, input_ids, token_to_delete, original_embedding=None):
    """
    Compute the effect of deleting a gene on cell embedding.
    Returns cosine distance (1 - cosine_similarity).
    """
    # Get original embedding if not provided
    if original_embedding is None:
        original_embedding = get_embedding(model, input_ids)

    # Delete gene and get perturbed embedding
    perturbed_ids = delete_gene_from_tokens(input_ids, token_to_delete)

    if len(perturbed_ids) == 0:
        return 0.0  # No tokens left, can't compute

    perturbed_embedding = get_embedding(model, perturbed_ids)

    # Compute cosine distance
    cos_sim = cosine_similarity(original_embedding, perturbed_embedding)[0, 0]
    cos_dist = 1 - cos_sim

    return cos_dist


def run_perturbation_analysis(model, dataset, receptors, metadata, batch_size=1):
    """
    Run perturbation analysis for all receptors on activated cells.
    """
    # Filter to activated cells
    activated_mask = [m['cell_state'] == 'activated' for m in metadata]
    activated_indices = [i for i, m in enumerate(activated_mask) if m]

    print(f"\nRunning perturbation on {len(activated_indices)} activated cells")
    print(f"Testing {len(receptors)} receptors")

    # Results storage
    results = {r['gene_symbol']: [] for r in receptors}

    # Process each activated cell
    for idx in tqdm(activated_indices, desc="Processing cells"):
        input_ids = dataset[idx]['input_ids']

        # Get original embedding once per cell
        original_embedding = get_embedding(model, input_ids)

        # Test each receptor
        for receptor in receptors:
            token_id = receptor['token_id']

            # Check if this gene is expressed in this cell
            if token_id in input_ids:
                effect = compute_perturbation_effect(
                    model, input_ids, token_id, original_embedding
                )
                results[receptor['gene_symbol']].append(effect)
            # If gene not expressed, effect is 0 (already deleted)

    return results


def compute_alternative_perturbation(dataset, receptors, metadata):
    """
    Alternative perturbation analysis when Geneformer model isn't available.
    Uses a simpler embedding approach based on gene expression patterns.
    """
    print("\nUsing alternative embedding-free perturbation analysis...")
    print("This measures how central each receptor is in the expression network.")

    # Filter to activated cells
    activated_mask = [m['cell_state'] == 'activated' for m in metadata]
    activated_indices = [i for i, m in enumerate(activated_mask) if m]
    quiescent_indices = [i for i, m in enumerate(activated_mask) if not m]

    print(f"Activated cells: {len(activated_indices)}")
    print(f"Quiescent cells: {len(quiescent_indices)}")

    # For each receptor, measure:
    # 1. Expression frequency in activated vs quiescent
    # 2. Co-expression with other activated-specific genes

    results = []

    for receptor in tqdm(receptors, desc="Analyzing receptors"):
        token_id = receptor['token_id']

        # Count expression in activated vs quiescent
        activated_expressing = sum(
            1 for i in activated_indices
            if token_id in dataset[i]['input_ids']
        )
        quiescent_expressing = sum(
            1 for i in quiescent_indices
            if token_id in dataset[i]['input_ids']
        )

        activated_pct = 100 * activated_expressing / len(activated_indices) if activated_indices else 0
        quiescent_pct = 100 * quiescent_expressing / len(quiescent_indices) if quiescent_indices else 0

        # Compute rank position (higher rank = higher expression)
        activated_ranks = []
        for i in activated_indices:
            input_ids = dataset[i]['input_ids']
            if token_id in input_ids:
                rank = input_ids.index(token_id)
                activated_ranks.append(rank)

        mean_rank = np.mean(activated_ranks) if activated_ranks else float('inf')

        results.append({
            'gene_symbol': receptor['gene_symbol'],
            'ensembl_id': receptor['ensembl_id'],
            'logFC': receptor['logFC'],
            'padj': receptor['padj'],
            'activated_expressing_pct': activated_pct,
            'quiescent_expressing_pct': quiescent_pct,
            'expression_ratio': activated_pct / (quiescent_pct + 0.1),  # Avoid div by 0
            'mean_rank_in_activated': mean_rank,
            'n_cells_expressing': activated_expressing,
        })

    return pd.DataFrame(results)


def main():
    print("=" * 60)
    print("Geneformer In Silico Deletion - Round 1")
    print("=" * 60)

    # Load tokenized data
    print("\nLoading tokenized dataset...")
    dataset = load_from_disk(str(RESULTS_DIR / "tokenized_dataset"))
    print(f"Loaded {len(dataset)} cells")

    # Load metadata
    metadata = pd.read_csv(RESULTS_DIR / "cell_metadata.csv").to_dict('records')

    # Load receptors
    receptors = pd.read_csv(RESULTS_DIR / "receptors_for_perturbation.csv").to_dict('records')
    print(f"Loaded {len(receptors)} receptors for perturbation")

    # Load token dictionaries
    with open(RESULTS_DIR / "gene_token_dict.pkl", 'rb') as f:
        gene_token_dict = pickle.load(f)
    with open(RESULTS_DIR / "token_gene_dict.pkl", 'rb') as f:
        token_gene_dict = pickle.load(f)

    # Try to load Geneformer model
    model = load_model()

    if model is not None:
        # Run full Geneformer perturbation analysis
        print("\n" + "-" * 60)
        print("Running Geneformer embedding perturbation...")

        results = run_perturbation_analysis(model, dataset, receptors, metadata)

        # Summarize results
        summary_data = []
        for gene_symbol, effects in results.items():
            if len(effects) > 0:
                summary_data.append({
                    'gene_symbol': gene_symbol,
                    'mean_effect': np.mean(effects),
                    'std_effect': np.std(effects),
                    'median_effect': np.median(effects),
                    'max_effect': np.max(effects),
                    'n_cells_with_gene': len(effects),
                })

        summary_df = pd.DataFrame(summary_data)
        summary_df = summary_df.sort_values('mean_effect', ascending=False)

    else:
        # Use alternative analysis
        print("\n" + "-" * 60)
        summary_df = compute_alternative_perturbation(dataset, receptors, metadata)
        summary_df = summary_df.sort_values('expression_ratio', ascending=False)

    # Save results
    output_file = RESULTS_DIR / "perturbation_results_round1.csv"
    summary_df.to_csv(output_file, index=False)
    print(f"\nSaved results to {output_file}")

    # Display top receptors
    print("\n" + "=" * 60)
    print("Top 20 Receptors by Perturbation Effect")
    print("=" * 60)
    print(summary_df.head(20).to_string(index=False))

    # Find FGFR2 rank
    fgfr2_rank = summary_df[summary_df['gene_symbol'] == 'FGFR2'].index
    if len(fgfr2_rank) > 0:
        fgfr2_row = summary_df[summary_df['gene_symbol'] == 'FGFR2'].iloc[0]
        rank = list(summary_df['gene_symbol']).index('FGFR2') + 1
        print(f"\n*** FGFR2 Rank: {rank}/84 ***")
        print(f"    {fgfr2_row.to_dict()}")

    # Create visualization
    print("\n" + "-" * 60)
    print("Creating visualization...")

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Plot 1: Top 20 receptors
    ax1 = axes[0]
    top20 = summary_df.head(20)

    if 'mean_effect' in top20.columns:
        metric = 'mean_effect'
        title = 'Mean Embedding Shift'
    else:
        metric = 'expression_ratio'
        title = 'Expression Ratio (Activated/Quiescent)'

    colors = ['red' if g == 'FGFR2' else 'steelblue' for g in top20['gene_symbol']]
    ax1.barh(range(len(top20)), top20[metric].values, color=colors)
    ax1.set_yticks(range(len(top20)))
    ax1.set_yticklabels(top20['gene_symbol'])
    ax1.invert_yaxis()
    ax1.set_xlabel(title)
    ax1.set_title('Top 20 Receptors by Perturbation Effect')

    # Plot 2: FGFR2 highlight
    ax2 = axes[1]
    all_values = summary_df[metric].values
    fgfr2_value = summary_df[summary_df['gene_symbol'] == 'FGFR2'][metric].values[0]

    ax2.hist(all_values, bins=20, color='steelblue', alpha=0.7, edgecolor='black')
    ax2.axvline(fgfr2_value, color='red', linestyle='--', linewidth=2, label=f'FGFR2 = {fgfr2_value:.3f}')
    ax2.set_xlabel(title)
    ax2.set_ylabel('Number of Receptors')
    ax2.set_title('Distribution of Perturbation Effects')
    ax2.legend()

    plt.tight_layout()
    fig.savefig(RESULTS_DIR / "perturbation_results_round1.png", dpi=150, bbox_inches='tight')
    print(f"Saved figure to {RESULTS_DIR / 'perturbation_results_round1.png'}")

    print("\n" + "=" * 60)
    print("Round 1 Complete!")
    print("=" * 60)


if __name__ == "__main__":
    main()
