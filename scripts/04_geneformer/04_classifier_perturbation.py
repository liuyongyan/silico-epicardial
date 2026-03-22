#!/usr/bin/env python3
"""
Step 4: In Silico Deletion - Classifier Perturbation Analysis

Round 2: Using fine-tuned cell state classifier
- For each activated cell, get P(activated) from classifier
- For each receptor, delete it and get new P(activated)
- Measure the drop in P(activated)
- Rank receptors by how much their deletion reduces activated probability

This directly answers: "Which receptors are most critical for maintaining
the activated epicardial phenotype?"
"""

import os
import pickle
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from pathlib import Path
from datasets import load_from_disk
from transformers import AutoModelForSequenceClassification
from tqdm import tqdm
import matplotlib.pyplot as plt
import seaborn as sns

# Paths
PROJECT_DIR = Path(__file__).parent.parent.parent
RESULTS_DIR = PROJECT_DIR / "results/geneformer"
MODEL_DIR = RESULTS_DIR / "finetuned_model"

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


def create_simple_classifier(vocab_size, embed_dim=128, hidden_dim=256, num_labels=2):
    """
    Recreate the simple classifier architecture.
    """
    class SimpleClassifier(nn.Module):
        def __init__(self, vocab_size, embed_dim, hidden_dim, num_labels):
            super().__init__()
            self.embedding = nn.Embedding(vocab_size, embed_dim, padding_idx=0)
            self.lstm = nn.LSTM(embed_dim, hidden_dim, batch_first=True, bidirectional=True)
            self.classifier = nn.Sequential(
                nn.Linear(hidden_dim * 2, hidden_dim),
                nn.ReLU(),
                nn.Dropout(0.1),
                nn.Linear(hidden_dim, num_labels)
            )

        def forward(self, input_ids, attention_mask=None, labels=None):
            x = self.embedding(input_ids)
            lstm_out, (h_n, c_n) = self.lstm(x)
            hidden = torch.cat((h_n[-2], h_n[-1]), dim=1)
            logits = self.classifier(hidden)

            loss = None
            if labels is not None:
                loss_fn = nn.CrossEntropyLoss()
                loss = loss_fn(logits, labels)

            return type('Output', (), {'loss': loss, 'logits': logits})()

    return SimpleClassifier(vocab_size, embed_dim, hidden_dim, num_labels)


def load_model():
    """
    Load the fine-tuned classifier.
    """
    print("\nLoading fine-tuned model...")

    # Load model info
    with open(MODEL_DIR / "model_info.pkl", 'rb') as f:
        model_info = pickle.load(f)

    if model_info['use_geneformer']:
        model = AutoModelForSequenceClassification.from_pretrained(
            str(MODEL_DIR / "best_model"),
            trust_remote_code=True
        )
    else:
        model = create_simple_classifier(
            model_info['vocab_size'],
            embed_dim=128,
            hidden_dim=256,
            num_labels=2
        )
        model.load_state_dict(torch.load(MODEL_DIR / "best_model.pt", map_location=DEVICE))

    model = model.to(DEVICE)
    model.eval()

    print(f"Loaded model (Geneformer: {model_info['use_geneformer']})")
    print(f"  Test accuracy: {model_info['test_accuracy']:.4f}")
    print(f"  Test F1: {model_info['test_f1']:.4f}")

    return model


def get_activated_probability(model, input_ids, max_length=2048):
    """
    Get the probability of being classified as 'activated'.
    """
    with torch.no_grad():
        # Prepare input
        if isinstance(input_ids, list):
            input_ids = input_ids[:max_length]
            input_ids = torch.tensor([input_ids], dtype=torch.long)

        input_ids = input_ids.to(DEVICE)

        # Get model output
        outputs = model(input_ids=input_ids)
        logits = outputs.logits

        # Softmax to get probabilities
        probs = torch.softmax(logits, dim=1)

        # Return P(activated) - class 1
        return probs[0, 1].item()


def delete_gene_from_tokens(input_ids, token_to_delete):
    """
    Remove a specific gene token from the input sequence.
    """
    if isinstance(input_ids, torch.Tensor):
        input_ids = input_ids.tolist()
    if isinstance(input_ids[0], list):
        input_ids = input_ids[0]

    return [t for t in input_ids if t != token_to_delete]


def run_perturbation_analysis(model, dataset, receptors, metadata):
    """
    Run perturbation analysis for all receptors on activated cells.

    For each receptor:
    1. Delete it from each activated cell
    2. Measure the drop in P(activated)
    3. Average across cells

    Higher drop = more important for maintaining activated state
    """
    # Filter to activated cells
    activated_indices = [
        i for i, m in enumerate(metadata)
        if m['cell_state'] == 'activated'
    ]

    print(f"\nRunning perturbation on {len(activated_indices)} activated cells")
    print(f"Testing {len(receptors)} receptors")

    # Results storage
    results = {r['gene_symbol']: {
        'prob_drops': [],
        'original_probs': [],
        'perturbed_probs': [],
        'cells_expressing': 0,
        'classification_flips': 0,
    } for r in receptors}

    # Create receptor token lookup
    receptor_tokens = {r['gene_symbol']: r['token_id'] for r in receptors}

    # Process each activated cell
    for idx in tqdm(activated_indices, desc="Processing cells"):
        input_ids = dataset[idx]['input_ids']

        # Get original P(activated)
        original_prob = get_activated_probability(model, input_ids)

        # Test each receptor
        for receptor in receptors:
            gene = receptor['gene_symbol']
            token_id = receptor['token_id']

            # Check if this gene is expressed in this cell
            if token_id in input_ids:
                results[gene]['cells_expressing'] += 1

                # Delete gene and get new probability
                perturbed_ids = delete_gene_from_tokens(input_ids, token_id)

                if len(perturbed_ids) > 0:
                    perturbed_prob = get_activated_probability(model, perturbed_ids)

                    # Calculate probability drop
                    prob_drop = original_prob - perturbed_prob

                    results[gene]['prob_drops'].append(prob_drop)
                    results[gene]['original_probs'].append(original_prob)
                    results[gene]['perturbed_probs'].append(perturbed_prob)

                    # Check if classification flipped
                    if original_prob >= 0.5 and perturbed_prob < 0.5:
                        results[gene]['classification_flips'] += 1

    return results


def summarize_results(results, receptors):
    """
    Summarize perturbation results into a DataFrame.
    """
    summary_data = []

    for receptor in receptors:
        gene = receptor['gene_symbol']
        r = results[gene]

        if len(r['prob_drops']) > 0:
            summary_data.append({
                'gene_symbol': gene,
                'ensembl_id': receptor['ensembl_id'],
                'logFC': receptor['logFC'],
                'padj': receptor['padj'],
                'mean_prob_drop': np.mean(r['prob_drops']),
                'std_prob_drop': np.std(r['prob_drops']),
                'median_prob_drop': np.median(r['prob_drops']),
                'max_prob_drop': np.max(r['prob_drops']),
                'mean_original_prob': np.mean(r['original_probs']),
                'mean_perturbed_prob': np.mean(r['perturbed_probs']),
                'cells_expressing': r['cells_expressing'],
                'classification_flips': r['classification_flips'],
                'flip_rate': r['classification_flips'] / r['cells_expressing'] if r['cells_expressing'] > 0 else 0,
            })
        else:
            summary_data.append({
                'gene_symbol': gene,
                'ensembl_id': receptor['ensembl_id'],
                'logFC': receptor['logFC'],
                'padj': receptor['padj'],
                'mean_prob_drop': 0,
                'std_prob_drop': 0,
                'median_prob_drop': 0,
                'max_prob_drop': 0,
                'mean_original_prob': 0,
                'mean_perturbed_prob': 0,
                'cells_expressing': 0,
                'classification_flips': 0,
                'flip_rate': 0,
            })

    return pd.DataFrame(summary_data)


def main():
    print("=" * 60)
    print("Geneformer In Silico Deletion - Round 2")
    print("(Classifier Perturbation Analysis)")
    print("=" * 60)

    # Load model
    model = load_model()

    # Load tokenized data
    print("\nLoading tokenized dataset...")
    dataset = load_from_disk(str(RESULTS_DIR / "tokenized_dataset"))
    print(f"Loaded {len(dataset)} cells")

    # Load metadata
    metadata = pd.read_csv(RESULTS_DIR / "cell_metadata.csv").to_dict('records')

    # Load receptors
    receptors = pd.read_csv(RESULTS_DIR / "receptors_for_perturbation.csv").to_dict('records')
    print(f"Loaded {len(receptors)} receptors for perturbation")

    # Run perturbation analysis
    print("\n" + "-" * 60)
    print("Running classifier perturbation analysis...")
    results = run_perturbation_analysis(model, dataset, receptors, metadata)

    # Summarize results
    summary_df = summarize_results(results, receptors)

    # Sort by mean probability drop (most important first)
    summary_df = summary_df.sort_values('mean_prob_drop', ascending=False)
    summary_df['rank'] = range(1, len(summary_df) + 1)

    # Save results
    output_file = RESULTS_DIR / "perturbation_results_round2.csv"
    summary_df.to_csv(output_file, index=False)
    print(f"\nSaved results to {output_file}")

    # Display results
    print("\n" + "=" * 60)
    print("Top 20 Receptors by Probability Drop")
    print("(Higher = more important for activated state)")
    print("=" * 60)

    display_cols = ['rank', 'gene_symbol', 'mean_prob_drop', 'flip_rate', 'cells_expressing', 'logFC']
    print(summary_df[display_cols].head(20).to_string(index=False))

    # Find FGFR2 rank
    fgfr2_row = summary_df[summary_df['gene_symbol'] == 'FGFR2']
    if len(fgfr2_row) > 0:
        fgfr2_rank = fgfr2_row['rank'].values[0]
        print(f"\n" + "=" * 60)
        print(f"*** FGFR2 Rank: {fgfr2_rank}/84 ***")
        print("=" * 60)
        print(fgfr2_row[display_cols].to_string(index=False))

    # Also show FGFR1 for comparison
    fgfr1_row = summary_df[summary_df['gene_symbol'] == 'FGFR1']
    if len(fgfr1_row) > 0:
        print(f"\nFGFR1 for comparison:")
        print(fgfr1_row[display_cols].to_string(index=False))

    # Create visualizations
    print("\n" + "-" * 60)
    print("Creating visualizations...")

    fig, axes = plt.subplots(2, 2, figsize=(14, 12))

    # Plot 1: Top 20 by probability drop
    ax1 = axes[0, 0]
    top20 = summary_df.head(20)
    colors = ['red' if g == 'FGFR2' else ('orange' if g == 'FGFR1' else 'steelblue')
              for g in top20['gene_symbol']]
    ax1.barh(range(len(top20)), top20['mean_prob_drop'].values, color=colors)
    ax1.set_yticks(range(len(top20)))
    ax1.set_yticklabels(top20['gene_symbol'])
    ax1.invert_yaxis()
    ax1.set_xlabel('Mean P(activated) Drop')
    ax1.set_title('Top 20 Receptors by Perturbation Effect')
    ax1.axvline(0, color='black', linestyle='-', linewidth=0.5)

    # Plot 2: Distribution with FGFR2 highlighted
    ax2 = axes[0, 1]
    all_drops = summary_df['mean_prob_drop'].values
    fgfr2_drop = summary_df[summary_df['gene_symbol'] == 'FGFR2']['mean_prob_drop'].values[0]

    ax2.hist(all_drops, bins=20, color='steelblue', alpha=0.7, edgecolor='black')
    ax2.axvline(fgfr2_drop, color='red', linestyle='--', linewidth=2,
                label=f'FGFR2 = {fgfr2_drop:.4f}')
    ax2.set_xlabel('Mean P(activated) Drop')
    ax2.set_ylabel('Number of Receptors')
    ax2.set_title('Distribution of Perturbation Effects')
    ax2.legend()

    # Plot 3: Probability drop vs logFC
    ax3 = axes[1, 0]
    ax3.scatter(summary_df['logFC'], summary_df['mean_prob_drop'],
                alpha=0.6, s=50, c='steelblue')

    # Highlight FGFR2 and FGFR1
    fgfr2_data = summary_df[summary_df['gene_symbol'] == 'FGFR2']
    fgfr1_data = summary_df[summary_df['gene_symbol'] == 'FGFR1']
    ax3.scatter(fgfr2_data['logFC'], fgfr2_data['mean_prob_drop'],
                color='red', s=100, zorder=5, label='FGFR2')
    ax3.scatter(fgfr1_data['logFC'], fgfr1_data['mean_prob_drop'],
                color='orange', s=100, zorder=5, label='FGFR1')

    ax3.set_xlabel('Log Fold Change')
    ax3.set_ylabel('Mean P(activated) Drop')
    ax3.set_title('Perturbation Effect vs Expression Change')
    ax3.legend()

    # Plot 4: Flip rate vs probability drop
    ax4 = axes[1, 1]
    ax4.scatter(summary_df['mean_prob_drop'], summary_df['flip_rate'],
                alpha=0.6, s=50, c='steelblue')
    ax4.scatter(fgfr2_data['mean_prob_drop'], fgfr2_data['flip_rate'],
                color='red', s=100, zorder=5, label='FGFR2')
    ax4.scatter(fgfr1_data['mean_prob_drop'], fgfr1_data['flip_rate'],
                color='orange', s=100, zorder=5, label='FGFR1')
    ax4.set_xlabel('Mean P(activated) Drop')
    ax4.set_ylabel('Classification Flip Rate')
    ax4.set_title('Flip Rate vs Perturbation Effect')
    ax4.legend()

    plt.tight_layout()
    fig.savefig(RESULTS_DIR / "perturbation_results_round2.png", dpi=150, bbox_inches='tight')
    print(f"Saved figure to {RESULTS_DIR / 'perturbation_results_round2.png'}")

    # Summary statistics
    print("\n" + "=" * 60)
    print("Summary Statistics")
    print("=" * 60)
    print(f"Receptors tested: {len(summary_df)}")
    print(f"Mean prob drop across all receptors: {summary_df['mean_prob_drop'].mean():.4f}")
    print(f"Receptors with positive effect (drop > 0): {(summary_df['mean_prob_drop'] > 0).sum()}")
    print(f"Receptors causing classification flips: {(summary_df['flip_rate'] > 0).sum()}")

    print("\n" + "=" * 60)
    print("Round 2 Complete!")
    print("=" * 60)


if __name__ == "__main__":
    main()
