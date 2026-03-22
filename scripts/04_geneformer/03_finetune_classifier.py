#!/usr/bin/env python3
"""
Step 3: Fine-tune Geneformer for Cell State Classification

Fine-tunes Geneformer to classify:
- Activated epicardial cells (label=1)
- Quiescent epicardial cells (label=0)

This teaches the model what distinguishes activated vs quiescent states,
which is essential for meaningful perturbation analysis.
"""

import os
import pickle
import numpy as np
import pandas as pd
import torch
from pathlib import Path
from datasets import load_from_disk, Dataset
from transformers import (
    AutoModel,
    AutoModelForSequenceClassification,
    AutoConfig,
    TrainingArguments,
    Trainer,
    DataCollatorWithPadding,
)
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, precision_recall_fscore_support, confusion_matrix
import matplotlib.pyplot as plt
import seaborn as sns

# Paths
PROJECT_DIR = Path(__file__).parent.parent.parent
RESULTS_DIR = PROJECT_DIR / "results/geneformer"
MODEL_DIR = RESULTS_DIR / "finetuned_model"
MODEL_DIR.mkdir(parents=True, exist_ok=True)

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


class CellStateDataCollator:
    """
    Custom data collator for variable-length token sequences.
    Pads sequences to the same length within a batch.
    """
    def __init__(self, pad_token_id=0, max_length=2048):
        self.pad_token_id = pad_token_id
        self.max_length = max_length

    def __call__(self, features):
        # Get max length in this batch
        max_len = min(
            max(len(f['input_ids']) for f in features),
            self.max_length
        )

        # Pad sequences
        batch = {
            'input_ids': [],
            'attention_mask': [],
            'labels': [],
        }

        for f in features:
            seq = f['input_ids'][:max_len]
            padding_length = max_len - len(seq)

            batch['input_ids'].append(seq + [self.pad_token_id] * padding_length)
            batch['attention_mask'].append([1] * len(seq) + [0] * padding_length)
            batch['labels'].append(f['label'])

        # Convert to tensors
        batch['input_ids'] = torch.tensor(batch['input_ids'], dtype=torch.long)
        batch['attention_mask'] = torch.tensor(batch['attention_mask'], dtype=torch.long)
        batch['labels'] = torch.tensor(batch['labels'], dtype=torch.long)

        return batch


def compute_metrics(eval_pred):
    """
    Compute evaluation metrics.
    """
    predictions, labels = eval_pred
    predictions = np.argmax(predictions, axis=1)

    accuracy = accuracy_score(labels, predictions)
    precision, recall, f1, _ = precision_recall_fscore_support(
        labels, predictions, average='binary'
    )

    return {
        'accuracy': accuracy,
        'precision': precision,
        'recall': recall,
        'f1': f1,
    }


def create_simple_classifier(input_dim, hidden_dim=256, num_labels=2):
    """
    Create a simple classifier for when Geneformer isn't available.
    Uses embedding of ranked genes.
    """
    import torch.nn as nn

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
            # Embedding
            x = self.embedding(input_ids)

            # LSTM
            lstm_out, (h_n, c_n) = self.lstm(x)

            # Use final hidden states
            hidden = torch.cat((h_n[-2], h_n[-1]), dim=1)

            # Classify
            logits = self.classifier(hidden)

            loss = None
            if labels is not None:
                loss_fn = nn.CrossEntropyLoss()
                loss = loss_fn(logits, labels)

            return type('Output', (), {'loss': loss, 'logits': logits})()

    return SimpleClassifier

def main():
    print("=" * 60)
    print("Fine-tuning Cell State Classifier")
    print("=" * 60)

    # Load tokenized data
    print("\nLoading tokenized dataset...")
    dataset = load_from_disk(str(RESULTS_DIR / "tokenized_dataset"))
    print(f"Loaded {len(dataset)} cells")

    # Load metadata
    metadata = pd.read_csv(RESULTS_DIR / "cell_metadata.csv")
    print(f"\nLabel distribution:")
    print(metadata['cell_state'].value_counts())

    # Load token dictionary for vocab size
    with open(RESULTS_DIR / "gene_token_dict.pkl", 'rb') as f:
        gene_token_dict = pickle.load(f)
    vocab_size = len(gene_token_dict) + 1  # +1 for padding token

    # Split into train/val/test
    print("\n" + "-" * 60)
    print("Splitting data...")

    # Get indices
    indices = list(range(len(dataset)))
    labels = [dataset[i]['label'] for i in indices]

    # Stratified split
    train_idx, temp_idx = train_test_split(
        indices, test_size=0.3, stratify=labels, random_state=42
    )
    temp_labels = [labels[i] for i in temp_idx]
    val_idx, test_idx = train_test_split(
        temp_idx, test_size=0.5, stratify=temp_labels, random_state=42
    )

    print(f"Train: {len(train_idx)} cells")
    print(f"Val: {len(val_idx)} cells")
    print(f"Test: {len(test_idx)} cells")

    # Create dataset splits
    train_data = dataset.select(train_idx)
    val_data = dataset.select(val_idx)
    test_data = dataset.select(test_idx)

    # Check class balance
    train_labels = [train_data[i]['label'] for i in range(len(train_data))]
    print(f"\nTrain label distribution: {sum(train_labels)} activated, {len(train_labels) - sum(train_labels)} quiescent")

    # Try to load Geneformer for fine-tuning
    print("\n" + "-" * 60)
    print("Loading model for fine-tuning...")

    try:
        model = AutoModelForSequenceClassification.from_pretrained(
            "ctheodoris/Geneformer",
            num_labels=2,
            trust_remote_code=True
        )
        model = model.to(DEVICE)
        print("Loaded Geneformer for fine-tuning")
        use_geneformer = True
    except Exception as e:
        print(f"Could not load Geneformer: {e}")
        print("Using simple LSTM classifier instead")

        # Create simple model
        SimpleClassifier = create_simple_classifier(vocab_size, hidden_dim=256, num_labels=2)
        model = SimpleClassifier(vocab_size, embed_dim=128, hidden_dim=256, num_labels=2)
        model = model.to(DEVICE)
        use_geneformer = False

    # Data collator
    collator = CellStateDataCollator(pad_token_id=0)

    # Training arguments
    # Note: use_mps_device is deprecated - transformers auto-detects MPS
    # Using small batch sizes due to MPS memory constraints
    training_args = TrainingArguments(
        output_dir=str(MODEL_DIR),
        num_train_epochs=5,  # Reduced epochs
        per_device_train_batch_size=4,  # Reduced from 32 for MPS memory
        per_device_eval_batch_size=8,   # Reduced from 64 for MPS memory
        gradient_accumulation_steps=8,  # Effective batch size = 32
        warmup_steps=50,
        weight_decay=0.01,
        logging_dir=str(MODEL_DIR / "logs"),
        logging_steps=50,
        eval_strategy="epoch",
        save_strategy="epoch",
        load_best_model_at_end=True,
        metric_for_best_model="f1",
        greater_is_better=True,
        dataloader_num_workers=0,  # MPS works better with 0 workers
        report_to="none",
        fp16=False,  # MPS doesn't support fp16
        bf16=False,
    )

    # Trainer
    trainer = Trainer(
        model=model,
        args=training_args,
        train_dataset=train_data,
        eval_dataset=val_data,
        data_collator=collator,
        compute_metrics=compute_metrics,
    )

    # Train
    print("\n" + "-" * 60)
    print("Training classifier...")
    trainer.train()

    # Evaluate on test set
    print("\n" + "-" * 60)
    print("Evaluating on test set...")
    test_results = trainer.evaluate(test_data)
    print(f"\nTest Results:")
    for key, value in test_results.items():
        print(f"  {key}: {value:.4f}")

    # Get predictions for confusion matrix
    predictions = trainer.predict(test_data)
    pred_labels = np.argmax(predictions.predictions, axis=1)
    true_labels = [test_data[i]['label'] for i in range(len(test_data))]

    # Confusion matrix
    cm = confusion_matrix(true_labels, pred_labels)
    print(f"\nConfusion Matrix:")
    print(f"  TN={cm[0,0]}, FP={cm[0,1]}")
    print(f"  FN={cm[1,0]}, TP={cm[1,1]}")

    # Save model
    print("\n" + "-" * 60)
    print("Saving model...")

    if use_geneformer:
        model.save_pretrained(str(MODEL_DIR / "best_model"))
    else:
        torch.save(model.state_dict(), MODEL_DIR / "best_model.pt")

    # Save model info
    model_info = {
        'use_geneformer': use_geneformer,
        'vocab_size': vocab_size,
        'test_accuracy': test_results['eval_accuracy'],
        'test_f1': test_results['eval_f1'],
    }
    with open(MODEL_DIR / "model_info.pkl", 'wb') as f:
        pickle.dump(model_info, f)

    print(f"Saved model to {MODEL_DIR}")

    # Create visualization
    print("\n" + "-" * 60)
    print("Creating visualizations...")

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Plot 1: Confusion matrix
    ax1 = axes[0]
    sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', ax=ax1,
                xticklabels=['Quiescent', 'Activated'],
                yticklabels=['Quiescent', 'Activated'])
    ax1.set_xlabel('Predicted')
    ax1.set_ylabel('True')
    ax1.set_title('Confusion Matrix (Test Set)')

    # Plot 2: Metrics bar chart
    ax2 = axes[1]
    metrics = ['Accuracy', 'Precision', 'Recall', 'F1']
    values = [
        test_results['eval_accuracy'],
        test_results['eval_precision'],
        test_results['eval_recall'],
        test_results['eval_f1'],
    ]
    colors = ['steelblue'] * 4
    ax2.bar(metrics, values, color=colors, edgecolor='black')
    ax2.set_ylim(0, 1)
    ax2.set_ylabel('Score')
    ax2.set_title('Classification Performance')
    for i, v in enumerate(values):
        ax2.text(i, v + 0.02, f'{v:.3f}', ha='center')

    plt.tight_layout()
    fig.savefig(RESULTS_DIR / "classifier_performance.png", dpi=150, bbox_inches='tight')
    print(f"Saved figure to {RESULTS_DIR / 'classifier_performance.png'}")

    print("\n" + "=" * 60)
    print("Fine-tuning Complete!")
    print(f"Test Accuracy: {test_results['eval_accuracy']:.4f}")
    print(f"Test F1: {test_results['eval_f1']:.4f}")
    print("=" * 60)


if __name__ == "__main__":
    main()
