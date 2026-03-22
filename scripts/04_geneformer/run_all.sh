#!/bin/bash
# =============================================================================
# Run Complete Geneformer In Silico Deletion Analysis
# =============================================================================

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

echo "=============================================="
echo "Geneformer In Silico Deletion Analysis"
echo "=============================================="
echo ""

# Check if conda environment exists
if ! conda env list | grep -q "geneformer"; then
    echo "Environment 'geneformer' not found."
    echo "Please run setup first: bash 00_setup_environment.sh"
    exit 1
fi

# Activate environment
echo "Activating geneformer environment..."
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate geneformer

echo ""
echo "=============================================="
echo "Step 1: Prepare Data"
echo "=============================================="
python 01_prepare_data.py

echo ""
echo "=============================================="
echo "Step 2: Round 1 - Embedding Perturbation"
echo "=============================================="
python 02_embedding_perturbation.py

echo ""
echo "=============================================="
echo "Step 3: Fine-tune Classifier"
echo "=============================================="
python 03_finetune_classifier.py

echo ""
echo "=============================================="
echo "Step 4: Round 2 - Classifier Perturbation"
echo "=============================================="
python 04_classifier_perturbation.py

echo ""
echo "=============================================="
echo "Step 5: Compare Results"
echo "=============================================="
python 05_compare_results.py

echo ""
echo "=============================================="
echo "All steps complete!"
echo ""
echo "Results saved to: results/geneformer/"
echo "  - perturbation_results_round1.csv"
echo "  - perturbation_results_round2.csv"
echo "  - combined_rankings.csv"
echo "  - *.png visualizations"
echo "=============================================="
