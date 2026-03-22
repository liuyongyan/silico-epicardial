#!/bin/bash
# =============================================================================
# Run Round 2 after Round 1 completes
# Waits for Round 1 results, then runs fine-tuning and classifier perturbation
# =============================================================================

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$(dirname "$SCRIPT_DIR")")"
RESULTS_DIR="$PROJECT_DIR/results/geneformer"

echo "=============================================="
echo "Waiting for Round 1 to complete..."
echo "=============================================="

# Wait for Round 1 results file to appear
ROUND1_FILE="$RESULTS_DIR/perturbation_results_round1.csv"
while [ ! -f "$ROUND1_FILE" ]; do
    echo "$(date): Waiting for $ROUND1_FILE..."
    sleep 60
done

echo ""
echo "Round 1 complete! Results found at $ROUND1_FILE"
echo ""

# Activate conda environment
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate geneformer

cd "$SCRIPT_DIR"

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
echo "All analysis complete!"
echo ""
echo "Final results at:"
echo "  - $RESULTS_DIR/perturbation_results_round1.csv"
echo "  - $RESULTS_DIR/perturbation_results_round2.csv"
echo "  - $RESULTS_DIR/combined_rankings.csv"
echo "=============================================="
