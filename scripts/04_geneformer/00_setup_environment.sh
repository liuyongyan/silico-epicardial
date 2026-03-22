#!/bin/bash
# =============================================================================
# Geneformer Environment Setup for Apple Silicon (M4 Max)
# =============================================================================

set -e

ENV_NAME="geneformer"

echo "=============================================="
echo "Setting up Geneformer environment: $ENV_NAME"
echo "=============================================="

# Check if conda is available
if ! command -v conda &> /dev/null; then
    echo "ERROR: conda not found. Please install Anaconda or Miniconda first."
    exit 1
fi

# Create conda environment
echo ""
echo "Step 1: Creating conda environment..."
conda create -n $ENV_NAME python=3.10 -y

# Activate environment
echo ""
echo "Step 2: Activating environment..."
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate $ENV_NAME

# Install PyTorch with MPS support (Apple Silicon)
echo ""
echo "Step 3: Installing PyTorch with MPS support..."
pip install torch torchvision torchaudio

# Install core dependencies
echo ""
echo "Step 4: Installing core dependencies..."
pip install transformers>=4.28.0
pip install datasets>=2.12.0
pip install accelerate>=0.20.0
pip install scikit-learn>=1.2.0
pip install scipy>=1.10.0
pip install numpy>=1.24.0
pip install pandas>=2.0.0
pip install anndata>=0.9.0
pip install scanpy>=1.9.0
pip install matplotlib>=3.7.0
pip install seaborn>=0.12.0
pip install tqdm>=4.65.0

# Install Geneformer from HuggingFace
echo ""
echo "Step 5: Installing Geneformer..."
pip install git+https://huggingface.co/ctheodoris/Geneformer.git

# Verify installation
echo ""
echo "Step 6: Verifying installation..."
python -c "
import torch
print(f'PyTorch version: {torch.__version__}')
print(f'MPS available: {torch.backends.mps.is_available()}')
print(f'MPS built: {torch.backends.mps.is_built()}')

from transformers import AutoModel, AutoTokenizer
print('Transformers imported successfully')

try:
    from geneformer import TranscriptomeTokenizer
    print('Geneformer imported successfully')
except ImportError as e:
    print(f'Geneformer import warning: {e}')
    print('Will use manual tokenization approach')
"

echo ""
echo "=============================================="
echo "Setup complete!"
echo ""
echo "To activate the environment, run:"
echo "  conda activate $ENV_NAME"
echo ""
echo "To verify MPS (Apple GPU) support:"
echo "  python -c \"import torch; print(torch.backends.mps.is_available())\""
echo "=============================================="
