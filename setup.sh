#!/bin/bash

# NanoporeToBED Pipeline - Fast Environment Setup
# Only sets up the conda environment, no other files

set -e

echo "=========================================="
echo "NanoporeToBED - Fast Environment Setup"
echo "=========================================="
echo ""

# Check if running on HPC (common module commands)
if command -v module &> /dev/null; then
    echo "Detected HPC environment with module system"
    echo "Attempting to load conda module..."
    module load conda 2>/dev/null || module load miniconda3 2>/dev/null || module load anaconda3 2>/dev/null || true
fi

# Function to detect best conda variant
detect_conda() {
    if command -v mamba &> /dev/null; then
        echo "mamba"
    elif command -v micromamba &> /dev/null; then
        echo "micromamba"
    elif command -v conda &> /dev/null; then
        echo "conda"
    else
        echo "none"
    fi
}

CONDA_CMD=$(detect_conda)

if [[ "$CONDA_CMD" == "none" ]]; then
    echo "[ERROR] No conda/mamba/micromamba found"
    echo ""
    echo "Quick install options:"
    echo ""
    echo "1. Install micromamba (fastest, recommended):"
    echo "   curl -L https://micro.mamba.pm/install.sh | bash"
    echo ""
    echo "2. Install mamba in existing conda:"
    echo "   conda install -n base -c conda-forge mamba"
    echo ""
    echo "3. Install miniconda:"
    echo "   wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh"
    echo "   bash Miniconda3-latest-Linux-x86_64.sh"
    exit 1
fi

echo "[OK] Using: $CONDA_CMD"
echo ""

# Environment name
ENV_NAME="nanopore_methylation"

# Check if environment already exists
if $CONDA_CMD env list 2>/dev/null | grep -q "^$ENV_NAME "; then
    echo "[WARNING] Environment '$ENV_NAME' already exists"
    read -p "Remove and recreate? (y/n): " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        echo "Removing existing environment..."
        $CONDA_CMD env remove -n $ENV_NAME -y 2>/dev/null || true
    else
        echo "Keeping existing environment"
        echo ""
        echo "[SUCCESS] Setup complete (using existing environment)"
        echo ""
        echo "Activate with: $CONDA_CMD activate $ENV_NAME"
        exit 0
    fi
fi

echo "=========================================="
echo "Creating Environment: $ENV_NAME"
echo "=========================================="
echo ""

# Install based on available tool
if [[ "$CONDA_CMD" == "mamba" ]] || [[ "$CONDA_CMD" == "micromamba" ]]; then
    echo ">> Fast installation with $CONDA_CMD"
    echo "   This should take 2-5 minutes..."
    echo ""
    
    $CONDA_CMD create -y -n $ENV_NAME \
        -c conda-forge -c bioconda -c defaults \
        python=3.9 \
        samtools=1.18 \
        minimap2=2.26 \
        ont-modkit=0.3.0 \
        qualimap=2.3 \
        pigz=2.8 \
        parallel=20230722 \
        bedtools=2.31.0 \
        bcftools=1.18 \
        --no-channel-priority
    
    if [ $? -eq 0 ]; then
        echo ""
        echo "[OK] Environment created successfully"
    else
        echo ""
        echo "[ERROR] Environment creation failed"
        echo ""
        echo "Try manual installation (see below)"
        exit 1
    fi

elif [[ "$CONDA_CMD" == "conda" ]]; then
    echo ">> Using conda (slower solver)"
    echo "   Consider installing mamba for faster setup:"
    echo "   conda install -n base -c conda-forge mamba"
    echo ""
    echo "   This may take 10-20 minutes..."
    echo ""
    
    # Create base environment
    conda create -y -n $ENV_NAME python=3.9
    
    # Activate and install in stages
    eval "$(conda shell.bash hook)"
    conda activate $ENV_NAME
    
    echo "Installing core tools..."
    conda install -y -c conda-forge -c bioconda \
        samtools=1.18 minimap2=2.26 pigz parallel
    
    echo "Installing bioinformatics tools..."
    conda install -y -c bioconda \
        ont-modkit bedtools bcftools qualimap
    
    conda deactivate
    
    echo ""
    echo "[OK] Environment created successfully"
fi

echo ""
echo "=========================================="
echo "[SUCCESS] Environment Setup Complete!"
echo "=========================================="
echo ""
echo "=========================================="
echo "NEXT STEP: Activate the environment"
echo "=========================================="
echo ""
echo "Run this command now:"
echo ""
echo "  $CONDA_CMD activate $ENV_NAME"
echo ""
echo "=========================================="
echo ""
echo "To verify installation:"
echo "  samtools --version"
echo "  minimap2 --version"
echo "  modkit --version"
echo ""
echo "=========================================="
echo ""

# Create quick reference file
cat > environment_setup_complete.txt << EOF
NanoporeToBED Environment Setup Complete
========================================

Environment name: $ENV_NAME
Installation date: $(date)
Conda tool used: $CONDA_CMD

Activate with:
  $CONDA_CMD activate $ENV_NAME

Installed tools:
  - samtools 1.18
  - minimap2 2.26
  - ont-modkit 0.3.0
  - qualimap 2.3
  - bedtools 2.31.0
  - bcftools 1.18
  - pigz
  - parallel

Citation:
  Drag et al. (2025) bioRxiv 2025.04.11.648151
EOF

echo "Setup info saved to: environment_setup_complete.txt"
echo ""
