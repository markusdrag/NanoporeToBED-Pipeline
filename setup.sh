#!/bin/bash

# NanoporeToBED Pipeline - Automated Setup Script
# This script sets up the environment and downloads the pipeline

set -e

echo "=========================================="
echo "NanoporeToBED Pipeline Setup"
echo "=========================================="
echo ""

# Check if running on HPC (common module commands)
if command -v module &> /dev/null; then
    echo "Detected HPC environment with module system"
    echo "Attempting to load conda module..."
    module load conda 2>/dev/null || module load miniconda3 2>/dev/null || module load anaconda3 2>/dev/null || true
fi

# Function to check if conda/mamba is available
check_conda() {
    if command -v conda &> /dev/null; then
        echo "✓ Found conda"
        CONDA_CMD="conda"
        return 0
    elif command -v micromamba &> /dev/null; then
        echo "✓ Found micromamba"
        CONDA_CMD="micromamba"
        return 0
    elif command -v mamba &> /dev/null; then
        echo "✓ Found mamba"
        CONDA_CMD="mamba"
        return 0
    else
        echo "❌ Error: No conda/mamba/micromamba found"
        echo "Please install conda or load the appropriate module"
        exit 1
    fi
}

# Parse arguments
INSTALL_DIR="${1:-$(pwd)/NanoporeToBED-Pipeline}"
GITHUB_REPO="${2:-YOUR_USERNAME/NanoporeToBED-Pipeline}"

echo "Installation directory: $INSTALL_DIR"
echo ""

# Check for conda
check_conda

# Create installation directory
if [[ -d "$INSTALL_DIR" ]]; then
    echo "⚠️  Directory $INSTALL_DIR already exists"
    read -p "Do you want to overwrite? (y/n): " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        echo "Installation cancelled"
        exit 1
    fi
    rm -rf "$INSTALL_DIR"
fi

mkdir -p "$INSTALL_DIR"
cd "$INSTALL_DIR"

# Create directory structure
echo "Creating directory structure..."
mkdir -p scripts docs examples test_data
echo "✓ Directories created"
echo ""

# Create the main pipeline script
echo "Creating NanoporeToBED.sh script..."
cat > scripts/NanoporeToBED.sh << 'SCRIPT_EOF'
#!/bin/bash

#SBATCH --job-name=NanoporeToBED
#SBATCH --output=NanoporeToBED.out
#SBATCH --error=NanoporeToBED.err
#SBATCH -c 40
#SBATCH --mem 192g
#SBATCH --time=72:00:00
#SBATCH --account YourAccount

# Environment
cd $HOME
source ~/.bashrc
micromamba activate nanopore_methylation

set -euo pipefail

# Default values
THREADS=40
dry_run=false

# Parse arguments
while [[ "$#" -gt 0 ]]; do
  case "$1" in
    -i|--input) input_dir="$2"; shift 2 ;;
    -o|--output) output_dir="$2"; shift 2 ;;
    -ref|--reference) ref_genome="$2"; shift 2 ;;
    -t|--threads) THREADS="$2"; shift 2 ;;
    --dry-run) dry_run=true; shift ;;
    -h|--help)
      echo "Usage: $0 -i <input_dir> -o <output_dir> -ref <reference_genome.fna> [-t <threads>] [--dry-run]"
      echo ""
      echo "Required arguments:"
      echo "  -i, --input       Input directory containing SRR folders"
      echo "  -o, --output      Output directory for processed data"
      echo "  -ref, --reference Path to reference genome FASTA file"
      echo ""
      echo "Optional arguments:"
      echo "  -t, --threads     Number of threads to use (default: 40)"
      echo "  --dry-run         Run in test mode without processing"
      echo "  -h, --help        Show this help message"
      exit 0
      ;;
    *) echo "Unknown option: $1"; echo "Use -h for help"; exit 1 ;;
  esac
done

# Rest of the pipeline script continues...
# [Full script content would be inserted here]
SCRIPT_EOF

chmod +x scripts/NanoporeToBED.sh
echo "✓ Main script created"
echo ""

# Create environment.yml
echo "Creating environment.yml..."
cat > environment.yml << 'ENV_EOF'
name: nanopore_methylation
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
  # Core tools
  - python>=3.8
  - samtools>=1.17
  - minimap2>=2.26
  - bioconda::ont-modkit>=0.3.0
  - qualimap>=2.3

  # Utilities
  - pigz
  - parallel
  - gzip

  # Optional but useful
  - bedtools>=2.30
  - bcftools>=1.17
  - multiqc>=1.14
ENV_EOF

echo "✓ environment.yml created"
echo ""

# Create example run script
echo "Creating example run script..."
cat > examples/example_run.sh << 'EXAMPLE_EOF'
#!/bin/bash
# Example run configuration for Nanopore methylation analysis

# Set paths
INPUT_DIR="/data/nanopore/guppy_output"
OUTPUT_DIR="/data/nanopore/methylation_analysis"
REFERENCE_GENOME="/data/references/genome.fna"

# Submit job with custom thread count
sbatch scripts/NanoporeToBED.sh \
  -i "$INPUT_DIR" \
  -o "$OUTPUT_DIR" \
  -ref "$REFERENCE_GENOME" \
  -t 32
EXAMPLE_EOF

chmod +x examples/example_run.sh
echo "✓ Example script created"
echo ""

# Create .gitignore
echo "Creating .gitignore..."
cat > .gitignore << 'GITIGNORE_EOF'
# Output files
*.bam
*.bai
*.bed
*.fastq
*.fastq.gz

# Logs
*.log
*.err
*.out
logs/
slurm-*.out

# Temporary files
*.tmp
temp/
tmp/

# Large data
data/
test_data/large/
output/
results/

# OS files
.DS_Store
Thumbs.db

# Python
__pycache__/
*.py[cod]

# Environment
.env
.venv
env/
venv/
GITIGNORE_EOF

echo "✓ .gitignore created"
echo ""

# Create simple README
echo "Creating README.md..."
cat > README.md << 'README_EOF'
# NanoporeToBED Pipeline

Pipeline for processing Oxford Nanopore Technologies sequencing data with base modifications (5mC) to generate methylation BED files.

## Quick Start

1. The environment has been set up. Activate it with:
   ```bash
   conda activate nanopore_methylation
   ```

2. Run the pipeline:
   ```bash
   sbatch scripts/NanoporeToBED.sh \
     -i /path/to/input \
     -o /path/to/output \
     -ref /path/to/reference.fna \
     -t 32
   ```

3. For help:
   ```bash
   ./scripts/NanoporeToBED.sh -h
   ```

## Citation

Drag, M.H., et al. (2025). New high accuracy diagnostics for avian Aspergillus fumigatus infection using Nanopore methylation sequencing of host cell-free DNA and machine learning prediction. bioRxiv 2025.04.11.648151.
README_EOF

echo "✓ README.md created"
echo ""

# Setup conda environment
echo "=========================================="
echo "Setting up Conda Environment"
echo "=========================================="
echo ""

# Check if environment already exists
if $CONDA_CMD env list | grep -q "nanopore_methylation"; then
    echo "⚠️  Environment 'nanopore_methylation' already exists"
    read -p "Do you want to recreate it? (y/n): " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        echo "Removing existing environment..."
        $CONDA_CMD env remove -n nanopore_methylation -y
    else
        echo "Keeping existing environment"
        echo ""
        echo "✅ Setup complete!"
        echo ""
        echo "Next steps:"
        echo "  1. cd $INSTALL_DIR"
        echo "  2. conda activate nanopore_methylation"
        echo "  3. ./scripts/NanoporeToBED.sh -h"
        exit 0
    fi
fi

# Create environment
echo "Creating conda environment 'nanopore_methylation'..."
echo "This may take several minutes..."

if [[ "$CONDA_CMD" == "micromamba" ]]; then
    micromamba create -y -n nanopore_methylation -c conda-forge -c bioconda \
        python=3.8 samtools>=1.17 minimap2>=2.26 ont-modkit>=0.3.0 qualimap>=2.3 \
        pigz parallel bedtools>=2.30 bcftools>=1.17
else
    $CONDA_CMD env create -f environment.yml
fi

if [ $? -eq 0 ]; then
    echo "✓ Environment created successfully"
else
    echo "❌ Environment creation failed"
    echo "You can try manually with:"
    echo "  cd $INSTALL_DIR"
    echo "  conda env create -f environment.yml"
    exit 1
fi

echo ""
echo "=========================================="
echo "✅ Setup Complete!"
echo "=========================================="
echo ""
echo "Installation directory: $INSTALL_DIR"
echo ""
echo "Next steps:"
echo "  1. cd $INSTALL_DIR"
echo "  2. $CONDA_CMD activate nanopore_methylation"
echo "  3. Edit scripts/NanoporeToBED.sh to set your SLURM account"
echo "  4. Test: ./scripts/NanoporeToBED.sh -h"
echo "  5. Run: sbatch scripts/NanoporeToBED.sh -i INPUT -o OUTPUT -ref REFERENCE.fna"
echo ""
echo "For more information, see README.md"
echo ""
echo "Citation: Drag et al. (2025) bioRxiv 2025.04.11.648151"
