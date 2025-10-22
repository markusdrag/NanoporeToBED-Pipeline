#!/bin/bash

# NanoporeToBED Pipeline - Unified Intelligent Setup
# Automatically detects platform and chooses optimal installation strategy
# Supports: Linux x86_64, macOS Intel (x86_64), macOS ARM64 (Apple Silicon)

set -e

echo "=========================================="
echo "NanoporeToBED - Intelligent Setup"
echo "=========================================="
echo "Version: 2.0 (Unified)"
echo "Pipeline: Drag et al. (2025) bioRxiv"
echo ""

# ============================================================================
# PLATFORM DETECTION
# ============================================================================

PLATFORM=$(uname -s)
ARCH=$(uname -m)
echo "[INFO] Detected Platform: $PLATFORM ($ARCH)"
echo ""

# Determine platform category
if [[ "$PLATFORM" == "Darwin" ]] && [[ "$ARCH" == "arm64" ]]; then
    PLATFORM_TYPE="macos_arm64"
    echo "[INFO] System: macOS with Apple Silicon (M1/M2/M3/M4)"
    NEEDS_SPECIAL_HANDLING=true
elif [[ "$PLATFORM" == "Darwin" ]] && [[ "$ARCH" == "x86_64" ]]; then
    PLATFORM_TYPE="macos_intel"
    echo "[INFO] System: macOS with Intel processor"
    NEEDS_SPECIAL_HANDLING=false
elif [[ "$PLATFORM" == "Linux" ]]; then
    PLATFORM_TYPE="linux"
    echo "[INFO] System: Linux"
    NEEDS_SPECIAL_HANDLING=false
else
    echo "[WARNING] Unsupported platform: $PLATFORM ($ARCH)"
    echo "Attempting standard installation..."
    PLATFORM_TYPE="unknown"
    NEEDS_SPECIAL_HANDLING=false
fi
echo ""

# ============================================================================
# CHECK FOR HPC ENVIRONMENT
# ============================================================================

if command -v module &> /dev/null; then
    echo "[INFO] HPC environment detected (module system available)"
    echo "Attempting to load conda module..."
    module load conda 2>/dev/null || module load miniconda3 2>/dev/null || module load anaconda3 2>/dev/null || true
    echo ""
fi

# ============================================================================
# DETECT CONDA VARIANT
# ============================================================================

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

# ============================================================================
# HELPER FUNCTION FOR ENVIRONMENT ACTIVATION
# ============================================================================

activate_env() {
    local env_name=$1

    if [[ "$CONDA_CMD" == "micromamba" ]]; then
        eval "$(micromamba shell hook --shell bash)"
        micromamba activate $env_name
    elif [[ "$CONDA_CMD" == "mamba" ]]; then
        eval "$(mamba shell hook --shell bash)"
        mamba activate $env_name
    else
        eval "$(conda shell.bash hook)"
        conda activate $env_name
    fi
}

deactivate_env() {
    # After proper activation via shell hook, conda deactivate works for all
    # The issue is if conda isn't in PATH, we use the tool-specific deactivate

    if [[ "$CONDA_CMD" == "micromamba" ]]; then
        micromamba deactivate 2>/dev/null || conda deactivate 2>/dev/null || true
    elif [[ "$CONDA_CMD" == "mamba" ]]; then
        conda deactivate 2>/dev/null || mamba deactivate 2>/dev/null || true
    else
        conda deactivate 2>/dev/null || true
    fi
}

# ============================================================================
# CHECK IF CONDA IS INSTALLED
# ============================================================================

if [[ "$CONDA_CMD" == "none" ]]; then
    echo "[ERROR] No conda/mamba/micromamba found"
    echo ""
    echo "=========================================="
    echo "Installation Options:"
    echo "=========================================="
    echo ""
    echo "1. Install micromamba (recommended - fastest, standalone):"
    echo "   curl -L https://micro.mamba.pm/install.sh | bash"
    echo ""
    echo "2. Install mamba in existing conda:"
    echo "   conda install -n base -c conda-forge mamba"
    echo ""
    echo "3. Install miniconda:"
    if [[ "$PLATFORM" == "Darwin" ]]; then
        echo "   curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-$ARCH.sh"
        echo "   bash Miniconda3-latest-MacOSX-$ARCH.sh"
    else
        echo "   wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh"
        echo "   bash Miniconda3-latest-Linux-x86_64.sh"
    fi
    exit 1
fi

echo "[OK] Using conda variant: $CONDA_CMD"

if [[ "$CONDA_CMD" == "mamba" ]] || [[ "$CONDA_CMD" == "micromamba" ]]; then
    echo "[INFO] Fast solver detected - installation will be quick!"
else
    echo "[WARNING] Using standard conda (slower)"
    echo "         Consider installing mamba for faster setup:"
    echo "         conda install -n base -c conda-forge mamba"
fi
echo ""

# ============================================================================
# ENVIRONMENT CONFIGURATION
# ============================================================================

ENV_NAME="nanopore_methylation"

# Check if environment already exists
if $CONDA_CMD env list 2>/dev/null | grep -q "^$ENV_NAME "; then
    echo "[WARNING] Environment '$ENV_NAME' already exists"
    echo ""
    read -p "Choose action: [K]eep existing, [R]emove and recreate, [A]bort: " -n 1 -r
    echo ""
    case $REPLY in
        [Kk])
            echo "[INFO] Keeping existing environment"
            echo ""
            echo "=========================================="
            echo "[SUCCESS] Using existing environment"
            echo "=========================================="
            echo ""
            echo "Activate with: $CONDA_CMD activate $ENV_NAME"
            exit 0
            ;;
        [Rr])
            echo "[INFO] Removing existing environment..."
            $CONDA_CMD env remove -n $ENV_NAME -y 2>/dev/null || true
            echo ""
            ;;
        *)
            echo "[INFO] Installation aborted by user"
            exit 0
            ;;
    esac
fi

# ============================================================================
# MACOS ARM64 INSTALLATION STRATEGY
# ============================================================================

if [[ "$PLATFORM_TYPE" == "macos_arm64" ]]; then
    echo "=========================================="
    echo "macOS ARM64 Installation Strategy"
    echo "=========================================="
    echo ""
    echo "Many bioinformatics tools lack native ARM64 builds."
    echo "Choose your preferred installation method:"
    echo ""
    echo "1. Native ARM64 + pip fallbacks    [RECOMMENDED - Fast]"
    echo "   - Uses native ARM64 packages where available"
    echo "   - Falls back to pip for missing tools (like modkit)"
    echo "   - Fastest performance, minimal compatibility issues"
    echo ""
    echo "2. x86_64 via Rosetta 2            [Full Compatibility]"
    echo "   - Uses Intel packages translated by Rosetta 2"
    echo "   - All tools available, slightly slower"
    echo "   - Requires Rosetta 2 (auto-installed by macOS)"
    echo ""
    echo "3. Minimal conda + Docker          [Isolation]"
    echo "   - Python environment only"
    echo "   - Heavy tools via Docker containers"
    echo "   - Requires Docker Desktop"
    echo ""
    read -p "Choose method (1/2/3) [default: 1]: " -n 1 -r
    echo ""

    if [[ -z "$REPLY" ]]; then
        INSTALL_METHOD=1
    else
        INSTALL_METHOD=$REPLY
    fi
    echo ""
fi

# ============================================================================
# INSTALLATION FUNCTIONS
# ============================================================================

install_standard() {
    echo "=========================================="
    echo "Standard Installation"
    echo "=========================================="
    echo ""
    echo "Creating environment: $ENV_NAME"
    echo "This may take 2-5 minutes with mamba/micromamba..."
    echo "or 10-20 minutes with standard conda"
    echo ""

    if [[ "$CONDA_CMD" == "mamba" ]] || [[ "$CONDA_CMD" == "micromamba" ]]; then
        $CONDA_CMD create -y -n $ENV_NAME \
            -c conda-forge -c bioconda -c defaults \
            python=3.9 \
            samtools \
            minimap2 \
            ont-modkit \
            qualimap \
            pigz \
            parallel \
            bedtools \
            bcftools \
            --no-channel-priority
    else
        conda create -y -n $ENV_NAME python=3.9

        activate_env $ENV_NAME


        echo "[INFO] Installing core tools..."
        conda install -y -c conda-forge -c bioconda \
            samtools minimap2 pigz parallel

        echo "[INFO] Installing bioinformatics tools..."
        conda install -y -c bioconda \
            ont-modkit bedtools bcftools qualimap

        deactivate_env
    fi
}

install_macos_arm64_native() {
    echo "=========================================="
    echo "macOS ARM64 - Native Installation"
    echo "=========================================="
    echo ""

    if [[ "$CONDA_CMD" == "mamba" ]] || [[ "$CONDA_CMD" == "micromamba" ]]; then
        echo "[INFO] Creating environment with ARM64-compatible packages..."

        $CONDA_CMD create -y -n $ENV_NAME \
            -c conda-forge -c bioconda \
            python=3.9 \
            samtools \
            minimap2 \
            bedtools \
            bcftools \
            pigz \
            parallel \
            --no-channel-priority

        echo ""
        echo "[INFO] Core environment created"
        echo "[INFO] Installing additional tools via pip..."

        # Activate environment
        activate_env $ENV_NAME

        # Install Python packages
        pip install --quiet pysam pandas numpy scipy matplotlib seaborn

        # Install modkit - try conda first with proper channel specification
        echo "[INFO] Installing modkit..."
        $CONDA_CMD install -y bioconda::ont-modkit 2>/dev/null || \
        $CONDA_CMD install -y -c bioconda ont-modkit 2>/dev/null || \
        pip install modkit || \
        echo "[WARNING] modkit installation failed - install manually with: conda install bioconda::ont-modkit"

        # Try qualimap from conda (may fail on ARM64)
        echo ""
        echo "[INFO] Attempting to install qualimap..."
        $CONDA_CMD install -y -c bioconda qualimap 2>/dev/null || {
            echo "[WARNING] Qualimap not available for ARM64"
            echo "          Install separately if needed:"
            echo "          wget https://bitbucket.org/kokonech/qualimap/downloads/qualimap_v2.3.zip"
        }

        deactivate_env
    else
        conda create -y -n $ENV_NAME python=3.9

        activate_env $ENV_NAME


        echo "[INFO] Installing core tools..."
        conda install -y -c conda-forge -c bioconda \
            samtools minimap2 bedtools bcftools pigz parallel

        echo "[INFO] Installing Python packages..."
        pip install --quiet pysam pandas numpy scipy matplotlib seaborn

        echo "[INFO] Installing modkit..."
        conda install -y bioconda::ont-modkit 2>/dev/null || \
        conda install -y -c bioconda ont-modkit 2>/dev/null || \
        pip install modkit || \
        echo "[WARNING] modkit installation failed - install manually"

        echo "[INFO] Attempting qualimap..."
        conda install -y -c bioconda qualimap 2>/dev/null || \
            echo "[WARNING] Qualimap not available for ARM64"

        deactivate_env
    fi
}

install_macos_arm64_rosetta() {
    echo "=========================================="
    echo "macOS ARM64 - Rosetta 2 Installation"
    echo "=========================================="
    echo ""
    echo "[INFO] Installing x86_64 packages via Rosetta 2"
    echo "[INFO] First run may be slower as Rosetta translates binaries"
    echo ""

    # Force x86_64 architecture
    export CONDA_SUBDIR=osx-64

    if [[ "$CONDA_CMD" == "mamba" ]] || [[ "$CONDA_CMD" == "micromamba" ]]; then
        $CONDA_CMD create -y -n $ENV_NAME \
            -c conda-forge -c bioconda -c defaults \
            python=3.9 \
            samtools \
            minimap2 \
            ont-modkit \
            qualimap \
            pigz \
            parallel \
            bedtools \
            bcftools \
            --no-channel-priority
    else
        conda create -y -n $ENV_NAME python=3.9

        activate_env $ENV_NAME


        # Set architecture permanently for this environment
        conda config --env --set subdir osx-64

        conda install -y -c conda-forge -c bioconda \
            samtools minimap2 pigz parallel

        conda install -y -c bioconda \
            ont-modkit bedtools bcftools qualimap

        deactivate_env
    fi

    unset CONDA_SUBDIR

    echo ""
    echo "[OK] x86_64 environment created (runs via Rosetta 2)"
}

install_minimal_docker() {
    echo "=========================================="
    echo "Minimal Environment + Docker Setup"
    echo "=========================================="
    echo ""

    echo "[INFO] Creating minimal Python environment..."

    $CONDA_CMD create -y -n $ENV_NAME \
        -c conda-forge \
        python=3.9 \
        pip \
        pandas \
        numpy \
        scipy \
        matplotlib \
        seaborn

    # Activate environment properly
    if [[ "$CONDA_CMD" == "micromamba" ]]; then
        eval "$(micromamba shell hook --shell bash)"
        micromamba activate $ENV_NAME
    elif [[ "$CONDA_CMD" == "mamba" ]]; then
        eval "$(mamba shell hook --shell bash)"
        mamba activate $ENV_NAME
    else
        activate_env $ENV_NAME

    fi

    echo "[INFO] Installing Python packages..."
    pip install --quiet pysam pandas numpy scipy matplotlib seaborn

    echo "[INFO] Installing modkit..."
    conda install -y bioconda::ont-modkit 2>/dev/null || \
    pip install modkit || \
    echo "[WARNING] modkit installation failed"

    deactivate_env

    echo ""
    echo "[OK] Minimal environment created"
    echo ""
    echo "=========================================="
    echo "Docker Setup Required"
    echo "=========================================="
    echo ""
    echo "Install Docker Desktop:"
    echo "  https://www.docker.com/products/docker-desktop"
    echo ""
    echo "After installing Docker, pull these containers:"
    echo ""
    echo "  docker pull ontresearch/modkit:latest"
    echo "  docker pull staphb/samtools:latest"
    echo "  docker pull staphb/minimap2:latest"
    echo "  docker pull staphb/bedtools:latest"
    echo ""
    echo "Usage example:"
    echo "  docker run -v \$(pwd):/data ontresearch/modkit:latest \\"
    echo "    modkit pileup /data/input.bam /data/output.bed"
    echo ""
}

# ============================================================================
# EXECUTE INSTALLATION
# ============================================================================

case $PLATFORM_TYPE in
    macos_arm64)
        case $INSTALL_METHOD in
            1)
                install_macos_arm64_native
                ;;
            2)
                install_macos_arm64_rosetta
                ;;
            3)
                install_minimal_docker
                ;;
            *)
                echo "[ERROR] Invalid choice"
                exit 1
                ;;
        esac
        ;;
    macos_intel|linux|unknown)
        install_standard
        ;;
esac

# ============================================================================
# VERIFICATION
# ============================================================================

echo ""
echo "=========================================="
echo "Installation Verification"
echo "=========================================="
echo ""

# Activate environment for testing
if [[ "$CONDA_CMD" == "micromamba" ]]; then
    eval "$(micromamba shell hook --shell bash)"
    micromamba activate $ENV_NAME
elif [[ "$CONDA_CMD" == "mamba" ]]; then
    eval "$(mamba shell hook --shell bash)"
    mamba activate $ENV_NAME
else
    activate_env $ENV_NAME

fi

INSTALLED_TOOLS=()
MISSING_TOOLS=()

# Function to check and record tool status
check_tool() {
    local tool=$1
    local check_cmd=$2

    if command -v $tool &> /dev/null; then
        local version=$($check_cmd 2>&1 | head -n1 || echo "version unknown")
        echo "  [✓] $tool: $version"
        INSTALLED_TOOLS+=("$tool")
        return 0
    else
        echo "  [✗] $tool: NOT FOUND"
        MISSING_TOOLS+=("$tool")
        return 1
    fi
}

echo "Checking installed tools:"
echo ""

check_tool python "python --version"
check_tool samtools "samtools --version"
check_tool minimap2 "minimap2 --version"
check_tool modkit "modkit --version"
check_tool bedtools "bedtools --version"
check_tool bcftools "bcftools --version"
check_tool pigz "pigz --version"
check_tool parallel "parallel --version"

# Special check for qualimap (Java-based)
if command -v qualimap &> /dev/null; then
    echo "  [✓] qualimap: installed"
    INSTALLED_TOOLS+=("qualimap")
else
    echo "  [✗] qualimap: NOT FOUND"
    MISSING_TOOLS+=("qualimap")
fi

deactivate_env

echo ""
echo "Installation Summary: ${#INSTALLED_TOOLS[@]}/9 tools installed"

# ============================================================================
# REPORT MISSING TOOLS
# ============================================================================

if [ ${#MISSING_TOOLS[@]} -gt 0 ]; then
    echo ""
    echo "=========================================="
    echo "Missing Tools - Installation Options"
    echo "=========================================="
    echo ""

    for tool in "${MISSING_TOOLS[@]}"; do
        case $tool in
            modkit)
                echo "→ modkit:"
                echo "  Option 1: pip install modkit"
                echo "  Option 2: conda install -n $ENV_NAME -c bioconda ont-modkit"
                echo "  Option 3: https://github.com/nanoporetech/modkit"
                ;;
            qualimap)
                echo "→ qualimap (optional - for QC reports):"
                echo "  Option 1: conda install -n $ENV_NAME -c bioconda qualimap"
                echo "  Option 2: Download from http://qualimap.conesalab.org/"
                echo "  Note: Requires Java 8+"
                ;;
            samtools|minimap2|bedtools|bcftools)
                echo "→ $tool:"
                echo "  conda install -n $ENV_NAME -c bioconda $tool"
                if [[ "$PLATFORM" == "Darwin" ]]; then
                    echo "  OR: brew install $tool"
                fi
                ;;
            *)
                echo "→ $tool:"
                echo "  conda install -n $ENV_NAME -c bioconda $tool"
                ;;
        esac
        echo ""
    done
fi

# ============================================================================
# CREATE REFERENCE FILE
# ============================================================================

cat > environment_setup_complete.txt << EOF
========================================
NanoporeToBED Environment Setup Complete
========================================

Installation Date: $(date)
Platform: $PLATFORM ($ARCH)
Platform Type: $PLATFORM_TYPE
Environment Name: $ENV_NAME
Conda Tool: $CONDA_CMD

========================================
Quick Start
========================================

Activate environment:
  $CONDA_CMD activate $ENV_NAME

Deactivate when done:
  deactivate_env

========================================
Installed Tools (${#INSTALLED_TOOLS[@]}/9)
========================================

${INSTALLED_TOOLS[*]}

Core Tools:
  • samtools    - BAM/CRAM/SAM file processing
  • minimap2    - Fast long-read aligner
  • modkit      - ONT methylation calling
  • bedtools    - Genomic interval operations
  • bcftools    - VCF/BCF variant file processing

Utilities:
  • pigz        - Parallel gzip compression
  • parallel    - GNU parallel for task distribution
  • qualimap    - Quality control and metrics

Python Packages:
  • pysam, pandas, numpy, scipy, matplotlib, seaborn

EOF

if [ ${#MISSING_TOOLS[@]} -gt 0 ]; then
    cat >> environment_setup_complete.txt << EOF
========================================
Missing Tools (${#MISSING_TOOLS[@]})
========================================

${MISSING_TOOLS[*]}

See installation options above.

EOF
fi

if [[ "$PLATFORM_TYPE" == "macos_arm64" ]]; then
    cat >> environment_setup_complete.txt << EOF
========================================
Platform Notes - macOS ARM64
========================================

EOF
    case $INSTALL_METHOD in
        1)
            cat >> environment_setup_complete.txt << EOF
Installation Method: Native ARM64 + pip fallbacks

This installation uses native ARM64 packages where
available and pip for tools without ARM64 builds.

Performance: Excellent (native execution)
Compatibility: Very Good (pip fills gaps)
EOF
            ;;
        2)
            cat >> environment_setup_complete.txt << EOF
Installation Method: x86_64 via Rosetta 2

This installation uses Intel (x86_64) packages
translated by Rosetta 2 on the fly.

Performance: Good (slight translation overhead)
Compatibility: Excellent (all tools available)
Note: First run of tools may be slower
EOF
            ;;
        3)
            cat >> environment_setup_complete.txt << EOF
Installation Method: Minimal + Docker

This installation provides a minimal Python
environment. Heavy bioinformatics tools should
be run via Docker containers.

Setup Docker containers:
  docker pull ontresearch/modkit:latest
  docker pull staphb/samtools:latest
  docker pull staphb/minimap2:latest

Usage:
  docker run -v \$(pwd):/data ontresearch/modkit:latest \\
    modkit pileup /data/input.bam /data/output.bed
EOF
            ;;
    esac
    echo "" >> environment_setup_complete.txt
fi

cat >> environment_setup_complete.txt << EOF
========================================
Citation
========================================

If using this pipeline, please cite:
  Drag et al. (2025) bioRxiv 2025.04.11.648151

Tool Citations:
  • modkit    - Oxford Nanopore Technologies
  • samtools  - Li et al. (2009) Bioinformatics 25:2078
  • minimap2  - Li (2018) Bioinformatics 34:3094
  • bedtools  - Quinlan & Hall (2010) Bioinformatics 26:841

========================================
Troubleshooting
========================================

Check environment:
  $CONDA_CMD env list

List installed packages:
  conda list -n $ENV_NAME

Reinstall tool:
  conda install -n $ENV_NAME -c bioconda <tool>

For more help:
  https://github.com/nanoporetech/modkit
  https://bioconda.github.io/

========================================
EOF

echo "=========================================="
echo "[SUCCESS] Setup Complete!"
echo "=========================================="
echo ""
echo "Environment: $ENV_NAME"
echo "Tools installed: ${#INSTALLED_TOOLS[@]}/9"
echo ""
echo "=========================================="
echo "Next Steps"
echo "=========================================="
echo ""
echo "1. Activate the environment:"
echo "   $CONDA_CMD activate $ENV_NAME"
echo ""
echo "2. Verify tools work:"
echo "   samtools --version"
echo "   minimap2 --version"
echo "   modkit --version"
echo ""
echo "3. When done, deactivate:"
echo "   conda deactivate"
echo ""
echo "=========================================="
echo ""
echo "Setup details saved to: environment_setup_complete.txt"
echo ""

if [[ "$PLATFORM_TYPE" == "macos_arm64" ]] && [[ "$INSTALL_METHOD" == "3" ]]; then
    echo "=========================================="
    echo "IMPORTANT: Docker Setup Required"
    echo "=========================================="
    echo ""
    echo "Don't forget to install Docker Desktop and"
    echo "pull the required containers (see above)"
    echo ""
fi

# Final platform-specific tips
if [[ "$PLATFORM_TYPE" == "macos_arm64" ]]; then
    echo "=========================================="
    echo "macOS ARM64 Tips"
    echo "=========================================="
    echo ""
    if [[ "$INSTALL_METHOD" == "1" ]]; then
        echo "✓ Using native ARM64 for best performance"
        echo "✓ pip installed modkit and Python packages"
        if [ ${#MISSING_TOOLS[@]} -gt 0 ]; then
            echo "⚠ Some tools missing - see options above"
        fi
    elif [[ "$INSTALL_METHOD" == "2" ]]; then
        echo "✓ Running via Rosetta 2 for full compatibility"
        echo "ℹ First tool execution may be slightly slower"
    else
        echo "✓ Minimal environment ready"
        echo "⚠ Remember to set up Docker for analysis tools"
    fi
    echo ""
fi

exit 0
