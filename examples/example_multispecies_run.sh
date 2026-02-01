#!/bin/bash

# Example: Multi-Species Run with NanoporeToBED Pipeline v1.5.0
# Scenario: 11 Minipigs (barcodes 1-11) + 5 Penguins (barcodes 12-16)

# SLURM headers (adjust for your cluster)
#SBATCH --job-name=MultiSpecies_NanoporeToBED
#SBATCH --output=multispecies_run.out
#SBATCH --error=multispecies_run.err
#SBATCH -c 40
#SBATCH --mem 192g
#SBATCH --time=72:00:00
#SBATCH --account YourAccount

# Activate environment
conda activate nanopore_methylation

# Define paths
INPUT_DIR="/path/to/fastq_gpu_hac_mod"
OUTPUT_DIR="/path/to/output/mixed_species"
REF_PIG="/path/to/references/Sus_scrofa.fna"
REF_PENGUIN="/path/to/references/Pygoscelis_adeliae.fna"
THREADS=40

# Run pipeline with multi-species mapping
bash NanoporeToBED.sh \
  -i "$INPUT_DIR" \
  -o "$OUTPUT_DIR" \
  --multi-mapping "1:11,12:16" \
  --multi-refs "$REF_PIG,$REF_PENGUIN" \
  -t "$THREADS" \
  --expanded-plots

# Expected output structure:
#
# output/mixed_species/
# ├── D01_Minipig_Control_b01/
# │   ├── D01_Minipig_Control_b01.CpG.bed  <- aligned to Sus_scrofa.fna
# │   └── ...
# ├── D02_Minipig_Treated_b02/
# │   ├── D02_Minipig_Treated_b02.CpG.bed  <- aligned to Sus_scrofa.fna
# │   └── ...
# ...
# ├── P01_Emperor_Control_b12/
# │   ├── P01_Emperor_Control_b12.CpG.bed  <- aligned to Pygoscelis_adeliae.fna
# │   └── ...
# └── logs/
#     └── pipeline_master_log_*.txt

echo "Multi-species processing complete!"
