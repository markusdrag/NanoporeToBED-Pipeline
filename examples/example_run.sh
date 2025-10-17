#!/bin/bash
# Example run configuration for Nanopore methylation analysis

# Set paths
INPUT_DIR="/data/nanopore/fastq_gpu_hac_mod"
OUTPUT_DIR="/data/nanopore/methylation_results"
REFERENCE_GENOME="/data/references/genome.fna"

# Submit job with custom thread count
sbatch scripts/NanoporeToBED.sh \
  -i "$INPUT_DIR" \
  -o "$OUTPUT_DIR" \
  -ref "$REFERENCE_GENOME" \
  -t 32

# Or use default 40 threads
# sbatch scripts/NanoporeToBED.sh \
#   -i "$INPUT_DIR" \
#   -o "$OUTPUT_DIR" \
#   -ref "$REFERENCE_GENOME"
