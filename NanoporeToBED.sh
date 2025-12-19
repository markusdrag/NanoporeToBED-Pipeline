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
micromamba activate nanopore_methylation 2>/dev/null || \
  conda activate nanopore_methylation 2>/dev/null || \
  echo "Warning: Could not activate nanopore_methylation environment"

set -euo pipefail

# Portable file size function (works on Linux and macOS)
file_size_bytes() {
  if [[ "$(uname)" == "Darwin" ]]; then
    stat -f%z "$1" 2>/dev/null || echo 0
  else
    stat -c%s "$1" 2>/dev/null || echo 0
  fi
}

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
      echo ""
      echo "Example:"
      echo "  $0 -i /data/nanopore -o /results -ref /ref/genome.fna -t 32"
      exit 0
      ;;
    *) echo "Unknown option: $1"; echo "Use -h for help"; exit 1 ;;
  esac
done

# Validate required arguments
if [[ -z "${input_dir:-}" || -z "${output_dir:-}" || -z "${ref_genome:-}" ]]; then
  echo "Error: Missing required arguments"
  echo "Usage: $0 -i <input_dir> -o <output_dir> -ref <reference_genome.fna> [-t <threads>] [--dry-run]"
  echo "Use -h for help"
  exit 1
fi

# Validate thread count
if ! [[ "$THREADS" =~ ^[0-9]+$ ]] || [[ "$THREADS" -lt 1 ]]; then
  echo "Error: Thread count must be a positive integer"
  exit 1
fi

# Adjust SLURM allocation if threads differ from default
SLURM_CPUS=${SLURM_CPUS_PER_TASK:-40}
if [[ "$THREADS" -gt "$SLURM_CPUS" ]]; then
  echo "Warning: Requested threads ($THREADS) exceeds SLURM allocation ($SLURM_CPUS)"
  echo "Reducing to $SLURM_CPUS threads"
  THREADS=$SLURM_CPUS
fi

timestamp=$(date +%Y%m%d_%H%M%S)

echo "=========================================="
echo "NanoporeToBED Pipeline"
echo "=========================================="
echo "Started at: $(date)"
echo ""
echo "Configuration:"
echo "  Input directory:    $input_dir"
echo "  Output directory:   $output_dir"
echo "  Reference genome:   $ref_genome"
echo "  Threads:            $THREADS"
echo "  Dry run mode:       $dry_run"
echo ""
echo "Citation: Drag et al. (2025) bioRxiv 2025.04.11.648151"
echo "          https://doi.org/10.1101/2025.04.11.648151"
echo ""

# Convert to absolute paths
input_dir=$(realpath "$input_dir")
output_dir=$(realpath "$output_dir")
ref_genome=$(realpath "$ref_genome")

echo "Absolute paths:"
echo "  Input:  $input_dir"
echo "  Output: $output_dir"
echo "  Ref:    $ref_genome"
echo ""

# Create output directory structure
echo "Creating output directory structure..."
mkdir -p "$output_dir/logs"
log_file="$output_dir/logs/pipeline_master_log_${timestamp}.txt"
echo "  ✓ Created: $output_dir/logs"
echo "  ✓ Master log: $log_file"
echo ""

# Simple logging - append output to log file
# Note: For full logging, run as: bash script.sh 2>&1 | tee logfile.txt

echo "Checking reference genome..."
if [[ ! -f "$ref_genome" ]]; then
  echo "❌ ERROR: Reference genome not found: $ref_genome"
  exit 1
fi
ref_size=$(du -h "$ref_genome" | cut -f1)
echo "  ✓ Reference found: $ref_genome ($ref_size)"
echo ""

echo "Scanning for sample directories..."

# Function to extract barcode number from BAM file
get_barcode() {
  local bam="$1"
  local bc_full=$(samtools view "$bam" 2>/dev/null | head -1 | grep -oE 'barcode[0-9]+' || echo "")
  if [[ -n "$bc_full" ]]; then
    # Extract just the number and format as b## (e.g., barcode03 -> b03)
    local bc_num=$(echo "$bc_full" | grep -oE '[0-9]+')
    printf "b%02d" "$bc_num"
  else
    echo "b00"
  fi
}

# Find all subdirectories containing BAM files
sample_dirs=""
for sample_dir in "$input_dir"/*/; do
  if [[ -d "$sample_dir" ]]; then
    # Check if directory contains BAM files
    bam_count=$(find "$sample_dir" -maxdepth 1 -name '*.bam' 2>/dev/null | wc -l)
    if [[ $bam_count -gt 0 ]]; then
      # Skip unclassified directories
      dir_name=$(basename "$sample_dir")
      if [[ ! "$dir_name" =~ unclassified ]]; then
        sample_dirs="$sample_dirs$sample_dir"$'\n'
      fi
    fi
  fi
done

sample_dirs=$(echo "$sample_dirs" | grep -v '^$' | sort)

if [[ -z "$sample_dirs" ]]; then
  echo "❌ ERROR: No sample directories found containing BAM files"
  echo ""
  echo "Searched in: $input_dir"
  echo "Looking for: subdirectories containing *.bam files"
  echo ""
  echo "Please check:"
  echo "  1. Input directory is correct"
  echo "  2. Sample directories contain BAM files directly"
  exit 1
fi

sample_count=$(echo "$sample_dirs" | wc -l)
echo "✓ Found $sample_count sample(s) to process"
echo ""

if [[ "$dry_run" == true ]]; then
  echo "=========================================="
  echo "DRY RUN MODE - No processing will occur"
  echo "=========================================="
  echo ""
fi

# Process each sample
sample_num=0
for sample_path in $sample_dirs; do
  sample_num=$((sample_num + 1))

  # Extract sample information
  input_dir_name=$(basename "$sample_path")
  
  # Get first BAM file to extract barcode
  first_bam=$(find "$sample_path" -maxdepth 1 -name '*.bam' 2>/dev/null | head -1)
  if [[ -n "$first_bam" ]]; then
    barcode=$(get_barcode "$first_bam")
  else
    barcode="b00"
  fi
  
  # Create sample name with barcode suffix (e.g., L07_HUMB_LAB_b03)
  sample_name="${input_dir_name}_${barcode}"

  echo "=========================================="
  echo "Sample $sample_num of $sample_count"
  echo "=========================================="
  echo "Input dir:     $input_dir_name"
  echo "Barcode:       $barcode"
  echo "Sample ID:     $sample_name"
  echo "Input path:    $sample_path"
  echo "Threads:       $THREADS"
  echo ""

  # Create output directory (flat structure)
  out_sample_dir="$output_dir/$sample_name"
  mkdir -p "$out_sample_dir"
  mkdir -p "$output_dir/logs"

  sample_log="$output_dir/logs/${sample_name}.log"

  echo "Output dir:    $out_sample_dir"
  echo "Sample log:    $sample_log"
  echo ""

  if [[ "$dry_run" == true ]]; then
    echo "[DRY RUN] Would process this sample with $THREADS threads"
    echo ""
    continue
  fi

  # Step 1: Merge BAM files (these contain methylation tags from basecalling)
  echo "Step 1/4: Merging BAM files with methylation tags"
  merged_bam="$out_sample_dir/${sample_name}.merged.bam"

  if [[ -s "$merged_bam" && $(file_size_bytes "$merged_bam") -gt 100000000 ]]; then
    merged_size=$(du -h "$merged_bam" | cut -f1)
    echo "  ✓ Already exists: $merged_bam ($merged_size)"
  else
    echo "  Searching for BAM files..."
    temp_bam_list="$out_sample_dir/bam_list.txt"
    find "$sample_path" -name '*.bam' > "$temp_bam_list"

    bam_count=$(wc -l < "$temp_bam_list")
    echo "    Found: $bam_count BAM files"

    if [[ $bam_count -eq 0 ]]; then
      echo "  ⚠️  WARNING: No BAM files found - skipping this sample"
      echo ""
      continue
    fi

    echo "  Merging BAM files (using $THREADS threads)..."
    {
      read firstbam
      if ! samtools view -@${THREADS} -H "$firstbam" > /dev/null 2>>"$sample_log"; then
        echo "  ❌ BAM header problem in $firstbam" | tee -a "$sample_log"
        continue
      fi
      samtools view -@${THREADS} -h "$firstbam"
      while read bam; do
        if samtools view -@${THREADS} "$bam" > /dev/null 2>>"$sample_log"; then
          samtools view -@${THREADS} "$bam"
        else
          echo "  ⚠️  Corrupt BAM: $bam" | tee -a "$sample_log"
        fi
      done
    } < "$temp_bam_list" | samtools view -@${THREADS} -ubS - | samtools sort -@${THREADS} -o "$merged_bam" -

    samtools index -@${THREADS} "$merged_bam"
    merged_size=$(du -h "$merged_bam" | cut -f1)
    echo "  ✓ BAMs merged: $merged_bam ($merged_size)"
  fi
  echo ""

  # Step 2: Alignment with minimap2 (preserving methylation tags)
  echo "Step 2/4: Aligning reads with minimap2"
  minimap_bam="$out_sample_dir/${sample_name}.minimap.bam"

  if [[ -s "$minimap_bam" && $(file_size_bytes "$minimap_bam") -gt 100000000 ]]; then
    bam_size=$(du -h "$minimap_bam" | cut -f1)
    echo "  ✓ Already exists: $minimap_bam ($bam_size)"
  else
    echo "  Running minimap2 alignment (preserving methylation tags)..."
    echo "    Threads: $THREADS"
    echo "    Mode: map-ont with -y (copy tags)"
    start_time=$(date +%s)

    samtools fastq -@${THREADS} -T MM,ML "$merged_bam" | \
      minimap2 -ax map-ont -t ${THREADS} -y --secondary=no "$ref_genome" - 2>>"$sample_log" | \
      samtools view -@${THREADS} -S -b - 2>>"$sample_log" | \
      samtools sort -@${THREADS} -o "$minimap_bam" -T "$out_sample_dir/reads.tmp" - 2>>"$sample_log"

    echo "  Creating BAM index..."
    samtools index -@${THREADS} "$minimap_bam" 2>>"$sample_log"

    end_time=$(date +%s)
    elapsed=$((end_time - start_time))
    bam_size=$(du -h "$minimap_bam" | cut -f1)

    echo "  ✓ Alignment complete ($elapsed seconds)"
    echo "    Output: $minimap_bam ($bam_size)"
  fi
  echo ""

  # Step 3: Methylation calling
  echo "Step 3/4: Calling methylation with modkit"
  methyl_bed="$out_sample_dir/${sample_name}.CpG.bed"

  if [[ -s "$methyl_bed" && $(file_size_bytes "$methyl_bed") -gt 100000 ]]; then
    bed_size=$(du -h "$methyl_bed" | cut -f1)
    echo "  ✓ Already exists: $methyl_bed ($bed_size)"
  else
    echo "  Running modkit pileup..."
    echo "    Mode: CpG methylation"
    echo "    Threads: $THREADS"
    start_time=$(date +%s)

    modkit pileup "$minimap_bam" "$methyl_bed" \
      --cpg --ref "$ref_genome" -t ${THREADS} --combine-mods 2>>"$sample_log"

    end_time=$(date +%s)
    elapsed=$((end_time - start_time))
    bed_size=$(du -h "$methyl_bed" | cut -f1)

    echo "  ✓ Methylation calling complete ($elapsed seconds)"
    echo "    Output: $methyl_bed ($bed_size)"
  fi
  echo ""

  # Step 4: Quality control
  echo "Step 4/4: Running Qualimap QC"
  qualimap_dir="$out_sample_dir/qualimap"

  if [[ -f "$qualimap_dir/qualimapReport.html" ]]; then
    echo "  ✓ Already exists: $qualimap_dir/qualimapReport.html"
  else
    echo "  Running Qualimap bamqc..."
    echo "    Window size: 5000"
    echo "    Threads: $THREADS"
    start_time=$(date +%s)

    # Note: Qualimap has a maximum thread limit, usually 32
    qualimap_threads=$THREADS
    if [[ $qualimap_threads -gt 32 ]]; then
      qualimap_threads=32
      echo "    (Qualimap limited to 32 threads)"
    fi

    qualimap bamqc -bam "$minimap_bam" -nw 5000 -nt ${qualimap_threads} \
      -c -outdir "$qualimap_dir" &>>"$sample_log"

    end_time=$(date +%s)
    elapsed=$((end_time - start_time))

    echo "  ✓ QC complete ($elapsed seconds)"
    echo "    Report: $qualimap_dir/qualimapReport.html"
  fi
  echo ""

  echo "✅ Sample $sample_name complete!"
  echo ""

done

echo "=========================================="
echo "Pipeline Complete!"
echo "=========================================="
echo "Finished at: $(date)"
echo ""
echo "Processed: $sample_count sample(s)"
echo "Threads used: $THREADS"
echo "Output directory: $output_dir"
echo "Master log: $log_file"
echo ""
echo "Citation: Drag et al. (2025) bioRxiv 2025.04.11.648151"
echo "🎉 All done!"
