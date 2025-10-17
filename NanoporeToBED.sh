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

# Log everything to file as well
exec > >(tee -a "$log_file") 2>&1

echo "Checking reference genome..."
if [[ ! -f "$ref_genome" ]]; then
  echo "❌ ERROR: Reference genome not found: $ref_genome"
  exit 1
fi
ref_size=$(du -h "$ref_genome" | cut -f1)
echo "  ✓ Reference found: $ref_genome ($ref_size)"
echo ""

echo "Scanning for sample directories..."
sample_dirs=""
for srr_dir in "$input_dir"/SRR*/fastq_gpu_hac_only_5mc_mod/pass/*-*_*_*_*/; do
  if [[ -d "$srr_dir" && ! "$srr_dir" =~ unclassified ]]; then
    sample_dirs="$sample_dirs$srr_dir"$'\n'
  fi
done

# Also check for alternative directory structures
if [[ -z "$sample_dirs" ]]; then
  echo "  Checking alternative structure (direct pass folders)..."
  for srr_dir in "$input_dir"/SRR*/pass/*_*/; do
    if [[ -d "$srr_dir" && ! "$srr_dir" =~ unclassified ]]; then
      sample_dirs="$sample_dirs$srr_dir"$'\n'
    fi
  done
fi

sample_dirs=$(echo "$sample_dirs" | sort)

if [[ -z "$sample_dirs" ]]; then
  echo "❌ ERROR: No sample directories found matching pattern"
  echo ""
  echo "Searched in: $input_dir"
  echo "Looking for patterns:"
  echo "  1. SRR*/fastq_gpu_hac_only_5mc_mod/pass/*_*/"
  echo "  2. SRR*/pass/*_*/"
  echo ""
  echo "Please check:"
  echo "  1. Input directory is correct"
  echo "  2. Sample directories exist in pass/ folders"
  echo "  3. Sample names contain underscores for metadata separation"
  exit 1
fi

sample_count=$(echo "$sample_dirs" | wc -l)
echo "✓ Found $sample_count sample(s) to process"
echo ""

# Group by SRR code for summary
declare -A srr_counts
for sample_path in $sample_dirs; do
  rel_path=$(realpath --relative-to="$input_dir" "$sample_path")
  srr_code=$(echo "$rel_path" | cut -d'/' -f1)
  srr_counts[$srr_code]=$((${srr_counts[$srr_code]:-0} + 1))
done

echo "Sample distribution by library:"
for srr in "${!srr_counts[@]}"; do
  echo "  $srr: ${srr_counts[$srr]} sample(s)"
done | sort
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
  sample_name=$(basename "$sample_path")
  rel_path=$(realpath --relative-to="$input_dir" "$sample_path")
  srr_code=$(echo "$rel_path" | cut -d'/' -f1)

  # Try to extract date from path (if present)
  if [[ "$rel_path" =~ ([0-9]{8}|[0-9]{4}-[0-9]{2}-[0-9]{2}) ]]; then
    date_code="${BASH_REMATCH[1]}"
  else
    date_code=$(echo "$rel_path" | cut -d'/' -f2)
  fi

  echo "=========================================="
  echo "Sample $sample_num of $sample_count"
  echo "=========================================="
  echo "Sample ID:     $sample_name"
  echo "Library (SRR): $srr_code"
  echo "Date:          $date_code"
  echo "Input path:    $sample_path"
  echo "Threads:       $THREADS"
  echo ""

  # Create output directory
  out_sample_dir="$output_dir/$srr_code/$date_code/$sample_name"
  mkdir -p "$out_sample_dir"
  mkdir -p "$output_dir/logs/$srr_code/$date_code"

  sample_log="$output_dir/logs/$srr_code/$date_code/${sample_name}.log"

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

  if [[ -s "$merged_bam" && $(stat -c%s "$merged_bam") -gt 100000000 ]]; then
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

  if [[ -s "$minimap_bam" && $(stat -c%s "$minimap_bam") -gt 100000000 ]]; then
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

  if [[ -s "$methyl_bed" && $(stat -c%s "$methyl_bed") -gt 100000 ]]; then
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
echo "Summary by library:"
for srr in "${!srr_counts[@]}"; do
  echo "  $srr: ${srr_counts[$srr]} sample(s)"
done | sort
echo ""
echo "Citation: Drag et al. (2025) bioRxiv 2025.04.11.648151"
echo "🎉 All done!"
