#!/bin/bash

# NanoporeToBED Pipeline
# Version: 1.7.0

#SBATCH --job-name=NanoporeToBED
#SBATCH --output=NanoporeToBED.out
#SBATCH --error=NanoporeToBED.err
#SBATCH -c 40
#SBATCH --mem 192g
#SBATCH --time=72:00:00
#SBATCH --account YourAccount

VERSION="1.7.0"

# Resolve the script directory before any later cd, so the summary script is found
# when this file is started with a relative path.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

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

# Multi-species mapping: Parse barcode ranges and build mapping table
# Creates global parallel arrays: BARCODES and REFERENCES
parse_multi_mapping() {
  local mapping="$1"
  local refs="$2"
  
  # Split comma-separated values
  IFS=',' read -ra RANGES <<< "$mapping"
  IFS=',' read -ra REFS <<< "$refs"
  
  # Validate matching counts
  if [[ ${#RANGES[@]} -ne ${#REFS[@]} ]]; then
    echo "Error: Number of barcode ranges (${#RANGES[@]}) doesn't match number of references (${#REFS[@]})"
    echo "  Ranges: $mapping"
    echo "  References: $refs"
    exit 1
  fi
  
  # Process each range and map to corresponding reference
  for i in "${!RANGES[@]}"; do
    range="${RANGES[$i]}"
    ref="${REFS[$i]}"
    
    # Trim whitespace
    range=$(echo "$range" | xargs)
    ref=$(echo "$ref" | xargs)
    
    # Parse range: "1:11" (range) or "5" (single)
    if [[ "$range" =~ ^([0-9]+):([0-9]+)$ ]]; then
      # Range format
      start="${BASH_REMATCH[1]}"
      end="${BASH_REMATCH[2]}"
      
      if [[ $start -gt $end ]]; then
        echo "Error: Invalid range '$range' (start > end)"
        exit 1
      fi
      
      # Map all barcodes in range to this reference
      for bc_num in $(seq $start $end); do
        bc_formatted=$(printf "b%02d" $bc_num)
        BARCODES+=("$bc_formatted")
        REFERENCES+=("$ref")
      done
      
    elif [[ "$range" =~ ^([0-9]+)$ ]]; then
      # Single barcode
      bc_num="${BASH_REMATCH[1]}"
      bc_formatted=$(printf "b%02d" $bc_num)
      BARCODES+=("$bc_formatted")
      REFERENCES+=("$ref")
      
    else
      echo "Error: Invalid range format '$range'"
      echo "  Expected: '1:11' (range) or '5' (single barcode)"
      exit 1
    fi
  done
  
  # Validate all reference files exist
  local checked_refs=()
  for ref in "${REFS[@]}"; do
    ref=$(echo "$ref" | xargs)
    # Skip if already checked
    local already_checked=false
    if [[ ${#checked_refs[@]} -gt 0 ]]; then
      for checked in "${checked_refs[@]}"; do
        if [[ "$checked" == "$ref" ]]; then
          already_checked=true
          break
        fi
      done
    fi
    if [[ "$already_checked" == true ]]; then
      continue
    fi
    checked_refs+=("$ref")
    
    if [[ ! -f "$ref" ]]; then
      echo "Error: Reference genome not found: $ref"
      exit 1
    fi
  done
  
  echo "  [OK] Multi-species mapping parsed successfully"
}

# Get reference genome for a given barcode using parallel arrays
get_reference_for_sample() {
  local barcode="$1"
  
  # Search through parallel arrays
  for i in "${!BARCODES[@]}"; do
    if [[ "${BARCODES[$i]}" == "$barcode" ]]; then
      echo "${REFERENCES[$i]}"
      return 0
    fi
  done
  
  # No match found
  echo "Error: No reference genome mapped for barcode $barcode" >&2
  echo "Available barcode mappings:" >&2
  for i in "${!BARCODES[@]}"; do
    echo "  ${BARCODES[$i]} -> ${REFERENCES[$i]}" >&2
  done | sort >&2
  return 1
}

# Default values
THREADS=40
dry_run=false
include_fail=false
include_barcodes=false
expanded_plots=false
combine_strands=false
force=false
aligned_bams=false
multi_mapping=""
multi_refs=""

# Global arrays for multi-species mapping (parallel arrays: barcodes and refs)
declare -a BARCODES=()
declare -a REFERENCES=()

# Parse arguments
while [[ "$#" -gt 0 ]]; do
  case "$1" in
    -i|--input) input_dir="$2"; shift 2 ;;
    -o|--output) output_dir="$2"; shift 2 ;;
    -ref|--reference) ref_genome="$2"; shift 2 ;;
    -t|--threads) THREADS="$2"; shift 2 ;;
    --multi-mapping) multi_mapping="$2"; shift 2 ;;
    --multi-refs) multi_refs="$2"; shift 2 ;;
    --dry-run) dry_run=true; shift ;;
    --include-fail) include_fail=true; shift ;;
    --include-empty-barcodes) include_barcodes=true; shift ;;
    --expanded-plots) expanded_plots=true; shift ;;
    --combine-strands) combine_strands=true; shift ;;
    --force) force=true; shift ;;
    --aligned-bams) aligned_bams=true; shift ;;
    -h|--help)
      echo "Usage: $0 -i <input_dir> -o <output_dir> [-ref <reference.fna> | --multi-mapping <ranges> --multi-refs <refs>] [options]"
      echo ""
      echo "Required arguments (single-species mode):"
      echo "  -i, --input       Input directory containing pass/fail subdirectories with barcoded samples"
      echo "  -o, --output      Output directory for processed data"
      echo "  -ref, --reference Path to reference genome FASTA file"
      echo ""
      echo "Required arguments (multi-species mode):"
      echo "  -i, --input         Input directory containing pass/fail subdirectories"
      echo "  -o, --output        Output directory for processed data"
      echo "  --multi-mapping     Barcode ranges (e.g., '1:11,12:16' or '1:5,10,15:20')"
      echo "  --multi-refs        Comma-separated reference genomes (e.g., 'pig.fna,penguin.fna')"
      echo ""
      echo "Optional arguments:"
      echo "  -t, --threads              Number of threads to use (default: 40)"
      echo "  --dry-run                  Run in test mode without processing"
      echo "  --include-fail             Also process samples from fail/ directory (default: only pass/)"
      echo "  --include-empty-barcodes   Include barcode## directories (default: skip, only process named samples)"
      echo "  --expanded-plots           Generate extended analysis plots (distribution, QC, comparative)"
      echo "  --combine-strands          Sum both cytosines of each CpG onto the plus-strand coordinate."
      echo "                             Default keeps the plus and minus rows separate."
      echo "                             Sites with few reads on each strand can then pass a later coverage cutoff."
      echo "                             A difference between the two strands is averaged."
      echo "                             On one shallow cfDNA cohort, coverage 10 gained 1,529 CpGs this way,"
      echo "                             and 18 of 35 samples still had no site."
      echo "                             An existing CpG bed is left unchanged unless --force is also set."
      echo "  --force                    Re-run merge, alignment, methylation calling, and Qualimap"
      echo "                             even when those outputs already exist."
      echo "                             Each step replaces its output only after that step succeeds."
      echo "                             With --aligned-bams, only modkit is re-run."
      echo "  --aligned-bams             Run modkit on BAMs that are already aligned to -ref."
      echo "                             -i is one BAM, or a directory of *.bam files."
      echo "                             Files in subdirectories are ignored."
      echo "                             Merge, alignment, and Qualimap are not run."
      echo "                             An unaligned merged BAM cannot be used in this mode."
      echo "  -h, --help                 Show this help message"
      echo ""
      echo "Expected input structure:"
      echo "  input_dir/"
      echo "  ├── pass/              <- processed by default"
      echo "  │   ├── barcode01/     <- skipped by default (use --include-empty-barcodes)"
      echo "  │   ├── D01_SAMPLE/    <- processed (named samples)"
      echo "  │   └── D02_SAMPLE/"
      echo "  └── fail/              <- only with --include-fail"
      echo ""
      echo "Examples:"
      echo "  # Single species"
      echo "  $0 -i /data/nanopore/fastq_gpu_hac_mod -o /results -ref /ref/genome.fna -t 32"
      echo ""
      echo "  # Multiple species (barcodes 1-11 = pigs, 12-16 = penguins)"
      echo "  $0 -i /data/mixed -o /results \\"
      echo "    --multi-mapping '1:11,12:16' \\"
      echo "    --multi-refs '/refs/pig.fna,/refs/penguin.fna' \\"
      echo "    -t 32"
      echo ""
      echo "  # Rebuild every sample, with one row per CpG"
      echo "  $0 -i /data/nanopore/fastq_gpu_hac_mod -o /results -ref /ref/genome.fna --combine-strands --force"
      echo ""
      echo "  # Modkit only, from BAMs that are already aligned"
      echo "  $0 --aligned-bams -i /data/mapped_bam_files -o /results -ref /ref/genome.fna"
      echo ""
      echo "  # Same, with one row per CpG, replacing any existing bed"
      echo "  $0 --aligned-bams -i /data/mapped_bam_files -o /results -ref /ref/genome.fna --combine-strands --force"
      exit 0
      ;;
    *) echo "Unknown option: $1"; echo "Use -h for help"; exit 1 ;;
  esac
done

# Validate required arguments
if [[ -z "${input_dir:-}" || -z "${output_dir:-}" ]]; then
  echo "Error: Missing required arguments"
  echo "Usage: $0 -i <input_dir> -o <output_dir> [-ref <reference.fna> | --multi-mapping <ranges> --multi-refs <refs>]"
  echo "Use -h for help"
  exit 1
fi

# Validate reference mode (single-species OR multi-species)
if [[ -n "${ref_genome:-}" && -n "$multi_mapping" ]]; then
  echo "Error: Cannot use both --reference and --multi-mapping/--multi-refs"
  echo "Use --reference for single-species, or --multi-mapping + --multi-refs for multi-species"
  exit 1
fi

if [[ -z "${ref_genome:-}" && -z "$multi_mapping" ]]; then
  echo "Error: Must provide either --reference or --multi-mapping + --multi-refs"
  echo "Use -h for help"
  exit 1
fi

if [[ -n "$multi_mapping" && -z "$multi_refs" ]]; then
  echo "Error: --multi-mapping requires --multi-refs"
  echo "Use -h for help"
  exit 1
fi

if [[ -z "$multi_mapping" && -n "$multi_refs" ]]; then
  echo "Error: --multi-refs requires --multi-mapping"
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
echo "NanoporeToBED Pipeline v${VERSION}"
echo "=========================================="
echo "Started at: $(date)"
echo ""
echo "Configuration:"
echo "  Input directory:    $input_dir"
echo "  Output directory:   $output_dir"

# Display reference configuration based on mode
if [[ -n "$multi_mapping" ]]; then
  echo "  Mode:               Multi-species"
  echo "  Mapping:            $multi_mapping"
  echo "  References:         $multi_refs"
else
  echo "  Mode:               Single-species"
  echo "  Reference genome:   $ref_genome"
fi

if [[ "$combine_strands" == true ]]; then
  strand_mode="summed onto the plus-strand coordinate"
else
  strand_mode="kept as separate rows"
fi

echo "  Threads:            $THREADS"
echo "  Dry run mode:       $dry_run"
echo "  Strands:            $strand_mode"
if [[ "$aligned_bams" == true ]]; then
  echo "  Pipeline:           aligned BAMs, modkit only"
else
  echo "  Pipeline:           full pipeline"
fi
if [[ "$force" == true && "$aligned_bams" == true ]]; then
  echo "  Existing outputs:   re-run modkit"
elif [[ "$force" == true ]]; then
  echo "  Existing outputs:   re-run merge, alignment, methylation calling, and Qualimap"
else
  echo "  Existing outputs:   reuse when already present"
fi
if [[ "$combine_strands" == true && "$force" != true ]]; then
  echo "  Note: existing CpG beds larger than 100 kb are kept. Use --force to call them again."
fi
echo ""
echo "Citation: Drag et al. (2025) bioRxiv 2025.04.11.648151"
echo "          https://doi.org/10.1101/2025.04.11.648151"
echo ""

# Convert to absolute paths
input_dir=$(realpath "$input_dir")
output_dir=$(realpath "$output_dir")

# Handle single-species reference path conversion
if [[ -n "${ref_genome:-}" ]]; then
  ref_genome=$(realpath "$ref_genome")
fi

# Parse multi-species mapping if provided
if [[ -n "$multi_mapping" ]]; then
  echo "Parsing multi-species configuration..."
  parse_multi_mapping "$multi_mapping" "$multi_refs"
  echo ""
  
  # Display mapping summary
  echo "=========================================="
  echo "Multi-Species Barcode Mapping"
  echo "=========================================="
  
  # Get unique references and group barcodes
  declare -a unique_refs=()
  for ref in "${REFERENCES[@]}"; do
    # Check if already in list
    already_added=false
    if [[ ${#unique_refs[@]} -gt 0 ]]; then
      for uref in "${unique_refs[@]}"; do
        if [[ "$uref" == "$ref" ]]; then
          already_added=true
          break
        fi
      done
    fi
    if [[ "$already_added" == false ]]; then
      unique_refs+=("$ref")
    fi
  done
  
  # Display grouped by reference
  for ref in "${unique_refs[@]}"; do
    ref_abs=$(realpath "$ref")
    ref_name=$(basename "$ref")
    
    if [[ -f "$ref" ]]; then
      ref_size=$(du -h "$ref" | cut -f1)
      echo "  $ref_name ($ref_size)"
    else
      echo "  $ref_name"
    fi
    
    # Collect barcodes for this reference
    barcodes_for_ref=""
    for i in "${!REFERENCES[@]}"; do
      if [[ "${REFERENCES[$i]}" == "$ref" ]]; then
        barcodes_for_ref+="${BARCODES[$i]} "
      fi
    done
    
    # Sort and display barcodes
    sorted_barcodes=$(echo "$barcodes_for_ref" | tr ' ' '\n' | sort | tr '\n' ' ')
    echo "    → Barcodes: $sorted_barcodes"
    
    # Update references to absolute paths
    for i in "${!REFERENCES[@]}"; do
      if [[ "${REFERENCES[$i]}" == "$ref" ]]; then
        REFERENCES[$i]="$ref_abs"
      fi
    done
  done
  echo ""
fi

echo "Absolute paths:"
echo "  Input:  $input_dir"
echo "  Output: $output_dir"
if [[ -n "${ref_genome:-}" ]]; then
  echo "  Ref:    $ref_genome"
fi
echo ""

# Create output directory structure
echo "Creating output directory structure..."
mkdir -p "$output_dir/logs"
log_file="$output_dir/logs/pipeline_master_log_${timestamp}.txt"
echo "  [OK] Created: $output_dir/logs"
echo "  [OK] Master log: $log_file"
echo ""

# Simple logging - append output to log file
# Note: For full logging, run as: bash script.sh 2>&1 | tee logfile.txt

# Validate reference genome(s)
if [[ -n "$multi_mapping" ]]; then
  echo "Validating reference genomes (multi-species mode)..."
  # References already checked in parse_multi_mapping()
  echo "  [OK] All reference genomes validated"
  echo ""
else
  echo "Checking reference genome..."
  if [[ ! -f "$ref_genome" ]]; then
    echo "[ERROR] ERROR: Reference genome not found: $ref_genome"
    exit 1
  fi
  ref_size=$(du -h "$ref_genome" | cut -f1)
  echo "  [OK] Reference found: $ref_genome ($ref_size)"
  echo ""
fi

# Sample name for an aligned BAM. 001_ROCK_EDI_b01.minimap.bam becomes 001_ROCK_EDI_b01.
aligned_sample_name() {
  local base
  base=$(basename "$1")
  base="${base%.bam}"
  base="${base%.minimap}"
  printf '%s\n' "$base"
}

# True when the BAM header contains a reference sequence.
bam_has_sq() {
  samtools view -H "$1" 2>/dev/null | grep -q '^@SQ'
}

echo "Scanning for samples..."

# Function to extract barcode number from BAM file
get_barcode() {
  local bam="$1"
  # Use subshell to avoid SIGPIPE from head closing pipe early
  local bc_full=$(set +o pipefail; samtools view "$bam" 2>/dev/null | head -1 | grep -oE 'barcode[0-9]+' || true)
  if [[ -n "$bc_full" ]]; then
    # Extract just the number and format as b## (e.g., barcode03 -> b03)
    local bc_num=$(echo "$bc_full" | grep -oE '[0-9]+')
    # Use 10# to force decimal interpretation (avoids octal issues with leading zeros)
    printf "b%02d" "$((10#$bc_num))"
  else
    echo "b00"
  fi
}

# Determine which directories to scan
sample_dirs=""
if [[ "$aligned_bams" == true ]]; then
  echo "  Pipeline: aligned BAMs, modkit only"
  echo "  Subdirectories are ignored"
  if [[ -f "$input_dir" ]]; then
    case "$input_dir" in
      *.bam)
        echo "    $(basename "$input_dir")"
        sample_dirs="$input_dir"
        ;;
      *)
        echo ""
        echo "[ERROR] --aligned-bams expected a .bam file or a directory of .bam files"
        echo "  Given: $input_dir"
        exit 1
        ;;
    esac
  elif [[ -d "$input_dir" ]]; then
    shopt -s nullglob
    for aligned_bam in "$input_dir"/*.bam; do
      echo "    $(basename "$aligned_bam")"
      sample_dirs="$sample_dirs$aligned_bam"$'\n'
    done
    shopt -u nullglob
    sample_dirs=$(printf '%s\n' "$sample_dirs" | grep -v '^$' | sort)
  else
    echo ""
    echo "[ERROR] Input not found: $input_dir"
    exit 1
  fi

  if [[ -z "$sample_dirs" ]]; then
    echo ""
    echo "[ERROR] No aligned BAM files found in $input_dir"
    echo "  Put *.bam files directly in that directory. Subdirectories are ignored."
    echo "  An unaligned merged BAM cannot be used with --aligned-bams."
    exit 1
  fi
else
echo "  Include fail directory: $include_fail"

# Build list of base directories to scan
scan_dirs=""
if [[ -d "$input_dir/pass" ]]; then
  scan_dirs="$input_dir/pass"
  echo "  Found pass/ directory"
fi
if [[ "$include_fail" == true && -d "$input_dir/fail" ]]; then
  scan_dirs="$scan_dirs $input_dir/fail"
  echo "  Including fail/ directory"
fi

# If no pass/fail structure, try scanning input_dir directly
if [[ -z "$scan_dirs" ]]; then
  echo "  No pass/fail structure found, scanning input directory directly"
  scan_dirs="$input_dir"
fi

# Find all sample subdirectories containing BAM files
sample_dirs=""
for base_dir in $scan_dirs; do
  echo "  Scanning: $base_dir"
  for sample_dir in "$base_dir"/*/; do
    if [[ -d "$sample_dir" ]]; then
      dir_name=$(basename "$sample_dir")
      # Skip unclassified, logs, and other non-sample directories
      if [[ "$dir_name" =~ ^(unclassified|logs|tmp)$ ]]; then
        continue
      fi
      # Skip barcode## directories unless --include-empty-barcodes is set
      if [[ "$include_barcodes" == false && "$dir_name" =~ ^barcode[0-9]+$ ]]; then
        echo "    $dir_name: skipped (use --include-empty-barcodes to include)"
        continue
      fi
      # Check if directory contains BAM files
      bam_count=$(find "$sample_dir" -maxdepth 1 -name '*.bam' 2>/dev/null | wc -l | tr -d ' ')
      if [[ "$bam_count" -gt 0 ]]; then
        echo "    $dir_name: $bam_count BAM files"
        sample_dirs="$sample_dirs$sample_dir"$'\n'
      fi
    fi
  done
done

sample_dirs=$(echo "$sample_dirs" | grep -v '^$' | sort)

if [[ -z "$sample_dirs" ]]; then
  echo ""
  echo "[ERROR] ERROR: No sample directories found containing BAM files"
  echo ""
  echo "Searched in: $input_dir"
  echo "Expected structure: pass/<barcode_dirs>/*.bam"
  echo ""
  echo "Directory contents:"
  ls -la "$input_dir" | head -20
  if [[ -d "$input_dir/pass" ]]; then
    echo ""
    echo "Contents of pass/:"
    ls -la "$input_dir/pass" | head -10
  fi
  echo ""
  echo "Please check:"
  echo "  1. Input directory is correct"
  echo "  2. pass/ directory exists with barcode subdirectories"
  echo "  3. BAM files are directly in barcode directories"
  exit 1
fi
fi

sample_count=$(echo "$sample_dirs" | grep -c .)
echo "[OK] Found $sample_count sample(s) to process"
echo ""

if [[ "$dry_run" == true ]]; then
  echo "=========================================="
  echo "DRY RUN MODE - No processing will occur"
  echo "=========================================="
  echo ""
fi

# Process each sample
sample_num=0
processed_samples=""
for sample_path in $sample_dirs; do
  sample_num=$((sample_num + 1))

  if [[ "$aligned_bams" == true ]]; then
    input_bam="$sample_path"
    sample_name=$(aligned_sample_name "$input_bam")
    input_dir_name=$(basename "$input_bam")
    if [[ -n "$multi_mapping" ]]; then
      barcode_tag=$(set +o pipefail; samtools view "$input_bam" 2>/dev/null | head -1 | grep -oE 'barcode[0-9]+' || true)
      if [[ -z "$barcode_tag" ]]; then
        echo "=========================================="
        echo "Sample $sample_num of $sample_count (v${VERSION})"
        echo "=========================================="
        echo "[WARN] No barcode tag in $input_dir_name"
        echo "  Multi-species mode needs a barcode in the BAM. Skipping."
        echo ""
        continue
      fi
      barcode=$(get_barcode "$input_bam")
    else
      barcode=""
    fi
  else
    input_bam=""
    # Extract sample information
    input_dir_name=$(basename "$sample_path")

    # Get first BAM file to extract barcode (disable pipefail to avoid SIGPIPE)
    first_bam=$(set +o pipefail; find "$sample_path" -maxdepth 1 -name '*.bam' 2>/dev/null | head -1 || true)
    if [[ -n "$first_bam" ]]; then
      barcode=$(get_barcode "$first_bam")
    else
      barcode="b00"
    fi

    # Create sample name with barcode suffix (e.g., L07_HUMB_LAB_b03)
    sample_name="${input_dir_name}_${barcode}"
  fi

  echo "=========================================="
  echo "Sample $sample_num of $sample_count (v${VERSION})"
  echo "=========================================="
  if [[ "$aligned_bams" == true ]]; then
    echo "Mode:          aligned BAMs, modkit only"
    echo "Aligned BAM:   $input_bam"
    echo "Sample ID:     $sample_name"
  else
    echo "Input dir:     $input_dir_name"
    echo "Barcode:       $barcode"
    echo "Sample ID:     $sample_name"
    echo "Input path:    $sample_path"
  fi
  echo "Threads:       $THREADS"
  
  # Determine which reference to use for this sample
  if [[ -n "$multi_mapping" ]]; then
    # Multi-species mode: get reference based on barcode
    sample_ref=$(get_reference_for_sample "$barcode")
    if [[ $? -ne 0 ]]; then
      echo "[ERROR] Failed to find reference for $sample_name ($barcode)"
      echo "  Skipping this sample..."
      echo ""
      continue
    fi
    ref_name=$(basename "$sample_ref")
    echo "Reference:     $ref_name (barcode $barcode)"
  else
    # Single-species mode: use global reference
    sample_ref="$ref_genome"
    ref_name=$(basename "$ref_genome")
    echo "Reference:     $ref_name"
  fi
  echo ""

  if [[ "$aligned_bams" == true ]] && ! bam_has_sq "$input_bam"; then
    echo "[WARN] No @SQ line in the header of $input_dir_name"
    echo "  This BAM is not aligned. An unaligned merged BAM cannot be used here."
    echo "  Skipping this sample."
    echo ""
    continue
  fi

  # Create output directory (flat structure)
  out_sample_dir="$output_dir/$sample_name"
  mkdir -p "$out_sample_dir"
  mkdir -p "$output_dir/logs"

  sample_log="$output_dir/logs/${sample_name}.log"

  echo "Output dir:    $out_sample_dir"
  echo "Sample log:    $sample_log"
  echo ""

  if [[ "$dry_run" == true ]]; then
    if [[ "$aligned_bams" == true ]]; then
      echo "[DRY RUN] Would run modkit on $input_bam using $ref_name"
    else
      echo "[DRY RUN] Would process this sample with $THREADS threads using $ref_name"
    fi
    echo ""
    continue
  fi

  if [[ "$aligned_bams" == true ]]; then
    echo "Using the supplied aligned BAM. Merge and alignment are skipped."
    minimap_bam="$input_bam"
    if [[ ! -f "${minimap_bam}.bai" && ! -f "${minimap_bam}.csi" ]]; then
      echo "  No index found. Building one with samtools index."
      samtools index -@${THREADS} "$minimap_bam" 2>>"$sample_log"
    fi
    echo ""
  else
  # Step 1: Merge BAM files (these contain methylation tags from basecalling)
  echo "Step 1/4: Merging BAM files with methylation tags"
  merged_bam="$out_sample_dir/${sample_name}.merged.bam"

  if [[ "$force" != true && -s "$merged_bam" && $(file_size_bytes "$merged_bam") -gt 100000000 ]]; then
    merged_size=$(du -h "$merged_bam" | cut -f1)
    echo "  [OK] Already exists: $merged_bam ($merged_size)"
  else
    echo "  Searching for BAM files..."
    temp_bam_list="$out_sample_dir/bam_list.txt"
    find "$sample_path" -name '*.bam' > "$temp_bam_list"

    bam_count=$(wc -l < "$temp_bam_list")
    echo "    Found: $bam_count BAM files"

    if [[ $bam_count -eq 0 ]]; then
      echo "  [WARN]  WARNING: No BAM files found - skipping this sample"
      echo ""
      continue
    fi

    echo "  Merging BAM files (using $THREADS threads)..."
    merged_partial="${merged_bam}.partial.${timestamp}.$$"
    {
      read firstbam
      if ! samtools view -@${THREADS} -H "$firstbam" > /dev/null 2>>"$sample_log"; then
        echo "  [ERROR] BAM header problem in $firstbam" | tee -a "$sample_log"
        continue
      fi
      samtools view -@${THREADS} -h "$firstbam"
      while read bam; do
        if samtools view -@${THREADS} "$bam" > /dev/null 2>>"$sample_log"; then
          samtools view -@${THREADS} "$bam"
        else
          echo "  [WARN]  Corrupt BAM: $bam" | tee -a "$sample_log"
        fi
      done
    } < "$temp_bam_list" | samtools view -@${THREADS} -ubS - | samtools sort -@${THREADS} -o "$merged_partial" -

    mv -f "$merged_partial" "$merged_bam"
    samtools index -@${THREADS} "$merged_bam"
    merged_size=$(du -h "$merged_bam" | cut -f1)
    echo "  [OK] BAMs merged: $merged_bam ($merged_size)"
  fi
  echo ""

  # Step 2: Alignment with minimap2 (preserving methylation tags)
  echo "Step 2/4: Aligning reads with minimap2"
  minimap_bam="$out_sample_dir/${sample_name}.minimap.bam"

  if [[ "$force" != true && -s "$minimap_bam" && $(file_size_bytes "$minimap_bam") -gt 100000000 ]]; then
    bam_size=$(du -h "$minimap_bam" | cut -f1)
    echo "  [OK] Already exists: $minimap_bam ($bam_size)"
  else
    echo "  Running minimap2 alignment (preserving methylation tags)..."
    echo "    Reference: $ref_name"
    echo "    Threads: $THREADS"
    echo "    Mode: map-ont with -y (copy tags)"
    start_time=$(date +%s)

    minimap_partial="${minimap_bam}.partial.${timestamp}.$$"
    samtools fastq -@${THREADS} -T MM,ML "$merged_bam" | \
      minimap2 -ax map-ont -t ${THREADS} -y --secondary=no "$sample_ref" - 2>>"$sample_log" | \
      samtools view -@${THREADS} -S -b - 2>>"$sample_log" | \
      samtools sort -@${THREADS} -o "$minimap_partial" -T "$out_sample_dir/reads.tmp.${timestamp}.$$" - 2>>"$sample_log"

    mv -f "$minimap_partial" "$minimap_bam"
    echo "  Creating BAM index..."
    samtools index -@${THREADS} "$minimap_bam" 2>>"$sample_log"

    end_time=$(date +%s)
    elapsed=$((end_time - start_time))
    bam_size=$(du -h "$minimap_bam" | cut -f1)

    echo "  [OK] Alignment complete ($elapsed seconds)"
    echo "    Output: $minimap_bam ($bam_size)"
  fi
  echo ""
  fi

  # Step 3: Methylation calling
  if [[ "$aligned_bams" == true ]]; then
    echo "Calling methylation with modkit"
    echo "  Source BAM: $minimap_bam"
  else
    echo "Step 3/4: Calling methylation with modkit"
  fi
  methyl_bed="$out_sample_dir/${sample_name}.CpG.bed"

  if [[ "$force" != true && -s "$methyl_bed" && $(file_size_bytes "$methyl_bed") -gt 100000 ]]; then
    bed_size=$(du -h "$methyl_bed" | cut -f1)
    echo "  [OK] Already exists: $methyl_bed ($bed_size)"
  else
    echo "  Running modkit pileup..."
    echo "    Reference: $ref_name"
    echo "    Mode: CpG methylation"
    echo "    Strands: $strand_mode"
    echo "    Threads: $THREADS"
    start_time=$(date +%s)

    # --modified-bases 5mC = 5-methylcytosine at CpG.
    # --combine-mods pools modification codes. It does not pool strands.
    # --combine-strands is added only when the user asks for one row per CpG.
    bed_partial="${methyl_bed}.partial.${timestamp}.$$"
    modkit_cmd=(
      modkit pileup "$minimap_bam" "$bed_partial"
      --cpg --ref "$sample_ref" -t "${THREADS}"
      --combine-mods --modified-bases 5mC
    )
    if [[ "$combine_strands" == true ]]; then
      modkit_cmd+=(--combine-strands)
    fi
    if "${modkit_cmd[@]}" 2>&1 | tee -a "$sample_log"; then
      mv -f "$bed_partial" "$methyl_bed"
      end_time=$(date +%s)
      elapsed=$((end_time - start_time))
      bed_size=$(du -h "$methyl_bed" | cut -f1)
      echo "  [OK] Methylation calling complete ($elapsed seconds)"
      echo "    Output: $methyl_bed ($bed_size)"
    else
      echo "  [WARN]  modkit pileup failed - check sample log for details"
      echo "  The previous CpG bed was left in place."
      echo "  Continuing to next sample..."
      continue
    fi
  fi
  echo ""

  if [[ "$aligned_bams" == true ]]; then
    echo "Qualimap is not run in aligned-BAM mode."
    echo ""
  else
  # Step 4: Quality control
  echo "Step 4/4: Running Qualimap QC"
  qualimap_dir="$out_sample_dir/qualimap"

  if [[ "$force" != true && -f "$qualimap_dir/qualimapReport.html" ]]; then
    echo "  [OK] Already exists: $qualimap_dir/qualimapReport.html"
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

    qualimap_rerun="$out_sample_dir/qualimap.rerun.${timestamp}.$$"
    if qualimap bamqc -bam "$minimap_bam" -nw 5000 -nt ${qualimap_threads} \
      -c -outdir "$qualimap_rerun" >>"$sample_log" 2>&1; then
      if [[ -d "$qualimap_dir" ]]; then
        qualimap_previous="$out_sample_dir/qualimap.previous-${timestamp}"
        if [[ -e "$qualimap_previous" ]]; then
          qualimap_previous="${qualimap_previous}-$$"
        fi
        mv "$qualimap_dir" "$qualimap_previous"
      fi
      mv "$qualimap_rerun" "$qualimap_dir"
      end_time=$(date +%s)
      elapsed=$((end_time - start_time))
      echo "  [OK] QC complete ($elapsed seconds)"
      echo "    Report: $qualimap_dir/qualimapReport.html"
    else
      echo "  [WARN]  Qualimap failed - check sample log for details"
      echo "  The previous Qualimap report was left in place."
      echo "  Continuing to next sample..."
      continue
    fi
  fi
  echo ""
  fi

  echo "[OK] Sample $sample_name complete!"
  echo ""
  
  # Track processed sample (for multi-species summary)
  processed_samples="$processed_samples$sample_name:$barcode "

done

# Multi-species processing summary
if [[ -n "$multi_mapping" && -n "$processed_samples" ]]; then
  echo "=========================================="
  echo "Multi-Species Processing Summary"
  echo "=========================================="
  
  # Get unique references processed
  declare -a processed_refs=()
  declare -a processed_ref_samples=()
  declare -a processed_ref_counts=()
  
  for sample_info in $processed_samples; do
    sample_name="${sample_info%:*}"
    barcode="${sample_info#*:}"
    
    # Find reference for this barcode
    sample_ref=""
    for i in "${!BARCODES[@]}"; do
      if [[ "${BARCODES[$i]}" == "$barcode" ]]; then
        sample_ref="${REFERENCES[$i]}"
        break
      fi
    done
    
    if [[ -n "$sample_ref" ]]; then
      ref_name=$(basename "$sample_ref")
      
      # Find or add reference in tracking arrays
      ref_idx=-1
      if [[ ${#processed_refs[@]} -gt 0 ]]; then
        for i in "${!processed_refs[@]}"; do
          if [[ "${processed_refs[$i]}" == "$ref_name" ]]; then
            ref_idx=$i
            break
          fi
        done
      fi
      
      if [[ $ref_idx -eq -1 ]]; then
        # New reference
        processed_refs+=("$ref_name")
        processed_ref_samples+=("$sample_name ")
        processed_ref_counts+=(1)
      else
        # Existing reference
        processed_ref_samples[$ref_idx]+="$sample_name "
        processed_ref_counts[$ref_idx]=$((${processed_ref_counts[$ref_idx]} + 1))
      fi
    fi
  done
  
  echo "Samples processed by reference genome:"
  echo ""
  for i in "${!processed_refs[@]}"; do
    ref_name="${processed_refs[$i]}"
    count="${processed_ref_counts[$i]}"
    samples="${processed_ref_samples[$i]}"
    
    echo "  $ref_name: $count sample(s)"
    # Display sample names (limit to 10 per line for readability)
    sample_array=($samples)
    for j in "${!sample_array[@]}"; do
      if [[ $((j % 10)) -eq 0 && $j -gt 0 ]]; then
        echo ""
        echo -n "    "
      fi
      if [[ $j -eq 0 ]]; then
        echo -n "    "
      fi
      echo -n "${sample_array[$j]} "
    done
    echo ""
    echo ""
  done
fi

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

# Step 5: Generate summary statistics and plots
echo "=========================================="
echo "Step 5/5: Generating Summary Report"
echo "=========================================="

SUMMARY_SCRIPT="$SCRIPT_DIR/generate_summary.R"

if [[ -f "$SUMMARY_SCRIPT" ]]; then
  if command -v Rscript &> /dev/null; then
    if [[ "$combine_strands" == true ]]; then
      echo "Strand-bias plots read column 6. With --combine-strands every row is the plus-strand coordinate, so that plot does not compare strands."
    fi
    echo "Running summary statistics and plot generation..."
    if [[ "$expanded_plots" == true ]]; then
      Rscript "$SUMMARY_SCRIPT" "$output_dir" --expanded-plots 2>&1 | tee -a "$log_file"
    else
      Rscript "$SUMMARY_SCRIPT" "$output_dir" 2>&1 | tee -a "$log_file"
    fi
  else
    echo "[WARN]  Rscript not found - skipping summary report"
    echo "   To generate summary, run: Rscript generate_summary.R $output_dir"
  fi
else
  echo "[WARN]  generate_summary.R not found in $SCRIPT_DIR"
  echo "   Download from: https://github.com/markusdrag/NanoporeToBED-Pipeline"
fi

echo ""
echo "Citation: Drag et al. (2025) bioRxiv 2025.04.11.648151"
echo "--- All done!"
