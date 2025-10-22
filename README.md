<div align="center">
  <img src="nanoporetobed.png" alt="NanoporeToBED Pipeline Logo" width="400">
</div>

#

A comprehensive pipeline for processing Oxford Nanopore Technologies sequencing data with base modifications (5mC) to generate methylation BED files and quality metrics.

## Overview

This pipeline processes modified BAM (modBAM) files from Oxford Nanopore Technologies sequencing runs with methylation calling. It takes the output from Guppy basecaller with methylation awareness and produces:
- Merged and aligned BAM files
- CpG methylation BED files  
- Quality control reports

**The output BED files are directly ready for use with the MethylSense package**, which provides powerful tools for differential methylation discovery and machine learning modelling of differentially methylated regions (DMRs). Once your BED files are generated, head over to the MethylSense repository to:
- Identify differentially methylated regions between conditions
- Build predictive models for diagnostic and prognostic testing
- Perform biomarker discovery using methylation patterns
- Create clinical classifiers based on cfDNA methylation signatures

MethylSense seamlessly integrates with this pipeline's output for comprehensive methylation analysis from raw Nanopore data to clinical insights.

## Prerequisites

### Input Data Structure
The pipeline expects data from a Guppy basecaller run with methylation calling:
```bash
guppy_basecaller --disable_pings --compress_fastq \
  -c dna_r10.4.1_e8.2_400bps_modbases_5mc_cg_hac.cfg \
  --num_callers 4 \
  -i pod5_skip \
  -s fastq_gpu_hac_mod \
  -x 'auto' --bam_out --recursive --min_qscore 7 \
  --barcode_kits 'SQK-NBD114-24'
```

### Expected Input Directory Structure
```
input_dir/
├── SRR000001/
│   └── fastq_gpu_hac_only_5mc_mod/
│       └── pass/
│           ├── sample1_metadata_info/
│           │   ├── *.bam           # Modified BAM files with MM/ML tags
│           │   └── *.fastq.gz      # FASTQ files (optional)
│           └── sample2_metadata_info/
│               ├── *.bam
│               └── *.fastq.gz
├── SRR000002/
│   └── fastq_gpu_hac_only_5mc_mod/
│       └── pass/
│           └── sample3_metadata_info/
│               └── *.bam
└── ...
```

**Directory structure notes:**
- Top level: SRR codes (sequencing run identifiers)
- Second level: Date or run identifier (automatically extracted from directory structure)
- Sample directories: Any naming convention with metadata separated by underscores

## Installation

### Quick Automated Setup (Recommended for HPC)

```bash
# Clone the repository
git clone https://github.com/markusdrag/NanoporeToBED-Pipeline.git
cd NanoporeToBED-Pipeline

# Run the setup script
bash setup.sh

# Or specify a custom installation directory
bash setup.sh /path/to/install/location
```

The setup script will:
- Create all necessary directories
- Install the pipeline script
- Create the conda environment automatically
- Set up example files and documentation

### Manual Installation

#### 1. Download the Pipeline

```bash
# Clone the repository
git clone https://github.com/markusdrag/NanoporeToBED-Pipeline.git
cd NanoporeToBED-Pipeline

# Make the script executable
chmod +x NanoporeToBED.sh
```

Alternatively, download just the script:
```bash
wget https://raw.githubusercontent.com/markusdrag/NanoporeToBED-Pipeline/main/NanoporeToBED.sh
chmod +x NanoporeToBED.sh
```

#### 2. Create Conda Environment

The setup script creates this automatically, but for manual setup:

```yaml
name: nanopore_methylation
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - python>=3.8
  - samtools>=1.17
  - minimap2>=2.26
  - bioconda::ont-modkit>=0.3.0
  - qualimap>=2.3
  - pigz
  - parallel
```

Install the environment:
```bash
# If you cloned the repository
conda env create -f environment.yml
conda activate nanopore_methylation

# Or if downloading script only, create environment directly
conda create -n nanopore_methylation -c conda-forge -c bioconda \
  samtools>=1.17 minimap2>=2.26 ont-modkit>=0.3.0 qualimap>=2.3 \
  pigz parallel python>=3.8
conda activate nanopore_methylation
```

#### 3. Alternative: Using Micromamba
```bash
micromamba create -n nanopore_methylation -c conda-forge -c bioconda \
  samtools minimap2 ont-modkit qualimap pigz parallel python>=3.8
micromamba activate nanopore_methylation
```

### Quick Start on HPC

For the fastest setup on your HPC, use the automated setup:
```bash
# Clone and setup
git clone https://github.com/markusdrag/NanoporeToBED-Pipeline.git
cd NanoporeToBED-Pipeline
bash setup.sh

# Then activate and run
conda activate nanopore_methylation
sbatch NanoporeToBED.sh \
  -i /data/nanopore/fastq_gpu_hac_mod \
  -o /data/nanopore/methylation_results \
  -ref /data/references/genome.fna \
  -t 32
```

For manual setup:
```bash
# 1. Get the pipeline
git clone https://github.com/markusdrag/NanoporeToBED-Pipeline.git
cd NanoporeToBED-Pipeline

# 2. Load your HPC's module system (if available)
module load conda  # or module load miniconda3

# 3. Create and activate environment
conda env create -f environment.yml
conda activate nanopore_methylation

# 4. Test the script
./NanoporeToBED.sh -h

# 5. Run on your data
sbatch NanoporeToBED.sh \
  -i /data/nanopore/fastq_gpu_hac_mod \
  -o /data/nanopore/methylation_results \
  -ref /data/references/genome.fna \
  -t 32
```

## Usage

### Basic Command
```bash
sbatch NanoporeToBED.sh \
  -i /data/nanopore/fastq_gpu_hac_mod \
  -o /data/nanopore/methylation_results \
  -ref /data/references/genome.fna \
  -t 40
```

Where the input directory should contain your SRR folders from the Guppy output:
```
/data/nanopore/fastq_gpu_hac_mod/
├── SRR000001/
│   └── fastq_gpu_hac_only_5mc_mod/pass/
├── SRR000002/
│   └── fastq_gpu_hac_only_5mc_mod/pass/
└── ...
```

### Parameters

| Flag | Long Form | Description | Required | Default |
|------|-----------|-------------|----------|---------|
| `-i` | `--input` | Input directory containing SRR folders with barcoded samples | Yes | - |
| `-o` | `--output` | Output directory for processed data | Yes | - |
| `-ref` | `--reference` | Path to reference genome FASTA file (.fna, .fa, or .fasta) | Yes | - |
| `-t` | `--threads` | Number of threads to use for processing | No | 40 |
| | `--dry-run` | Run in test mode without processing | No | false |
| `-h` | `--help` | Show help message | No | - |

### Reference Genome Format
The reference genome should be a single FASTA file:
```
/path/to/reference/
└── genome.fna    # Or genome.fa, genome.fasta
```

## Output Structure

```
output_dir/
├── SRR000001/
│   ├── 20240315/                              # Date folder
│   │   ├── sample1_metadata_info/
│   │   │   ├── sample1_metadata_info.merged.bam         # Merged BAM
│   │   │   ├── sample1_metadata_info.merged.bam.bai
│   │   │   ├── sample1_metadata_info.minimap.bam        # Aligned BAM
│   │   │   ├── sample1_metadata_info.minimap.bam.bai
│   │   │   ├── sample1_metadata_info.CpG.bed           # Methylation calls
│   │   │   ├── qualimap/                               # QC reports
│   │   │   │   ├── qualimapReport.html
│   │   │   │   ├── genome_results.txt
│   │   │   │   └── raw_data_qualimapReport/
│   │   │   └── bam_list.txt                            # Processing manifest
│   │   └── sample2_metadata_info/
│   │       └── ...
├── SRR000002/
│   └── ...
└── logs/
    ├── pipeline_master_log_YYYYMMDD_HHMMSS.txt        # Master pipeline log
    ├── SRR000001/
    │   └── 20240315/
    │       ├── sample1_metadata_info.log               # Sample-specific log
    │       └── sample2_metadata_info.log
    └── ...
```

### Output File Descriptions

- **`.merged.bam`**: Concatenated BAM files from all sequencing chunks with methylation tags preserved
- **`.minimap.bam`**: Re-aligned BAM files to reference genome
- **`.CpG.bed`**: BED file with CpG methylation frequencies
  - Format: chromosome, start, end, modification_frequency, coverage, strand
  - Compatible with MethylSense pipeline (see separate repository) and other methylation analysis tools
- **`qualimap/`**: HTML quality reports with coverage statistics and alignment metrics

## Pipeline Steps

1. **BAM Merging**: Combines multiple BAM files whilst preserving MM/ML methylation tags
2. **Alignment**: Re-aligns reads to reference genome using minimap2 with tag preservation
3. **Methylation Calling**: Extracts CpG methylation frequencies using modkit
4. **Quality Control**: Generates comprehensive QC reports using Qualimap

## Running Tips

### SLURM Configuration
Adjust SLURM parameters in the script header based on your data:
```bash
#SBATCH -c 40          # CPU cores (adjust based on availability)
#SBATCH --mem 192g     # Memory (scale with data size)
#SBATCH --time=72:00:00 # Time limit (depends on dataset size)
#SBATCH --account YourAccount # Your HPC account
```

### Performance Optimisation
- **Thread usage**: The `-t` parameter sets thread count (default: 40). Will automatically reduce if exceeds SLURM allocation
- **Memory**: ~4-8 GB per thread is recommended
- **Storage**: Ensure 3-5x input data size for temporary files
- **Large datasets**: Consider processing in batches or increasing time allocation

### Troubleshooting

1. **Insufficient memory**: Reduce thread count or increase memory allocation
2. **Corrupted BAM files**: Script automatically skips problematic BAMs with warnings
3. **Missing reference**: Verify reference genome path and file exists
4. **Timeout issues**: Extend time limit or process fewer samples
5. **No samples found**: Check input directory structure matches expected pattern (see error message for searched patterns)

### Monitoring Progress
```bash
# Check job status
squeue -u $USER

# Monitor SLURM output log
tail -f NanoporeToBED.out

# Check master log for overall progress
tail -f output_dir/logs/pipeline_master_log_*.txt

# Check individual sample logs
tail -f output_dir/logs/SRR*/*/sample_name.log
```

## Quality Control Metrics

The pipeline generates several QC checkpoints:
- File size validation (>100 MB threshold for merged files)
- BAM header integrity checks
- Alignment statistics via Qualimap
- Methylation coverage in BED files

### Key Metrics to Review
- **Coverage depth**: Check in Qualimap reports
- **Mapping rate**: Verify alignment efficiency
- **Methylation sites**: Number of CpG sites covered
- **Read length distribution**: Assess data quality

## Advanced Configuration

### Custom Reference Genomes
The pipeline works with any reference genome in FASTA format:
```bash
# Index your reference (optional, minimap2 will do this automatically)
minimap2 -d reference.mmi reference.fna

# Use in pipeline
sbatch NanoporeToBED.sh -ref /path/to/reference.fna ...
```

### Batch Processing
For multiple libraries, create a wrapper script:
```bash
#!/bin/bash
# Process multiple Guppy output directories
for lib in SRR00000{1..5}; do
  sbatch NanoporeToBED.sh \
    -i /data/nanopore/fastq_gpu_hac_mod \
    -o /data/nanopore/methylation_results/$lib \
    -ref /data/references/genome.fna \
    -t 40
done
```

## Citation

If you use this pipeline, please cite:

**Our methodology paper:**
- Drag, M.H., Hvilsom, C., Poulsen, L.L., Jensen, H.E., Tahas, S.A., Leineweber, C., Cray, C., Bertelsen, M.F., Bojesen, A.M. (2025). New high accuracy diagnostics for avian Aspergillus fumigatus infection using Nanopore methylation sequencing of host cell-free DNA and machine learning prediction. *bioRxiv* 2025.04.11.648151. https://doi.org/10.1101/2025.04.11.648151

**Software tools:**
- Oxford Nanopore Technologies modkit
- Minimap2: Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences. *Bioinformatics*, 34(18), 3094-3100.
- Samtools: Danecek, P., et al. (2021). Twelve years of SAMtools and BCFtools. *GigaScience*, 10(2), giab008.
- Qualimap: Okonechnikov, K., et al. (2016). Qualimap 2: advanced multi-sample quality control for high-throughput sequencing data. *Bioinformatics*, 32(2), 292-294.

## Licence

MIT Licence (see LICENCE file)

## Contact

- **Lead Developer**: Markus Hodal Drag
- **Email**: markus.drag@sund.ku.dk
- **Institution**: University of Copenhagen
- **GitHub**: https://github.com/markusdrag
- **ORCID**: https://orcid.org/0000-0002-7412-6402

For questions, bug reports, or feature requests, please:
- Open an issue on GitHub: https://github.com/markusdrag/NanoporeToBED-Pipeline/issues
- Or contact via email for collaboration enquiries
