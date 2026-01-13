#!/usr/bin/env Rscript

# NanoporeToBED Pipeline - Summary Statistics and Plots
# Version: 1.6.0
# Author: Markus Hodal Drag
#
# Changes in 1.6.0:
#   - Added comprehensive qualimap BAM statistics (27 fields)
#   - Improved aggregation for re-sequenced samples
#   - All plots use aggregated data when duplicate sample IDs exist
#   - New QC plots: mapping quality, error rate, GC content

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
})

# =============================================================================
# QUALIMAP PARSER - Extract all BAM statistics from genome_results.txt
# =============================================================================

parse_qualimap <- function(qualimap_file) {
  # Initialize all fields as NA
  result <- list(
    number_of_reads = NA_real_,
    number_of_mapped_reads = NA_real_,
    mapped_reads_percentage = NA_real_,
    number_of_supplementary_alignments = NA_real_,
    number_of_mapped_bases = NA_real_,
    number_of_sequenced_bases = NA_real_,
    number_of_duplicated_reads = NA_real_,
    duplication_rate = NA_real_,
    mean_insert_size = NA_real_,
    std_insert_size = NA_real_,
    median_insert_size = NA_real_,
    mean_mapping_quality = NA_real_,
    number_of_As = NA_real_,
    number_of_Cs = NA_real_,
    number_of_Ts = NA_real_,
    number_of_Gs = NA_real_,
    number_of_Ns = NA_real_,
    GC_percentage = NA_real_,
    general_error_rate = NA_real_,
    number_of_mismatches = NA_real_,
    number_of_insertions = NA_real_,
    mapped_reads_with_insertion_percentage = NA_real_,
    number_of_deletions = NA_real_,
    mapped_reads_with_deletion_percentage = NA_real_,
    homopolymer_indels = NA_real_,
    mean_coverage_bam = NA_real_,
    std_coverage_bam = NA_real_
  )

  if (!file.exists(qualimap_file)) {
    return(result)
  }

  lines <- tryCatch(readLines(qualimap_file, warn = FALSE), error = function(e) NULL)
  if (is.null(lines)) {
    return(result)
  }

  # Helper function to extract numeric value from a line
  extract_numeric <- function(pattern, lines, remove_commas = TRUE) {
    line <- grep(pattern, lines, value = TRUE, ignore.case = TRUE)
    if (length(line) == 0) {
      return(NA_real_)
    }
    # Extract the value after the = sign
    value_str <- sub("^.*=\\s*", "", line[1])
    value_str <- sub("\\s*\\(.*$", "", value_str) # Remove parenthetical info
    value_str <- trimws(value_str)
    if (remove_commas) value_str <- gsub(",", "", value_str)
    value_str <- gsub("[^0-9.\\-]", "", value_str)
    as.numeric(value_str)
  }

  # Helper for percentage extraction
  extract_percentage <- function(pattern, lines) {
    line <- grep(pattern, lines, value = TRUE, ignore.case = TRUE)
    if (length(line) == 0) {
      return(NA_real_)
    }
    # Look for pattern like (XX.XX%) or XX.XX%
    pct_match <- regmatches(line[1], regexpr("[0-9]+\\.?[0-9]*%", line[1]))
    if (length(pct_match) > 0) {
      return(as.numeric(gsub("%", "", pct_match)))
    }
    NA_real_
  }

  # Parse all fields
  result$number_of_reads <- extract_numeric("^\\s*number of reads\\s*=", lines)
  result$number_of_mapped_reads <- extract_numeric("^\\s*number of mapped reads\\s*=", lines)
  result$mapped_reads_percentage <- extract_percentage("number of mapped reads", lines)
  result$number_of_supplementary_alignments <- extract_numeric("number of supplementary alignments", lines)
  result$number_of_mapped_bases <- extract_numeric("number of mapped bases", lines)
  result$number_of_sequenced_bases <- extract_numeric("number of sequenced bases", lines)
  result$number_of_duplicated_reads <- extract_numeric("number of duplicated reads", lines)
  result$duplication_rate <- extract_percentage("duplication rate", lines)

  # Insert size stats
  result$mean_insert_size <- extract_numeric("mean insert size", lines)
  result$std_insert_size <- extract_numeric("std insert size", lines)
  result$median_insert_size <- extract_numeric("median insert size", lines)

  # Mapping quality
  result$mean_mapping_quality <- extract_numeric("mean mapping quality", lines)

  # Nucleotide content
  result$number_of_As <- extract_numeric("number of A's", lines)
  result$number_of_Cs <- extract_numeric("number of C's", lines)
  result$number_of_Ts <- extract_numeric("number of T's", lines)
  result$number_of_Gs <- extract_numeric("number of G's", lines)
  result$number_of_Ns <- extract_numeric("number of N's", lines)
  result$GC_percentage <- extract_numeric("GC percentage", lines)

  # Error rates
  result$general_error_rate <- extract_numeric("general error rate", lines)
  result$number_of_mismatches <- extract_numeric("number of mismatches", lines)

  # Insertions and deletions
  result$number_of_insertions <- extract_numeric("number of insertions", lines)
  result$mapped_reads_with_insertion_percentage <- extract_percentage("mapped reads with insertion percentage", lines)
  result$number_of_deletions <- extract_numeric("number of deletions", lines)
  result$mapped_reads_with_deletion_percentage <- extract_percentage("mapped reads with deletion percentage", lines)
  result$homopolymer_indels <- extract_numeric("homopolymer indels", lines)

  # Coverage (from BAM, not BED)
  result$mean_coverage_bam <- extract_numeric("mean coverageData", lines)
  if (is.na(result$mean_coverage_bam)) {
    result$mean_coverage_bam <- extract_numeric("mean coverage", lines)
  }
  result$std_coverage_bam <- extract_numeric("std coverageData", lines)
  if (is.na(result$std_coverage_bam)) {
    result$std_coverage_bam <- extract_numeric("std coverage", lines)
  }

  return(result)
}

# =============================================================================
# AGGREGATION FUNCTION - Combine re-sequenced samples
# =============================================================================

aggregate_samples <- function(df) {
  # Check for duplicate sample IDs
  if (!any(duplicated(df$sample))) {
    return(NULL) # No duplicates, no aggregation needed
  }

  cat("\n--- Detected re-sequenced samples ---\n")
  dup_samples <- unique(df$sample[duplicated(df$sample)])
  cat("Samples with multiple runs:", paste(dup_samples, collapse = ", "), "\n")

  # Define column types for aggregation
  sum_cols <- c(
    "cpg_sites", "hyper_meth_count", "hypo_meth_count",
    "number_of_reads", "number_of_mapped_reads",
    "number_of_supplementary_alignments", "number_of_mapped_bases",
    "number_of_sequenced_bases", "number_of_duplicated_reads",
    "number_of_As", "number_of_Cs", "number_of_Ts", "number_of_Gs", "number_of_Ns",
    "number_of_mismatches", "number_of_insertions", "number_of_deletions",
    "homopolymer_indels"
  )

  # Weighted average columns (weight by number_of_reads)
  wavg_cols <- c(
    "mean_methylation", "median_methylation", "mean_coverage", "low_cov_percent",
    "mapped_reads_percentage", "duplication_rate",
    "mean_insert_size", "median_insert_size", "mean_mapping_quality",
    "GC_percentage", "general_error_rate",
    "mapped_reads_with_insertion_percentage", "mapped_reads_with_deletion_percentage",
    "mean_coverage_bam"
  )

  # Perform aggregation
  agg_df <- df %>%
    group_by(sample) %>%
    summarise(
      n_runs = n(),
      # Sum columns
      across(all_of(intersect(sum_cols, names(df))), ~ sum(.x, na.rm = TRUE)),
      # Weighted average columns (weight by reads)
      across(
        all_of(intersect(wavg_cols, names(df))),
        ~ if (all(is.na(.x))) {
          NA_real_
        } else {
          weights <- number_of_reads
          weights[is.na(weights)] <- 1
          weighted.mean(.x, w = weights, na.rm = TRUE)
        }
      ),
      # Pooled std
      std_insert_size = if (all(is.na(std_insert_size))) NA_real_ else sqrt(mean(std_insert_size^2, na.rm = TRUE)),
      std_coverage_bam = if (all(is.na(std_coverage_bam))) NA_real_ else sqrt(mean(std_coverage_bam^2, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    # Recalculate percentages based on summed counts
    mutate(
      mapped_reads_percentage = ifelse(
        number_of_reads > 0,
        (number_of_mapped_reads / number_of_reads) * 100,
        NA_real_
      ),
      duplication_rate = ifelse(
        number_of_mapped_reads > 0,
        (number_of_duplicated_reads / number_of_mapped_reads) * 100,
        NA_real_
      )
    )

  return(agg_df)
}

# =============================================================================
# MAIN SCRIPT
# =============================================================================

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check for --summary_file flag (for regenerating from existing CSV)
summary_file_idx <- which(args == "--summary_file")
if (length(summary_file_idx) > 0 && (summary_file_idx + 1) <= length(args)) {
  input_summary_file <- args[summary_file_idx + 1]
  args <- args[-c(summary_file_idx, summary_file_idx + 1)]
} else {
  input_summary_file <- NULL
}

if (length(args) < 1) {
  cat("Usage: Rscript generate_summary.R <output_dir> [OPTIONS]\n")
  cat("\nOptions:\n")
  cat("  --expanded-plots     Generate additional distribution, QC, and comparative plots\n")
  cat("  --summary_file FILE  Use existing CSV instead of parsing BED files\n")
  cat("                       Note: Duplicate sample IDs are aggregated automatically\n")
  quit(status = 1)
}

output_dir <- args[1]
expanded_plots <- "--expanded-plots" %in% args

cat("==========================================\n")
cat("NanoporeToBED Summary Report v1.6.0\n")
cat("==========================================\n")
cat("Output directory:", output_dir, "\n")
cat("\n")
cat("Mode:\n")
if (expanded_plots) {
  cat("  [x] Basic plots (basic/)\n")
  cat("  [x] Distribution plots (distribution/)\n")
  cat("  [x] QC plots (qc/)\n")
  cat("  [x] Biological plots (biological/)\n")
  cat("  [x] Comparative plots (comparative/)\n")
} else {
  cat("  [x] Basic plots (basic/)\n")
  cat("  [ ] Distribution plots (use --expanded-plots)\n")
  cat("  [ ] QC plots (use --expanded-plots)\n")
  cat("  [ ] Biological plots (use --expanded-plots)\n")
  cat("  [ ] Comparative plots (use --expanded-plots)\n")
}
cat("\n")

# Initialize results list
results_list <- list()
all_meth_data <- list()
all_chr_data <- list()
all_strand_data <- list()

# Mode: Load from existing CSV or parse BED files
if (!is.null(input_summary_file)) {
  cat("Loading from existing summary file:", input_summary_file, "\n")
  results <- read.csv(input_summary_file, stringsAsFactors = FALSE)
  cat("Loaded", nrow(results), "samples\n")
} else {
  # Find all sample directories (those containing .CpG.bed files)
  bed_files <- list.files(output_dir, pattern = "\\.CpG\\.bed$", recursive = TRUE, full.names = TRUE)

  if (length(bed_files) == 0) {
    cat("ERROR: No .CpG.bed files found in", output_dir, "\n")
    quit(status = 1)
  }

  cat("Found", length(bed_files), "sample(s)\n\n")

  # Process each sample
  for (bed_file in bed_files) {
    sample_name <- gsub("\\.CpG\\.bed$", "", basename(bed_file))
    sample_dir <- dirname(bed_file)

    cat("Processing:", sample_name, "\n")

    # Initialize row
    row <- list(sample = sample_name)

    # Read BED file
    bed_lines <- tryCatch(readLines(bed_file, warn = FALSE), error = function(e) NULL)

    if (!is.null(bed_lines) && length(bed_lines) > 0) {
      bed_lines <- bed_lines[!grepl("^#", bed_lines) & nchar(bed_lines) > 0]

      if (length(bed_lines) > 0) {
        bed_data <- tryCatch(
          read.table(text = bed_lines, header = FALSE, stringsAsFactors = FALSE, fill = TRUE),
          error = function(e) NULL
        )

        if (!is.null(bed_data) && nrow(bed_data) > 0) {
          cat("  BED columns:", ncol(bed_data), "\n")

          # Extract methylation data
          if (ncol(bed_data) >= 12) {
            chr <- bed_data[, 1]
            strand <- bed_data[, 6]
            coverage <- as.numeric(bed_data[, 10])
            percent_mod <- as.numeric(bed_data[, 11])
            if (max(percent_mod, na.rm = TRUE) > 1) {
              methylation <- percent_mod / 100
            } else {
              methylation <- percent_mod
            }
          } else if (ncol(bed_data) >= 6) {
            chr <- bed_data[, 1]
            strand <- bed_data[, 6]
            score <- as.numeric(bed_data[, 5])
            methylation <- if (max(score, na.rm = TRUE) > 1) score / 100 else score
            coverage <- rep(NA, nrow(bed_data))
          } else {
            methylation <- NULL
          }

          if (!is.null(methylation)) {
            valid_meth <- !is.na(methylation) & methylation >= 0 & methylation <= 1
            if (sum(valid_meth) > 0) {
              methylation_clean <- methylation[valid_meth]
              coverage_clean <- coverage[valid_meth]
              chr_clean <- chr[valid_meth]
              strand_clean <- strand[valid_meth]

              row$cpg_sites <- length(methylation_clean)
              row$mean_methylation <- mean(methylation_clean, na.rm = TRUE)
              row$median_methylation <- median(methylation_clean, na.rm = TRUE)
              row$mean_coverage <- mean(coverage_clean, na.rm = TRUE)
              row$low_cov_percent <- sum(coverage_clean < 10, na.rm = TRUE) / row$cpg_sites * 100
              row$hyper_meth_count <- sum(methylation_clean > 0.8, na.rm = TRUE)
              row$hypo_meth_count <- sum(methylation_clean < 0.2, na.rm = TRUE)

              cat(sprintf(
                "  CpG sites: %d, Mean methylation: %.1f%%\n",
                row$cpg_sites, row$mean_methylation * 100
              ))

              # Store raw data for expanded plots
              if (expanded_plots) {
                all_meth_data[[sample_name]] <- methylation_clean
                all_chr_data[[sample_name]] <- data.frame(chr = chr_clean, methylation = methylation_clean)
                all_strand_data[[sample_name]] <- data.frame(strand = strand_clean, methylation = methylation_clean, coverage = coverage_clean)
              }
            }
          }
        }
      }
    }

    # Get qualimap statistics (all 27 fields)
    qualimap_file <- file.path(sample_dir, "qualimap", "genome_results.txt")
    qualimap_stats <- parse_qualimap(qualimap_file)
    for (stat_name in names(qualimap_stats)) {
      row[[stat_name]] <- qualimap_stats[[stat_name]]
    }

    if (!is.na(qualimap_stats$number_of_reads)) {
      cat(sprintf(
        "  Qualimap: %s reads, %.1f%% mapped\n",
        format(qualimap_stats$number_of_reads, big.mark = ","),
        qualimap_stats$mapped_reads_percentage
      ))
    }

    results_list[[length(results_list) + 1]] <- as.data.frame(row, stringsAsFactors = FALSE)
  }

  if (length(results_list) == 0) {
    cat("\nERROR: No valid data could be extracted\n")
    quit(status = 1)
  }

  results <- bind_rows(results_list)
}

# Save raw summary table
summary_file <- file.path(output_dir, "pipeline_summary.csv")
write.csv(results, summary_file, row.names = FALSE)
cat("\nRaw summary saved to:", summary_file, "\n")
cat("Columns:", ncol(results), "\n")

# Check for duplicates and create aggregated summary
agg_results <- aggregate_samples(results)
if (!is.null(agg_results)) {
  agg_file <- file.path(output_dir, "pipeline_summary_aggregated.csv")
  write.csv(agg_results, agg_file, row.names = FALSE)
  cat("Aggregated summary saved to:", agg_file, "\n")
  cat("Aggregated", nrow(results), "runs into", nrow(agg_results), "unique samples\n")
  plot_data <- agg_results
  cat("\n[INFO] Using aggregated data for all plots\n")
} else {
  plot_data <- results
}

# Print summary table
cat("\n==========================================\n")
cat("Summary Statistics\n")
cat("==========================================\n")
print_cols <- intersect(c(
  "sample", "n_runs", "cpg_sites", "mean_methylation", "mean_coverage",
  "number_of_reads", "mapped_reads_percentage"
), names(plot_data))
print(plot_data[, print_cols], row.names = FALSE)

# Sanity checks
cat("\n--- Sanity Checks ---\n")
if ("mean_methylation" %in% names(plot_data)) {
  cat(
    "Mean methylation range:", min(plot_data$mean_methylation, na.rm = TRUE), "-",
    max(plot_data$mean_methylation, na.rm = TRUE), "\n"
  )
}
if ("mean_coverage" %in% names(plot_data)) {
  cat(
    "Mean coverage range:", min(plot_data$mean_coverage, na.rm = TRUE), "-",
    max(plot_data$mean_coverage, na.rm = TRUE), "\n"
  )
}

# =============================================================================
# PLOTTING
# =============================================================================

plots_dir <- file.path(output_dir, "plots")
basic_dir <- file.path(plots_dir, "basic")
dir.create(basic_dir, recursive = TRUE, showWarnings = FALSE)

if (expanded_plots) {
  dist_dir <- file.path(plots_dir, "distribution")
  qc_dir <- file.path(plots_dir, "qc")
  bio_dir <- file.path(plots_dir, "biological")
  comp_dir <- file.path(plots_dir, "comparative")
  dir.create(dist_dir, showWarnings = FALSE)
  dir.create(qc_dir, showWarnings = FALSE)
  dir.create(bio_dir, showWarnings = FALSE)
  dir.create(comp_dir, showWarnings = FALSE)
}

theme_nanopore <- theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

cat("\n--- Generating Basic Plots ---\n")

# Plot 1: CpG sites
if ("cpg_sites" %in% names(plot_data) && any(!is.na(plot_data$cpg_sites) & plot_data$cpg_sites > 0)) {
  p1 <- ggplot(plot_data, aes(x = reorder(sample, -cpg_sites), y = cpg_sites)) +
    geom_bar(stat = "identity", fill = "#4e79a7", alpha = 0.8) +
    geom_text(aes(label = comma(cpg_sites)), vjust = -0.5, size = 3) +
    scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.1))) +
    labs(title = "CpG Sites per Sample", x = "Sample", y = "Number of CpG Sites") +
    theme_nanopore
  ggsave(file.path(basic_dir, "cpg_sites_per_sample.png"), p1, width = 10, height = 6, dpi = 300)
  ggsave(file.path(basic_dir, "cpg_sites_per_sample.pdf"), p1, width = 10, height = 6)
  cat("[OK] basic/cpg_sites_per_sample\n")
}

# Plot 2: Mean methylation
if ("mean_methylation" %in% names(plot_data) && any(!is.na(plot_data$mean_methylation))) {
  p2 <- ggplot(plot_data, aes(x = reorder(sample, -mean_methylation), y = mean_methylation * 100)) +
    geom_bar(stat = "identity", fill = "#e15759", alpha = 0.8) +
    geom_text(aes(label = sprintf("%.1f%%", mean_methylation * 100)), vjust = -0.5, size = 3) +
    scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.05))) +
    labs(title = "Mean CpG Methylation per Sample", x = "Sample", y = "Mean Methylation (%)") +
    theme_nanopore
  ggsave(file.path(basic_dir, "mean_methylation_per_sample.png"), p2, width = 10, height = 6, dpi = 300)
  ggsave(file.path(basic_dir, "mean_methylation_per_sample.pdf"), p2, width = 10, height = 6)
  cat("[OK] basic/mean_methylation_per_sample\n")
}

# Plot 3: Mean coverage
if ("mean_coverage" %in% names(plot_data) && any(!is.na(plot_data$mean_coverage) & plot_data$mean_coverage > 0)) {
  p3 <- ggplot(plot_data, aes(x = reorder(sample, -mean_coverage), y = mean_coverage)) +
    geom_bar(stat = "identity", fill = "#59a14f", alpha = 0.8) +
    geom_text(aes(label = sprintf("%.1f", mean_coverage)), vjust = -0.5, size = 3) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    labs(title = "Mean Coverage per Sample", x = "Sample", y = "Mean Coverage (reads)") +
    theme_nanopore
  ggsave(file.path(basic_dir, "mean_coverage_per_sample.png"), p3, width = 10, height = 6, dpi = 300)
  ggsave(file.path(basic_dir, "mean_coverage_per_sample.pdf"), p3, width = 10, height = 6)
  cat("[OK] basic/mean_coverage_per_sample\n")
}

# Plot 4: Total reads (from qualimap)
if ("number_of_reads" %in% names(plot_data) && !all(is.na(plot_data$number_of_reads))) {
  pd <- plot_data[!is.na(plot_data$number_of_reads), ]
  if (nrow(pd) > 0) {
    p4 <- ggplot(pd, aes(x = reorder(sample, -number_of_reads), y = number_of_reads)) +
      geom_bar(stat = "identity", fill = "#76b7b2", alpha = 0.8) +
      geom_text(aes(label = comma(number_of_reads)), vjust = -0.5, size = 3) +
      scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.1))) +
      labs(title = "Total Reads per Sample", x = "Sample", y = "Number of Reads") +
      theme_nanopore
    ggsave(file.path(basic_dir, "total_reads_per_sample.png"), p4, width = 10, height = 6, dpi = 300)
    ggsave(file.path(basic_dir, "total_reads_per_sample.pdf"), p4, width = 10, height = 6)
    cat("[OK] basic/total_reads_per_sample\n")
  }
}

# Plot 5: Mapped reads percentage
if ("mapped_reads_percentage" %in% names(plot_data) && !all(is.na(plot_data$mapped_reads_percentage))) {
  pd <- plot_data[!is.na(plot_data$mapped_reads_percentage), ]
  if (nrow(pd) > 0) {
    p5 <- ggplot(pd, aes(x = reorder(sample, -mapped_reads_percentage), y = mapped_reads_percentage)) +
      geom_bar(stat = "identity", fill = "#f28e2c", alpha = 0.8) +
      geom_text(aes(label = sprintf("%.1f%%", mapped_reads_percentage)), vjust = -0.5, size = 3) +
      scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.05))) +
      labs(title = "Mapped Reads Percentage", x = "Sample", y = "Mapped Reads (%)") +
      theme_nanopore
    ggsave(file.path(basic_dir, "mapped_reads_percentage.png"), p5, width = 10, height = 6, dpi = 300)
    ggsave(file.path(basic_dir, "mapped_reads_percentage.pdf"), p5, width = 10, height = 6)
    cat("[OK] basic/mapped_reads_percentage\n")
  }
}

# =============================================================================
# EXPANDED PLOTS
# =============================================================================

if (expanded_plots) {
  cat("\n--- Generating Expanded Plots ---\n")

  # QC: Mapping quality
  if ("mean_mapping_quality" %in% names(plot_data) && any(!is.na(plot_data$mean_mapping_quality))) {
    pq <- ggplot(plot_data, aes(x = reorder(sample, -mean_mapping_quality), y = mean_mapping_quality)) +
      geom_bar(stat = "identity", fill = "#9c755f", alpha = 0.8) +
      geom_text(aes(label = sprintf("%.1f", mean_mapping_quality)), vjust = -0.5, size = 3) +
      labs(title = "Mean Mapping Quality", x = "Sample", y = "Mapping Quality") +
      theme_nanopore
    ggsave(file.path(qc_dir, "mean_mapping_quality.png"), pq, width = 10, height = 6, dpi = 300)
    ggsave(file.path(qc_dir, "mean_mapping_quality.pdf"), pq, width = 10, height = 6)
    cat("[OK] qc/mean_mapping_quality\n")
  }

  # QC: Error rate
  if ("general_error_rate" %in% names(plot_data) && any(!is.na(plot_data$general_error_rate))) {
    pe <- ggplot(plot_data, aes(x = reorder(sample, general_error_rate), y = general_error_rate * 100)) +
      geom_bar(stat = "identity", fill = "#e15759", alpha = 0.8) +
      geom_text(aes(label = sprintf("%.2f%%", general_error_rate * 100)), hjust = -0.1, size = 3) +
      coord_flip() +
      labs(title = "General Error Rate", x = "Sample", y = "Error Rate (%)") +
      theme_nanopore +
      theme(axis.text.x = element_text(angle = 0))
    ggsave(file.path(qc_dir, "general_error_rate.png"), pe, width = 10, height = 6, dpi = 300)
    ggsave(file.path(qc_dir, "general_error_rate.pdf"), pe, width = 10, height = 6)
    cat("[OK] qc/general_error_rate\n")
  }

  # Biological: GC content
  if ("GC_percentage" %in% names(plot_data) && any(!is.na(plot_data$GC_percentage))) {
    pg <- ggplot(plot_data, aes(x = reorder(sample, -GC_percentage), y = GC_percentage)) +
      geom_bar(stat = "identity", fill = "#59a14f", alpha = 0.8) +
      geom_text(aes(label = sprintf("%.1f%%", GC_percentage)), vjust = -0.5, size = 3) +
      scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.05))) +
      labs(title = "GC Content", x = "Sample", y = "GC (%)") +
      theme_nanopore
    ggsave(file.path(bio_dir, "gc_content.png"), pg, width = 10, height = 6, dpi = 300)
    ggsave(file.path(bio_dir, "gc_content.pdf"), pg, width = 10, height = 6)
    cat("[OK] biological/gc_content\n")
  }

  # Biological: Hyper/Hypo methylated
  if (all(c("hyper_meth_count", "hypo_meth_count") %in% names(plot_data))) {
    hh <- plot_data %>%
      select(sample, hyper_meth_count, hypo_meth_count) %>%
      pivot_longer(cols = c(hyper_meth_count, hypo_meth_count), names_to = "type", values_to = "count") %>%
      mutate(type = ifelse(type == "hyper_meth_count", "Hyper (>80%)", "Hypo (<20%)"))

    ph <- ggplot(hh, aes(x = sample, y = count, fill = type)) +
      geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
      scale_fill_manual(values = c("Hyper (>80%)" = "#e15759", "Hypo (<20%)" = "#4e79a7")) +
      scale_y_continuous(labels = comma) +
      labs(title = "Hyper/Hypo-Methylated CpG", x = "Sample", y = "Count", fill = "State") +
      theme_nanopore
    ggsave(file.path(bio_dir, "hyper_hypo_methylated.png"), ph, width = 10, height = 6, dpi = 300)
    ggsave(file.path(bio_dir, "hyper_hypo_methylated.pdf"), ph, width = 10, height = 6)
    cat("[OK] biological/hyper_hypo_methylated\n")
  }
}

# Final summary
cat("\n==========================================\n")
cat("Plots saved to:", plots_dir, "\n")
cat("==========================================\n")
cat("--- Summary report complete!\n")
