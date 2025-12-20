#!/usr/bin/env Rscript

# NanoporeToBED Pipeline - Summary Statistics and Plots
# Version: 1.1.0
# Author: Markus Hodal Drag

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  cat("Usage: Rscript generate_summary.R <output_dir>\n")
  quit(status = 1)
}

output_dir <- args[1]

cat("==========================================\n")
cat("NanoporeToBED Summary Report\n")
cat("==========================================\n")
cat("Output directory:", output_dir, "\n\n")

# Find all sample directories (those containing .CpG.bed files)
bed_files <- list.files(output_dir, pattern = "\\.CpG\\.bed$", recursive = TRUE, full.names = TRUE)

if (length(bed_files) == 0) {
  cat("ERROR: No .CpG.bed files found in", output_dir, "\n")
  quit(status = 1)
}

cat("Found", length(bed_files), "sample(s)\n\n")

# Initialize results dataframe
results <- data.frame(
  sample = character(),
  cpg_sites = numeric(),
  mean_methylation = numeric(),
  median_methylation = numeric(),
  mean_coverage = numeric(),
  total_reads = numeric(),
  stringsAsFactors = FALSE
)

# Process each sample
for (bed_file in bed_files) {
  sample_name <- gsub("\\.CpG\\.bed$", "", basename(bed_file))
  sample_dir <- dirname(bed_file)

  cat("Processing:", sample_name, "\n")

  # Read BED file - modkit pileup bedMethyl format
  # Uses mixed delimiters by default (tabs for 1-9, spaces for rest)
  # Read raw lines and parse manually for robustness
  bed_lines <- tryCatch(
    {
      readLines(bed_file, warn = FALSE)
    },
    error = function(e) {
      cat("  Warning: Could not read", bed_file, "\n")
      return(NULL)
    }
  )

  if (is.null(bed_lines) || length(bed_lines) == 0) {
    cat("  Warning: Empty BED file\n")
    next
  }

  # Remove any header lines (start with #) or empty lines
  bed_lines <- bed_lines[!grepl("^#", bed_lines) & nchar(bed_lines) > 0]

  if (length(bed_lines) == 0) {
    cat("  Warning: No data lines in BED file\n")
    next
  }

  # Parse bedMethyl format: split on whitespace (tabs and spaces)
  # Standard bedMethyl columns:
  # 0:chr 1:start 2:end 3:mod_code 4:score 5:strand 6:start2 7:end2 8:color
  # 9:Nvalid_cov 10:percent_modified 11:Nmod 12:Ncanonical 13:Nother 14:Ndel 15:Nfail 16:Ndiff 17:Nnocall
  # Column 10 is the percent methylation (0-100) or fraction (0-1) depending on modkit version

  bed_data <- tryCatch(
    {
      data <- read.table(text = bed_lines, header = FALSE, stringsAsFactors = FALSE, fill = TRUE)
      data
    },
    error = function(e) {
      cat("  Warning: Could not parse BED data:", conditionMessage(e), "\n")
      return(NULL)
    }
  )

  if (is.null(bed_data) || nrow(bed_data) == 0) {
    next
  }

  # Debug: print column count and sample of data
  cat("  BED columns:", ncol(bed_data), "\n")

  # Extract methylation and coverage based on available columns
  if (ncol(bed_data) >= 12) {
    # Standard bedMethyl format
    # Column 10 (index 10) = Nvalid_cov (coverage)
    # Column 11 (index 11) = percent modified (can be 0-100 or 0-1)
    # Column 12 (index 12) = Nmod
    # Column 13 (index 13) = Ncanonical
    coverage <- as.numeric(bed_data[, 10])
    percent_mod <- as.numeric(bed_data[, 11])

    # Check if percent_mod is 0-100 or 0-1 scale
    if (max(percent_mod, na.rm = TRUE) > 1) {
      # It's 0-100 scale, convert to 0-1
      methylation <- percent_mod / 100
      cat("  Methylation scale: 0-100 (converted to fraction)\n")
    } else {
      # Already 0-1 scale
      methylation <- percent_mod
      cat("  Methylation scale: 0-1 (fraction)\n")
    }
  } else if (ncol(bed_data) >= 5) {
    # Minimal BED format - use score column (column 5)
    # Some tools put methylation percentage in score
    score <- as.numeric(bed_data[, 5])
    if (max(score, na.rm = TRUE) > 1) {
      methylation <- score / 100
    } else {
      methylation <- score
    }
    coverage <- rep(NA, nrow(bed_data))
    cat("  Using minimal BED format (score column)\n")
  } else {
    cat("  Warning: Unexpected BED format (only", ncol(bed_data), "columns)\n")
    next
  }

  # Sanity checks on methylation values
  valid_meth <- !is.na(methylation) & methylation >= 0 & methylation <= 1
  if (sum(valid_meth) == 0) {
    cat("  Warning: No valid methylation values found\n")
    cat("  Sample values:", head(methylation, 5), "\n")
    next
  }

  methylation_clean <- methylation[valid_meth]
  coverage_clean <- coverage[valid_meth]

  # Calculate statistics
  cpg_sites <- length(methylation_clean)
  mean_meth <- mean(methylation_clean, na.rm = TRUE)
  median_meth <- median(methylation_clean, na.rm = TRUE)
  mean_cov <- mean(coverage_clean, na.rm = TRUE)

  # Sanity check: print statistics
  cat(sprintf("  CpG sites: %d\n", cpg_sites))
  cat(sprintf("  Mean methylation: %.3f (%.1f%%)\n", mean_meth, mean_meth * 100))
  cat(sprintf("  Mean coverage: %.1f\n", mean_cov))

  # Try to get read count from qualimap
  qualimap_file <- file.path(sample_dir, "qualimap", "genome_results.txt")
  total_reads <- NA
  if (file.exists(qualimap_file)) {
    qualimap_lines <- readLines(qualimap_file)
    reads_line <- grep("number of reads", qualimap_lines, value = TRUE, ignore.case = TRUE)
    if (length(reads_line) > 0) {
      total_reads <- as.numeric(gsub("[^0-9]", "", reads_line[1]))
    }
  }

  # Add to results
  results <- rbind(results, data.frame(
    sample = sample_name,
    cpg_sites = cpg_sites,
    mean_methylation = mean_meth,
    median_methylation = median_meth,
    mean_coverage = mean_cov,
    total_reads = total_reads,
    stringsAsFactors = FALSE
  ))
}

# Check if we have any valid results
if (nrow(results) == 0) {
  cat("\nERROR: No valid data could be extracted from BED files\n")
  quit(status = 1)
}

# Save summary table
summary_file <- file.path(output_dir, "pipeline_summary.csv")
write.csv(results, summary_file, row.names = FALSE)
cat("\nSummary saved to:", summary_file, "\n")

# Print summary table
cat("\n==========================================\n")
cat("Summary Statistics\n")
cat("==========================================\n")
print(results, row.names = FALSE)

# Sanity check: verify all numeric columns have valid values
cat("\n--- Sanity Checks ---\n")
cat(
  "Mean methylation range:", min(results$mean_methylation, na.rm = TRUE), "-",
  max(results$mean_methylation, na.rm = TRUE), "\n"
)
cat(
  "Mean coverage range:", min(results$mean_coverage, na.rm = TRUE), "-",
  max(results$mean_coverage, na.rm = TRUE), "\n"
)
cat(
  "CpG sites range:", min(results$cpg_sites, na.rm = TRUE), "-",
  max(results$cpg_sites, na.rm = TRUE), "\n"
)

if (all(is.na(results$mean_methylation)) || all(results$mean_methylation == 0)) {
  cat("WARNING: All methylation values are NA or 0 - check BED file format!\n")
}

# Create plots directory
plots_dir <- file.path(output_dir, "plots")
dir.create(plots_dir, showWarnings = FALSE)

# Theme for plots
theme_nanopore <- theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

# Plot 1: CpG sites per sample
if (any(!is.na(results$cpg_sites) & results$cpg_sites > 0)) {
  p1 <- ggplot(results, aes(x = reorder(sample, -cpg_sites), y = cpg_sites)) +
    geom_bar(stat = "identity", fill = "#4e79a7", alpha = 0.8) +
    geom_text(aes(label = comma(cpg_sites)), vjust = -0.5, size = 3) +
    scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.1))) +
    labs(
      title = "CpG Sites per Sample",
      x = "Sample",
      y = "Number of CpG Sites"
    ) +
    theme_nanopore

  ggsave(file.path(plots_dir, "cpg_sites_per_sample.png"), p1, width = 10, height = 6, dpi = 300)
  ggsave(file.path(plots_dir, "cpg_sites_per_sample.pdf"), p1, width = 10, height = 6)
  cat("✓ Created: cpg_sites_per_sample.png/pdf\n")
} else {
  cat("⚠ Skipped CpG sites plot (no valid data)\n")
}

# Plot 2: Mean methylation per sample
if (any(!is.na(results$mean_methylation) & results$mean_methylation > 0)) {
  p2 <- ggplot(results, aes(x = reorder(sample, -mean_methylation), y = mean_methylation * 100)) +
    geom_bar(stat = "identity", fill = "#e15759", alpha = 0.8) +
    geom_text(aes(label = sprintf("%.1f%%", mean_methylation * 100)), vjust = -0.5, size = 3) +
    scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.05))) +
    labs(
      title = "Mean CpG Methylation per Sample",
      x = "Sample",
      y = "Mean Methylation (%)"
    ) +
    theme_nanopore

  ggsave(file.path(plots_dir, "mean_methylation_per_sample.png"), p2, width = 10, height = 6, dpi = 300)
  ggsave(file.path(plots_dir, "mean_methylation_per_sample.pdf"), p2, width = 10, height = 6)
  cat("✓ Created: mean_methylation_per_sample.png/pdf\n")
} else {
  cat("⚠ Skipped methylation plot (no valid data)\n")
}

# Plot 3: Mean coverage per sample
if (any(!is.na(results$mean_coverage) & results$mean_coverage > 0)) {
  p3 <- ggplot(results, aes(x = reorder(sample, -mean_coverage), y = mean_coverage)) +
    geom_bar(stat = "identity", fill = "#59a14f", alpha = 0.8) +
    geom_text(aes(label = sprintf("%.1f", mean_coverage)), vjust = -0.5, size = 3) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    labs(
      title = "Mean Coverage per Sample",
      x = "Sample",
      y = "Mean Coverage (reads)"
    ) +
    theme_nanopore

  ggsave(file.path(plots_dir, "mean_coverage_per_sample.png"), p3, width = 10, height = 6, dpi = 300)
  ggsave(file.path(plots_dir, "mean_coverage_per_sample.pdf"), p3, width = 10, height = 6)
  cat("✓ Created: mean_coverage_per_sample.png/pdf\n")
} else {
  cat("⚠ Skipped coverage plot (no valid data)\n")
}

# Plot 4: Total reads per sample (if available)
if (!all(is.na(results$total_reads))) {
  results_reads <- results[!is.na(results$total_reads), ]

  if (nrow(results_reads) > 0) {
    p4 <- ggplot(results_reads, aes(x = reorder(sample, -total_reads), y = total_reads)) +
      geom_bar(stat = "identity", fill = "#76b7b2", alpha = 0.8) +
      geom_text(aes(label = comma(total_reads)), vjust = -0.5, size = 3) +
      scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.1))) +
      labs(
        title = "Total Reads per Sample",
        x = "Sample",
        y = "Number of Reads"
      ) +
      theme_nanopore

    ggsave(file.path(plots_dir, "total_reads_per_sample.png"), p4, width = 10, height = 6, dpi = 300)
    ggsave(file.path(plots_dir, "total_reads_per_sample.pdf"), p4, width = 10, height = 6)
    cat("✓ Created: total_reads_per_sample.png/pdf\n")
  }
}

# Plot 5: Individual metric panels (NOT combined - each on own scale)
# This replaces the problematic "overview" that mixed incompatible scales

# Create a faceted plot with free scales for each metric
results_long <- results %>%
  select(sample, cpg_sites, mean_methylation, mean_coverage) %>%
  mutate(
    `CpG Sites` = cpg_sites,
    `Methylation (%)` = mean_methylation * 100,
    `Mean Coverage` = mean_coverage
  ) %>%
  select(sample, `CpG Sites`, `Methylation (%)`, `Mean Coverage`) %>%
  pivot_longer(cols = -sample, names_to = "Metric", values_to = "Value")

# Check if we have valid data for the faceted plot
if (any(!is.na(results_long$Value))) {
  p5 <- ggplot(results_long, aes(x = sample, y = Value, fill = Metric)) +
    geom_bar(stat = "identity", alpha = 0.8) +
    facet_wrap(~Metric, scales = "free_y", ncol = 1) +
    scale_fill_manual(values = c(
      "CpG Sites" = "#4e79a7",
      "Methylation (%)" = "#e15759",
      "Mean Coverage" = "#59a14f"
    )) +
    labs(
      title = "Sample Quality Metrics Overview",
      x = "Sample",
      y = "Value"
    ) +
    theme_nanopore +
    theme(
      legend.position = "none",
      strip.text = element_text(size = 11, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )

  ggsave(file.path(plots_dir, "summary_overview.png"), p5, width = 10, height = 10, dpi = 300)
  ggsave(file.path(plots_dir, "summary_overview.pdf"), p5, width = 10, height = 10)
  cat("✓ Created: summary_overview.png/pdf\n")
} else {
  cat("⚠ Skipped overview plot (no valid data)\n")
}

# Plot 6: Methylation distribution (if we have raw data - optional box/violin plot)
# This would require re-reading the BED files, so skip for now

cat("\n==========================================\n")
cat("Plots saved to:", plots_dir, "\n")
cat("==========================================\n")
cat("\n🎉 Summary report complete!\n")
