#!/usr/bin/env Rscript

# NanoporeToBED Pipeline - Summary Statistics and Plots
# Version: 1.0.0
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
  
  # Read BED file (modkit pileup format)
  # Columns: chr, start, end, mod_code, score, strand, start2, end2, color, n_valid, fraction, n_mod, n_canon, n_other, n_delete, n_fail, n_diff, n_nocall
  bed_data <- tryCatch({
    read.table(bed_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  }, error = function(e) {
    cat("  Warning: Could not read", bed_file, "\n")
    return(NULL)
  })
  
  if (is.null(bed_data) || nrow(bed_data) == 0) {
    next
  }
  
  # Parse standard modkit pileup BED format
  # Column 11 is fraction (methylation level), column 10 is n_valid (coverage)
  if (ncol(bed_data) >= 11) {
    methylation <- as.numeric(bed_data[, 11])
    coverage <- as.numeric(bed_data[, 10])
  } else {
    cat("  Warning: Unexpected BED format\n")
    next
  }
  
  # Calculate statistics
  cpg_sites <- nrow(bed_data)
  mean_meth <- mean(methylation, na.rm = TRUE)
  median_meth <- median(methylation, na.rm = TRUE)
  mean_cov <- mean(coverage, na.rm = TRUE)
  
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

# Save summary table
summary_file <- file.path(output_dir, "pipeline_summary.csv")
write.csv(results, summary_file, row.names = FALSE)
cat("\nSummary saved to:", summary_file, "\n")

# Print summary table
cat("\n==========================================\n")
cat("Summary Statistics\n")
cat("==========================================\n")
print(results, row.names = FALSE)

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

# Plot 2: Mean methylation per sample
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

# Plot 3: Mean coverage per sample
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

# Plot 4: Total reads per sample (if available)
if (!all(is.na(results$total_reads))) {
  results_reads <- results[!is.na(results$total_reads), ]
  
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
}

# Plot 5: Combined summary plot
results_long <- results %>%
  select(sample, cpg_sites, mean_methylation, mean_coverage) %>%
  mutate(
    cpg_sites_scaled = cpg_sites / max(cpg_sites, na.rm = TRUE),
    methylation_scaled = mean_methylation,
    coverage_scaled = mean_coverage / max(mean_coverage, na.rm = TRUE)
  ) %>%
  select(sample, cpg_sites_scaled, methylation_scaled, coverage_scaled) %>%
  pivot_longer(cols = -sample, names_to = "metric", values_to = "value") %>%
  mutate(metric = case_when(
    metric == "cpg_sites_scaled" ~ "CpG Sites (scaled)",
    metric == "methylation_scaled" ~ "Methylation",
    metric == "coverage_scaled" ~ "Coverage (scaled)"
  ))

p5 <- ggplot(results_long, aes(x = sample, y = value, fill = metric)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
  scale_fill_manual(values = c("#4e79a7", "#e15759", "#59a14f")) +
  scale_y_continuous(labels = percent_format(), expand = expansion(mult = c(0, 0.05))) +
  labs(
    title = "Sample Quality Metrics Overview",
    x = "Sample",
    y = "Normalized Value",
    fill = "Metric"
  ) +
  theme_nanopore

ggsave(file.path(plots_dir, "summary_overview.png"), p5, width = 12, height = 6, dpi = 300)
ggsave(file.path(plots_dir, "summary_overview.pdf"), p5, width = 12, height = 6)

cat("\n==========================================\n")
cat("Plots saved to:", plots_dir, "\n")
cat("==========================================\n")
cat("Generated plots:\n")
cat("  - cpg_sites_per_sample.png/pdf\n")
cat("  - mean_methylation_per_sample.png/pdf\n")
cat("  - mean_coverage_per_sample.png/pdf\n")
if (!all(is.na(results$total_reads))) {
  cat("  - total_reads_per_sample.png/pdf\n")
}
cat("  - summary_overview.png/pdf\n")
cat("\n🎉 Summary report complete!\n")
