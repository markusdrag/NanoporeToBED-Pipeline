#!/usr/bin/env Rscript

# NanoporeToBED Pipeline - Summary Statistics and Plots
# Version: 1.4.1
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
  cat("Usage: Rscript generate_summary.R <output_dir> [--expanded-plots]\n")
  cat("  --expanded-plots  Generate additional distribution, QC, and comparative plots\n")
  quit(status = 1)
}

output_dir <- args[1]
expanded_plots <- "--expanded-plots" %in% args

cat("==========================================\n")
cat("NanoporeToBED Summary Report v1.4.0\n")
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
  low_cov_percent = numeric(),
  hyper_meth_count = numeric(),
  hypo_meth_count = numeric(),
  stringsAsFactors = FALSE
)

# Store raw data for expanded plots
all_meth_data <- list()
all_chr_data <- list()
all_strand_data <- list()

# Process each sample
for (bed_file in bed_files) {
  sample_name <- gsub("\\.CpG\\.bed$", "", basename(bed_file))
  sample_dir <- dirname(bed_file)

  cat("Processing:", sample_name, "\n")

  # Read BED file
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

  bed_lines <- bed_lines[!grepl("^#", bed_lines) & nchar(bed_lines) > 0]

  if (length(bed_lines) == 0) {
    cat("  Warning: No data lines in BED file\n")
    next
  }

  bed_data <- tryCatch(
    {
      read.table(text = bed_lines, header = FALSE, stringsAsFactors = FALSE, fill = TRUE)
    },
    error = function(e) {
      cat("  Warning: Could not parse BED data:", conditionMessage(e), "\n")
      return(NULL)
    }
  )

  if (is.null(bed_data) || nrow(bed_data) == 0) next

  cat("  BED columns:", ncol(bed_data), "\n")

  # Extract data based on columns available
  if (ncol(bed_data) >= 12) {
    chr <- bed_data[, 1]
    strand <- bed_data[, 6]
    coverage <- as.numeric(bed_data[, 10])
    percent_mod <- as.numeric(bed_data[, 11])

    if (max(percent_mod, na.rm = TRUE) > 1) {
      methylation <- percent_mod / 100
      cat("  Methylation scale: 0-100 (converted to fraction)\n")
    } else {
      methylation <- percent_mod
      cat("  Methylation scale: 0-1 (fraction)\n")
    }
  } else if (ncol(bed_data) >= 6) {
    chr <- bed_data[, 1]
    strand <- bed_data[, 6]
    score <- as.numeric(bed_data[, 5])
    if (max(score, na.rm = TRUE) > 1) {
      methylation <- score / 100
    } else {
      methylation <- score
    }
    coverage <- rep(NA, nrow(bed_data))
    cat("  Using minimal BED format\n")
  } else {
    cat("  Warning: Unexpected BED format (only", ncol(bed_data), "columns)\n")
    next
  }

  # Clean data
  valid_meth <- !is.na(methylation) & methylation >= 0 & methylation <= 1
  if (sum(valid_meth) == 0) {
    cat("  Warning: No valid methylation values found\n")
    next
  }

  methylation_clean <- methylation[valid_meth]
  coverage_clean <- coverage[valid_meth]
  chr_clean <- chr[valid_meth]
  strand_clean <- strand[valid_meth]

  # Calculate statistics
  cpg_sites <- length(methylation_clean)
  mean_meth <- mean(methylation_clean, na.rm = TRUE)
  median_meth <- median(methylation_clean, na.rm = TRUE)
  mean_cov <- mean(coverage_clean, na.rm = TRUE)

  # Extended stats for expanded plots
  low_cov_percent <- sum(coverage_clean < 10, na.rm = TRUE) / cpg_sites * 100
  hyper_meth_count <- sum(methylation_clean > 0.8, na.rm = TRUE)
  hypo_meth_count <- sum(methylation_clean < 0.2, na.rm = TRUE)

  cat(sprintf("  CpG sites: %d\n", cpg_sites))
  cat(sprintf("  Mean methylation: %.3f (%.1f%%)\n", mean_meth, mean_meth * 100))
  cat(sprintf("  Mean coverage: %.1f\n", mean_cov))

  # Get total reads from qualimap
  qualimap_file <- file.path(sample_dir, "qualimap", "genome_results.txt")
  total_reads <- NA
  if (file.exists(qualimap_file)) {
    qualimap_lines <- readLines(qualimap_file)
    reads_line <- grep("number of reads", qualimap_lines, value = TRUE, ignore.case = TRUE)
    if (length(reads_line) > 0) {
      total_reads <- as.numeric(gsub("[^0-9]", "", reads_line[1]))
    }
  }

  # Store raw data for expanded plots
  if (expanded_plots) {
    all_meth_data[[sample_name]] <- methylation_clean
    all_chr_data[[sample_name]] <- data.frame(
      chr = chr_clean,
      methylation = methylation_clean,
      stringsAsFactors = FALSE
    )
    all_strand_data[[sample_name]] <- data.frame(
      strand = strand_clean,
      methylation = methylation_clean,
      coverage = coverage_clean,
      stringsAsFactors = FALSE
    )
  }

  # Add to results
  results <- rbind(results, data.frame(
    sample = sample_name,
    cpg_sites = cpg_sites,
    mean_methylation = mean_meth,
    median_methylation = median_meth,
    mean_coverage = mean_cov,
    total_reads = total_reads,
    low_cov_percent = low_cov_percent,
    hyper_meth_count = hyper_meth_count,
    hypo_meth_count = hypo_meth_count,
    stringsAsFactors = FALSE
  ))
}

# Check results
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
print(results[, c("sample", "cpg_sites", "mean_methylation", "mean_coverage", "total_reads")], row.names = FALSE)

# Sanity checks
cat("\n--- Sanity Checks ---\n")
cat(
  "Mean methylation range:", min(results$mean_methylation, na.rm = TRUE), "-",
  max(results$mean_methylation, na.rm = TRUE), "\n"
)
cat(
  "Mean coverage range:", min(results$mean_coverage, na.rm = TRUE), "-",
  max(results$mean_coverage, na.rm = TRUE), "\n"
)

# Create directory structure
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

# Theme for plots
theme_nanopore <- theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

cat("\n--- Generating Basic Plots ---\n")

# ============================================================================
# BASIC PLOTS (always generated)
# ============================================================================

# Plot 1: CpG sites per sample
if (any(!is.na(results$cpg_sites) & results$cpg_sites > 0)) {
  p1 <- ggplot(results, aes(x = reorder(sample, -cpg_sites), y = cpg_sites)) +
    geom_bar(stat = "identity", fill = "#4e79a7", alpha = 0.8) +
    geom_text(aes(label = comma(cpg_sites)), vjust = -0.5, size = 3) +
    scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.1))) +
    labs(title = "CpG Sites per Sample", x = "Sample", y = "Number of CpG Sites") +
    theme_nanopore

  ggsave(file.path(basic_dir, "cpg_sites_per_sample.png"), p1, width = 10, height = 6, dpi = 300)
  ggsave(file.path(basic_dir, "cpg_sites_per_sample.pdf"), p1, width = 10, height = 6)
  cat("[OK] basic/cpg_sites_per_sample.png/pdf\n")
}

# Plot 2: Mean methylation per sample
if (any(!is.na(results$mean_methylation) & results$mean_methylation > 0)) {
  p2 <- ggplot(results, aes(x = reorder(sample, -mean_methylation), y = mean_methylation * 100)) +
    geom_bar(stat = "identity", fill = "#e15759", alpha = 0.8) +
    geom_text(aes(label = sprintf("%.1f%%", mean_methylation * 100)), vjust = -0.5, size = 3) +
    scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.05))) +
    labs(title = "Mean CpG Methylation per Sample", x = "Sample", y = "Mean Methylation (%)") +
    theme_nanopore

  ggsave(file.path(basic_dir, "mean_methylation_per_sample.png"), p2, width = 10, height = 6, dpi = 300)
  ggsave(file.path(basic_dir, "mean_methylation_per_sample.pdf"), p2, width = 10, height = 6)
  cat("[OK] basic/mean_methylation_per_sample.png/pdf\n")
}

# Plot 3: Mean coverage per sample
if (any(!is.na(results$mean_coverage) & results$mean_coverage > 0)) {
  p3 <- ggplot(results, aes(x = reorder(sample, -mean_coverage), y = mean_coverage)) +
    geom_bar(stat = "identity", fill = "#59a14f", alpha = 0.8) +
    geom_text(aes(label = sprintf("%.1f", mean_coverage)), vjust = -0.5, size = 3) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    labs(title = "Mean Coverage per Sample", x = "Sample", y = "Mean Coverage (reads)") +
    theme_nanopore

  ggsave(file.path(basic_dir, "mean_coverage_per_sample.png"), p3, width = 10, height = 6, dpi = 300)
  ggsave(file.path(basic_dir, "mean_coverage_per_sample.pdf"), p3, width = 10, height = 6)
  cat("[OK] basic/mean_coverage_per_sample.png/pdf\n")
}

# Plot 4: Total reads per sample
if (!all(is.na(results$total_reads))) {
  results_reads <- results[!is.na(results$total_reads), ]
  if (nrow(results_reads) > 0) {
    p4 <- ggplot(results_reads, aes(x = reorder(sample, -total_reads), y = total_reads)) +
      geom_bar(stat = "identity", fill = "#76b7b2", alpha = 0.8) +
      geom_text(aes(label = comma(total_reads)), vjust = -0.5, size = 3) +
      scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.1))) +
      labs(title = "Total Reads per Sample", x = "Sample", y = "Number of Reads") +
      theme_nanopore

    ggsave(file.path(basic_dir, "total_reads_per_sample.png"), p4, width = 10, height = 6, dpi = 300)
    ggsave(file.path(basic_dir, "total_reads_per_sample.pdf"), p4, width = 10, height = 6)
    cat("[OK] basic/total_reads_per_sample.png/pdf\n")
  }
}

# Plot 5: Summary overview (faceted)
results_long <- results %>%
  select(sample, cpg_sites, mean_methylation, total_reads) %>%
  mutate(`CpG Sites` = cpg_sites, `Methylation (%)` = mean_methylation * 100, `Total Reads` = total_reads) %>%
  select(sample, `CpG Sites`, `Methylation (%)`, `Total Reads`) %>%
  pivot_longer(cols = -sample, names_to = "Metric", values_to = "Value")

if (any(!is.na(results_long$Value))) {
  p5 <- ggplot(results_long, aes(x = sample, y = Value, fill = Metric)) +
    geom_bar(stat = "identity", alpha = 0.8) +
    facet_wrap(~Metric, scales = "free_y", ncol = 1) +
    scale_fill_manual(values = c("CpG Sites" = "#4e79a7", "Methylation (%)" = "#e15759", "Total Reads" = "#76b7b2")) +
    labs(title = "Sample Quality Metrics Overview", x = "Sample", y = "Value") +
    theme_nanopore +
    theme(legend.position = "none", strip.text = element_text(size = 11, face = "bold"))

  ggsave(file.path(basic_dir, "summary_overview.png"), p5, width = 10, height = 10, dpi = 300)
  ggsave(file.path(basic_dir, "summary_overview.pdf"), p5, width = 10, height = 10)
  cat("[OK] basic/summary_overview.png/pdf\n")
}

# ============================================================================
# EXPANDED PLOTS (only when --expanded-plots flag is used)
# ============================================================================

if (expanded_plots) {
  cat("\n--- Generating Expanded Plots ---\n")

  # Distribution plots
  cat("\n[Distribution Plots]\n")

  # 1. Methylation distribution violin plot
  if (length(all_meth_data) > 0) {
    meth_dist_df <- do.call(rbind, lapply(names(all_meth_data), function(s) {
      data.frame(sample = s, methylation = all_meth_data[[s]], stringsAsFactors = FALSE)
    }))

    # Subsample for plotting if too large
    if (nrow(meth_dist_df) > 100000) {
      set.seed(42)
      meth_dist_df <- meth_dist_df[sample(nrow(meth_dist_df), 100000), ]
    }

    p_dist1 <- ggplot(meth_dist_df, aes(x = sample, y = methylation * 100, fill = sample)) +
      geom_violin(alpha = 0.7, scale = "width") +
      geom_boxplot(width = 0.1, outlier.size = 0.5) +
      scale_y_continuous(limits = c(0, 100)) +
      labs(title = "Methylation Distribution per Sample", x = "Sample", y = "Methylation (%)") +
      theme_nanopore +
      theme(legend.position = "none")

    ggsave(file.path(dist_dir, "methylation_distribution.png"), p_dist1, width = 12, height = 6, dpi = 300)
    ggsave(file.path(dist_dir, "methylation_distribution.pdf"), p_dist1, width = 12, height = 6)
    cat("[OK] distribution/methylation_distribution.png/pdf\n")
  }

  # 2. Coverage distribution histogram
  if (length(all_strand_data) > 0 && any(sapply(all_strand_data, function(x) any(!is.na(x$coverage))))) {
    cov_dist_df <- do.call(rbind, lapply(names(all_strand_data), function(s) {
      data.frame(sample = s, coverage = all_strand_data[[s]]$coverage, stringsAsFactors = FALSE)
    }))
    cov_dist_df <- cov_dist_df[!is.na(cov_dist_df$coverage), ]

    if (nrow(cov_dist_df) > 100000) {
      set.seed(42)
      cov_dist_df <- cov_dist_df[sample(nrow(cov_dist_df), 100000), ]
    }

    p_dist2 <- ggplot(cov_dist_df, aes(x = coverage, fill = sample)) +
      geom_histogram(bins = 50, alpha = 0.7, position = "identity") +
      facet_wrap(~sample, scales = "free_y") +
      scale_x_continuous(limits = c(0, quantile(cov_dist_df$coverage, 0.99, na.rm = TRUE))) +
      labs(title = "Coverage Distribution per Sample", x = "Coverage (reads)", y = "Count") +
      theme_nanopore +
      theme(legend.position = "none")

    ggsave(file.path(dist_dir, "coverage_distribution.png"), p_dist2, width = 12, height = 8, dpi = 300)
    ggsave(file.path(dist_dir, "coverage_distribution.pdf"), p_dist2, width = 12, height = 8)
    cat("[OK] distribution/coverage_distribution.png/pdf\n")
  }

  # QC plots
  cat("\n[QC Plots]\n")

  # 3. Low-coverage CpG percentage
  if (any(!is.na(results$low_cov_percent))) {
    p_qc1 <- ggplot(results, aes(x = reorder(sample, low_cov_percent), y = low_cov_percent)) +
      geom_bar(stat = "identity", fill = "#f28e2c", alpha = 0.8) +
      geom_text(aes(label = sprintf("%.1f%%", low_cov_percent)), hjust = -0.1, size = 3) +
      coord_flip() +
      scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
      labs(title = "Low Coverage CpG Sites (<10x)", x = "Sample", y = "Percentage of CpG Sites") +
      theme_nanopore +
      theme(axis.text.x = element_text(angle = 0))

    ggsave(file.path(qc_dir, "low_coverage_cpg_percent.png"), p_qc1, width = 10, height = 6, dpi = 300)
    ggsave(file.path(qc_dir, "low_coverage_cpg_percent.pdf"), p_qc1, width = 10, height = 6)
    cat("[OK] qc/low_coverage_cpg_percent.png/pdf\n")
  }

  # 4. Strand bias plot
  if (length(all_strand_data) > 0) {
    strand_summary <- do.call(rbind, lapply(names(all_strand_data), function(s) {
      d <- all_strand_data[[s]]
      data.frame(
        sample = s,
        strand = c("+", "-"),
        mean_meth = c(
          mean(d$methylation[d$strand == "+"], na.rm = TRUE),
          mean(d$methylation[d$strand == "-"], na.rm = TRUE)
        ),
        stringsAsFactors = FALSE
      )
    }))

    p_qc2 <- ggplot(strand_summary, aes(x = sample, y = mean_meth * 100, fill = strand)) +
      geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
      scale_fill_manual(values = c("+" = "#4e79a7", "-" = "#e15759")) +
      labs(title = "Strand Bias: Mean Methylation by Strand", x = "Sample", y = "Mean Methylation (%)", fill = "Strand") +
      theme_nanopore

    ggsave(file.path(qc_dir, "strand_bias.png"), p_qc2, width = 10, height = 6, dpi = 300)
    ggsave(file.path(qc_dir, "strand_bias.pdf"), p_qc2, width = 10, height = 6)
    cat("[OK] qc/strand_bias.png/pdf\n")
  }

  # Biological plots
  cat("\n[Biological Plots]\n")

  # 5. Hyper/Hypo methylated counts
  hyper_hypo <- results %>%
    select(sample, hyper_meth_count, hypo_meth_count) %>%
    pivot_longer(cols = c(hyper_meth_count, hypo_meth_count), names_to = "type", values_to = "count") %>%
    mutate(type = ifelse(type == "hyper_meth_count", "Hyper (>80%)", "Hypo (<20%)"))

  p_bio1 <- ggplot(hyper_hypo, aes(x = sample, y = count, fill = type)) +
    geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
    scale_fill_manual(values = c("Hyper (>80%)" = "#e15759", "Hypo (<20%)" = "#4e79a7")) +
    scale_y_continuous(labels = comma) +
    labs(title = "Hyper/Hypo-Methylated CpG Counts", x = "Sample", y = "Number of CpG Sites", fill = "Methylation State") +
    theme_nanopore

  ggsave(file.path(bio_dir, "hyper_hypo_methylated_counts.png"), p_bio1, width = 10, height = 6, dpi = 300)
  ggsave(file.path(bio_dir, "hyper_hypo_methylated_counts.pdf"), p_bio1, width = 10, height = 6)
  cat("[OK] biological/hyper_hypo_methylated_counts.png/pdf\n")

  # 6. Methylation by chromosome
  if (length(all_chr_data) > 0) {
    chr_df <- do.call(rbind, lapply(names(all_chr_data), function(s) {
      d <- all_chr_data[[s]]
      d$sample <- s
      d
    }))

    # Filter to main chromosomes only
    main_chr <- paste0("chr", c(1:22, "X", "Y", "M"))
    chr_df_main <- chr_df[chr_df$chr %in% main_chr, ]

    if (nrow(chr_df_main) > 200000) {
      set.seed(42)
      chr_df_main <- chr_df_main[sample(nrow(chr_df_main), 200000), ]
    }

    if (nrow(chr_df_main) > 0) {
      chr_df_main$chr <- factor(chr_df_main$chr, levels = main_chr)

      p_bio2 <- ggplot(chr_df_main, aes(x = chr, y = methylation * 100)) +
        geom_boxplot(aes(fill = chr), outlier.size = 0.3, alpha = 0.7) +
        facet_wrap(~sample, ncol = 2) +
        labs(title = "Methylation by Chromosome", x = "Chromosome", y = "Methylation (%)") +
        theme_nanopore +
        theme(legend.position = "none", axis.text.x = element_text(angle = 90, vjust = 0.5, size = 7))

      ggsave(file.path(bio_dir, "methylation_by_chromosome.png"), p_bio2, width = 14, height = max(6, length(unique(chr_df_main$sample)) * 2), dpi = 300)
      ggsave(file.path(bio_dir, "methylation_by_chromosome.pdf"), p_bio2, width = 14, height = max(6, length(unique(chr_df_main$sample)) * 2))
      cat("[OK] biological/methylation_by_chromosome.png/pdf\n")
    }
  }

  # Comparative plots
  cat("\n[Comparative Plots]\n")

  # 7. Sample correlation heatmap
  if (length(all_meth_data) >= 2) {
    tryCatch(
      {
        # Create a matrix of methylation at shared positions (simplified: use random sample)
        n_samples <- length(all_meth_data)
        cor_matrix <- matrix(NA, n_samples, n_samples)
        rownames(cor_matrix) <- names(all_meth_data)
        colnames(cor_matrix) <- names(all_meth_data)

        # Use sample means for correlation (a simplified approach)
        for (i in 1:n_samples) {
          for (j in 1:n_samples) {
            if (i == j) {
              cor_matrix[i, j] <- 1
            } else {
              # Subsample and correlate
              n_min <- min(length(all_meth_data[[i]]), length(all_meth_data[[j]]), 10000)
              set.seed(i * j)
              cor_matrix[i, j] <- cor(
                sample(all_meth_data[[i]], n_min),
                sample(all_meth_data[[j]], n_min),
                use = "complete.obs"
              )
            }
          }
        }

        cor_df <- as.data.frame(as.table(cor_matrix))
        names(cor_df) <- c("Sample1", "Sample2", "Correlation")

        p_comp1 <- ggplot(cor_df, aes(x = Sample1, y = Sample2, fill = Correlation)) +
          geom_tile() +
          geom_text(aes(label = sprintf("%.2f", Correlation)), size = 3) +
          scale_fill_gradient2(low = "#4e79a7", mid = "white", high = "#e15759", midpoint = 0.5, limits = c(0, 1)) +
          labs(title = "Sample Correlation Heatmap (Methylation)", x = "", y = "") +
          theme_nanopore +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))

        ggsave(file.path(comp_dir, "sample_correlation_heatmap.png"), p_comp1, width = 8, height = 7, dpi = 300)
        ggsave(file.path(comp_dir, "sample_correlation_heatmap.pdf"), p_comp1, width = 8, height = 7)
        cat("[OK] comparative/sample_correlation_heatmap.png/pdf\n")
      },
      error = function(e) {
        cat("[WARN] Skipped correlation heatmap:", conditionMessage(e), "\n")
      }
    )
  }

  # 8. PCA plot
  if (length(all_meth_data) >= 3) {
    tryCatch(
      {
        # Create a simplified PCA using summary statistics
        pca_data <- data.frame(
          sample = results$sample,
          mean_meth = results$mean_methylation,
          median_meth = results$median_methylation,
          mean_cov = results$mean_coverage,
          cpg_sites = results$cpg_sites,
          hyper = results$hyper_meth_count,
          hypo = results$hypo_meth_count
        )
        pca_data <- pca_data[complete.cases(pca_data), ]

        if (nrow(pca_data) >= 3) {
          pca_result <- prcomp(pca_data[, -1], scale. = TRUE)
          pca_coords <- as.data.frame(pca_result$x[, 1:2])
          pca_coords$sample <- pca_data$sample
          var_explained <- round(100 * summary(pca_result)$importance[2, 1:2], 1)

          p_comp2 <- ggplot(pca_coords, aes(x = PC1, y = PC2, label = sample)) +
            geom_point(size = 4, color = "#4e79a7", alpha = 0.8) +
            geom_text(vjust = -1, size = 3) +
            labs(
              title = "PCA of Sample Methylation Metrics",
              x = paste0("PC1 (", var_explained[1], "% variance)"),
              y = paste0("PC2 (", var_explained[2], "% variance)")
            ) +
            theme_nanopore

          ggsave(file.path(comp_dir, "pca_plot.png"), p_comp2, width = 8, height = 7, dpi = 300)
          ggsave(file.path(comp_dir, "pca_plot.pdf"), p_comp2, width = 8, height = 7)
          cat("[OK] comparative/pca_plot.png/pdf\n")
        }
      },
      error = function(e) {
        cat("[WARN] Skipped PCA plot:", conditionMessage(e), "\n")
      }
    )
  }
}

# Final summary
cat("\n==========================================\n")
cat("Plots saved to:", plots_dir, "\n")
if (expanded_plots) {
  cat("Subdirectories: basic/, distribution/, qc/, biological/, comparative/\n")
} else {
  cat("Subdirectory: basic/\n")
  cat("Use --expanded-plots for additional analysis plots\n")
}
cat("==========================================\n")
cat("\n--- Summary report complete!\n")
