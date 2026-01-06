#!/usr/bin/env Rscript
# Generate diverse BED files for NanoporeToBED demo
# Creates files with varying methylation levels, coverage, and CpG counts

set.seed(42)

output_dir <- "demo_data"
dir.create(output_dir, showWarnings = FALSE)

# Sample configuration with diversity
samples <- data.frame(
    name = c(
        "healthy_01", "healthy_02", "healthy_03", "healthy_04",
        "treated_01", "treated_02", "treated_03", "treated_04",
        "control_01", "control_02", "control_03", "control_04"
    ),
    mean_meth = c(
        0.25, 0.28, 0.22, 0.30, # healthy: low methylation
        0.65, 0.72, 0.58, 0.68, # treated: high methylation
        0.42, 0.48, 0.38, 0.45
    ), # control: medium methylation
    sd_meth = c(
        0.15, 0.18, 0.12, 0.16,
        0.20, 0.18, 0.22, 0.19,
        0.18, 0.20, 0.15, 0.17
    ),
    mean_cov = c(
        35, 28, 42, 30,
        22, 18, 25, 20,
        45, 38, 50, 42
    ),
    cpg_count = c(
        18000, 22000, 15000, 20000,
        25000, 20000, 28000, 23000,
        30000, 35000, 28000, 32000
    ),
    stringsAsFactors = FALSE
)

# Chromosomes
chromosomes <- paste0("chr", 1:22)

cat("Generating diverse BED files for NanoporeToBED demo...\n\n")

for (i in 1:nrow(samples)) {
    s <- samples[i, ]
    cat(sprintf(
        "Creating %s: ~%.0f%% meth, ~%.0fx cov, %d CpGs\n",
        s$name, s$mean_meth * 100, s$mean_cov, s$cpg_count
    ))

    # Generate positions
    n_cpg <- s$cpg_count
    chr <- sample(chromosomes, n_cpg,
        replace = TRUE,
        prob = c(rep(0.06, 10), rep(0.04, 12))
    ) # Weighted by chr size

    # Generate positions within chromosomes
    start <- sort(sample(1:250000000, n_cpg, replace = TRUE))
    end <- start + 1

    # Generate methylation values (beta distribution-like via truncated normal)
    meth_raw <- rnorm(n_cpg, mean = s$mean_meth, sd = s$sd_meth)
    meth <- pmax(0, pmin(1, meth_raw)) * 100 # Clamp to 0-100%

    # Generate coverage (negative binomial-like)
    cov <- pmax(1, rnbinom(n_cpg, size = 5, mu = s$mean_cov))

    # Calculate Cs and Ts from coverage and methylation
    n_mod <- round(cov * meth / 100)
    n_canonical <- cov - n_mod

    # Strand
    strand <- sample(c("+", "-"), n_cpg, replace = TRUE)

    # Create 12-column bedMethyl format
    bed_data <- data.frame(
        chr = chr,
        start = start,
        end = end,
        name = "m",
        score = round(meth),
        strand = strand,
        thickStart = start,
        thickEnd = end,
        itemRgb = "255,0,0",
        coverage = cov,
        percent_mod = round(meth, 1),
        n_mod = n_mod,
        n_canonical = n_canonical,
        n_other = 0,
        n_delete = 0,
        n_fail = 0,
        n_diff = 0,
        n_nocall = 0,
        stringsAsFactors = FALSE
    )

    # Sort by chr and position
    bed_data <- bed_data[order(bed_data$chr, bed_data$start), ]

    # Write BED file
    output_file <- file.path(output_dir, paste0(s$name, ".CpG.bed"))
    write.table(bed_data, output_file,
        sep = "\t", quote = FALSE,
        row.names = FALSE, col.names = FALSE
    )
}

cat("\n✓ Created", nrow(samples), "diverse BED files in", output_dir, "\n")
cat("\nSummary:\n")
print(samples[, c("name", "mean_meth", "mean_cov", "cpg_count")])
