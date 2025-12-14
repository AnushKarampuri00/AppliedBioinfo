#!/usr/bin/env Rscript
# DESeq2 Differential Expression Analysis
# Author: Generated for RNA-Seq Pipeline
# Description: Performs differential expression analysis using DESeq2

# Load required libraries
suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript deseq2_analysis.R <count_matrix.txt> <design.csv> <output_dir>")
}

count_file <- args[1]
design_file <- args[2]
output_dir <- args[3]

# Create output directory if it doesn't exist
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("=== DESeq2 Differential Expression Analysis ===\n")
cat("Count matrix:", count_file, "\n")
cat("Design file:", design_file, "\n")
cat("Output directory:", output_dir, "\n\n")

# Read count matrix
cat("Reading count matrix...\n")
counts <- read.table(count_file, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
cat("  Dimensions:", nrow(counts), "genes x", ncol(counts), "samples\n")

# Read design file
cat("Reading design file...\n")
design <- read.csv(design_file, stringsAsFactors = FALSE)
coldata <- data.frame(
  row.names = design$name,
  condition = factor(design$condition)
)
cat("  Conditions:", paste(levels(coldata$condition), collapse = ", "), "\n")

# Ensure column order matches
counts <- counts[, rownames(coldata)]

# Filter low count genes (at least 10 counts in at least 2 samples)
cat("\nFiltering low-count genes...\n")
keep <- rowSums(counts >= 10) >= 2
counts_filtered <- counts[keep, ]
cat("  Genes before filtering:", nrow(counts), "\n")
cat("  Genes after filtering:", nrow(counts_filtered), "\n")

# Create DESeq2 dataset
cat("\nCreating DESeq2 dataset...\n")
dds <- DESeqDataSetFromMatrix(
  countData = counts_filtered,
  colData = coldata,
  design = ~ condition
)

# Set reference level (control)
dds$condition <- relevel(dds$condition, ref = "control")

# Run DESeq2
cat("Running DESeq2 analysis...\n")
dds <- DESeq(dds)

# Get results
cat("Extracting results...\n")
res <- results(dds, contrast = c("condition", "treatment", "control"))
res_ordered <- res[order(res$padj), ]

# Summary statistics
cat("\n=== Results Summary ===\n")
summary(res)

# Count significant genes
sig_genes <- sum(res$padj < 0.05, na.rm = TRUE)
up_genes <- sum(res$padj < 0.05 & res$log2FoldChange > 0, na.rm = TRUE)
down_genes <- sum(res$padj < 0.05 & res$log2FoldChange < 0, na.rm = TRUE)

cat("\nSignificant genes (padj < 0.05):", sig_genes, "\n")
cat("  Upregulated:", up_genes, "\n")
cat("  Downregulated:", down_genes, "\n")

# Save results
cat("\nSaving results...\n")

# Full results
write.csv(as.data.frame(res_ordered), 
          file = file.path(output_dir, "deseq2_results.csv"),
          row.names = TRUE)

# Significant genes only
sig_results <- as.data.frame(res_ordered[which(res_ordered$padj < 0.05), ])
write.csv(sig_results,
          file = file.path(output_dir, "significant_genes.csv"),
          row.names = TRUE)

# Save normalized counts
normalized_counts <- counts(dds, normalized = TRUE)
write.csv(normalized_counts,
          file = file.path(output_dir, "normalized_counts.csv"),
          row.names = TRUE)

# Save DESeq2 object for downstream analysis
saveRDS(dds, file = file.path(output_dir, "dds_object.rds"))

cat("  deseq2_results.csv - Full results\n")
cat("  significant_genes.csv - Significant genes (padj < 0.05)\n")
cat("  normalized_counts.csv - Normalized count matrix\n")
cat("  dds_object.rds - DESeq2 object\n")

cat("\n=== Analysis Complete ===\n")
