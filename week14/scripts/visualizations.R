#!/usr/bin/env Rscript
# PCA and Heatmap Visualization
# Author: Generated for RNA-Seq Pipeline
# Description: Creates PCA plot and heatmap from DESeq2 results

# Load required libraries
suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript visualizations.R <deseq_output_dir> <plots_dir>")
}

deseq_dir <- args[1]
plots_dir <- args[2]

# Create output directory
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

cat("=== Generating Visualizations ===\n")

# Load DESeq2 object
cat("Loading DESeq2 object...\n")
dds <- readRDS(file.path(deseq_dir, "dds_object.rds"))

# Variance stabilizing transformation for visualization
cat("Applying variance stabilizing transformation...\n")
vsd <- vst(dds, blind = FALSE)

# ===== PCA Plot =====
cat("\nGenerating PCA plot...\n")

# Get PCA data
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))

# Create PCA plot
pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = condition, label = name)) +
  geom_point(size = 4) +
  geom_text(vjust = -0.8, hjust = 0.5, size = 3) +
  xlab(paste0("PC1: ", percent_var[1], "% variance")) +
  ylab(paste0("PC2: ", percent_var[2], "% variance")) +
  ggtitle("PCA Plot - RNA-Seq Samples") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "bottom"
  ) +
  scale_color_manual(values = c("control" = "#4575b4", "treatment" = "#d73027"))

ggsave(file.path(plots_dir, "pca_plot.png"), pca_plot, width = 8, height = 6, dpi = 300)
ggsave(file.path(plots_dir, "pca_plot.pdf"), pca_plot, width = 8, height = 6)

cat("  Saved: pca_plot.png, pca_plot.pdf\n")

# ===== Sample Distance Heatmap =====
cat("\nGenerating sample distance heatmap...\n")

# Calculate sample distances
sample_dists <- dist(t(assay(vsd)))
sample_dist_matrix <- as.matrix(sample_dists)
rownames(sample_dist_matrix) <- colnames(vsd)
colnames(sample_dist_matrix) <- colnames(vsd)

# Annotation for heatmap
annotation_col <- data.frame(
  Condition = colData(vsd)$condition,
  row.names = colnames(vsd)
)
ann_colors <- list(Condition = c(control = "#4575b4", treatment = "#d73027"))

# Create sample distance heatmap
png(file.path(plots_dir, "sample_distance_heatmap.png"), width = 800, height = 700, res = 120)
pheatmap(
  sample_dist_matrix,
  clustering_distance_rows = sample_dists,
  clustering_distance_cols = sample_dists,
  color = colorRampPalette(rev(brewer.pal(9, "Blues")))(100),
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  main = "Sample Distance Heatmap"
)
dev.off()

pdf(file.path(plots_dir, "sample_distance_heatmap.pdf"), width = 8, height = 7)
pheatmap(
  sample_dist_matrix,
  clustering_distance_rows = sample_dists,
  clustering_distance_cols = sample_dists,
  color = colorRampPalette(rev(brewer.pal(9, "Blues")))(100),
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  main = "Sample Distance Heatmap"
)
dev.off()

cat("  Saved: sample_distance_heatmap.png, sample_distance_heatmap.pdf\n")

# ===== Gene Expression Heatmap (Top DE Genes) =====
cat("\nGenerating gene expression heatmap...\n")

# Load results
res <- read.csv(file.path(deseq_dir, "deseq2_results.csv"), row.names = 1)

# Get top variable genes (or significant genes if available)
sig_genes <- rownames(res[which(res$padj < 0.05), ])

if (length(sig_genes) >= 2) {
  # Use significant genes
  genes_for_heatmap <- head(sig_genes, min(50, length(sig_genes)))
  heatmap_title <- paste0("Top ", length(genes_for_heatmap), " Differentially Expressed Genes")
} else {
  # Fall back to top variable genes
  cat("  Note: Few significant genes found, using top variable genes\n")
  rv <- rowVars(assay(vsd))
  top_genes <- order(rv, decreasing = TRUE)[1:min(50, length(rv))]
  genes_for_heatmap <- rownames(vsd)[top_genes]
  heatmap_title <- "Top 50 Variable Genes"
}

# Extract expression matrix for selected genes
mat <- assay(vsd)[genes_for_heatmap, ]
mat <- mat - rowMeans(mat)  # Center the data

# Create gene expression heatmap
png(file.path(plots_dir, "gene_expression_heatmap.png"), width = 900, height = 1000, res = 120)
pheatmap(
  mat,
  color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = (nrow(mat) <= 50),
  fontsize_row = 8,
  main = heatmap_title
)
dev.off()

pdf(file.path(plots_dir, "gene_expression_heatmap.pdf"), width = 9, height = 10)
pheatmap(
  mat,
  color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = (nrow(mat) <= 50),
  fontsize_row = 8,
  main = heatmap_title
)
dev.off()

cat("  Saved: gene_expression_heatmap.png, gene_expression_heatmap.pdf\n")

# ===== Volcano Plot =====
cat("\nGenerating volcano plot...\n")

# Prepare data for volcano plot
volcano_data <- as.data.frame(res)
volcano_data$gene <- rownames(volcano_data)
volcano_data$significant <- ifelse(
  !is.na(volcano_data$padj) & volcano_data$padj < 0.05,
  ifelse(volcano_data$log2FoldChange > 0, "Up", "Down"),
  "NS"
)

# Create volcano plot
volcano_plot <- ggplot(volcano_data, aes(x = log2FoldChange, y = -log10(pvalue), color = significant)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(
    values = c("Up" = "#d73027", "Down" = "#4575b4", "NS" = "grey60"),
    labels = c("Up" = "Upregulated", "Down" = "Downregulated", "NS" = "Not Significant")
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
  labs(
    title = "Volcano Plot",
    x = "log2 Fold Change",
    y = "-log10(p-value)",
    color = "Status"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "bottom"
  )

ggsave(file.path(plots_dir, "volcano_plot.png"), volcano_plot, width = 8, height = 6, dpi = 300)
ggsave(file.path(plots_dir, "volcano_plot.pdf"), volcano_plot, width = 8, height = 6)

cat("  Saved: volcano_plot.png, volcano_plot.pdf\n")

# ===== MA Plot =====
cat("\nGenerating MA plot...\n")

png(file.path(plots_dir, "ma_plot.png"), width = 800, height = 600, res = 120)
plotMA(results(dds), main = "MA Plot", ylim = c(-5, 5))
dev.off()

pdf(file.path(plots_dir, "ma_plot.pdf"), width = 8, height = 6)
plotMA(results(dds), main = "MA Plot", ylim = c(-5, 5))
dev.off()

cat("  Saved: ma_plot.png, ma_plot.pdf\n")

cat("\n=== Visualization Complete ===\n")
cat("All plots saved to:", plots_dir, "\n")
