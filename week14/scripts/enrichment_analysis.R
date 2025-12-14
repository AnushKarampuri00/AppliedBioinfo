#!/usr/bin/env Rscript
# Functional Enrichment Analysis
# Author: Generated for RNA-Seq Pipeline
# Description: Performs functional enrichment analysis on differentially expressed genes
# Note: For bacterial genomes, this uses a simplified approach since standard GO/KEGG
# databases are designed for well-annotated model organisms

# Load required libraries
suppressPackageStartupMessages({
  library(ggplot2)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript enrichment_analysis.R <deseq_output_dir> <gff_file> <output_dir>")
}

deseq_dir <- args[1]
gff_file <- args[2]
output_dir <- args[3]

# Create output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("=== Functional Enrichment Analysis ===\n")
cat("DESeq2 results:", deseq_dir, "\n")
cat("GFF annotation:", gff_file, "\n")
cat("Output directory:", output_dir, "\n\n")

# Load significant genes
cat("Loading differentially expressed genes...\n")
sig_genes_file <- file.path(deseq_dir, "significant_genes.csv")

if (!file.exists(sig_genes_file)) {
  cat("Warning: No significant_genes.csv found. Using full results.\n")
  sig_genes <- read.csv(file.path(deseq_dir, "deseq2_results.csv"), row.names = 1)
  sig_genes <- sig_genes[which(sig_genes$padj < 0.05), ]
} else {
  sig_genes <- read.csv(sig_genes_file, row.names = 1)
}

cat("Number of significant genes:", nrow(sig_genes), "\n")

if (nrow(sig_genes) == 0) {
  cat("\nNo significant genes found for enrichment analysis.\n")
  cat("Creating placeholder results...\n")
  
  # Create placeholder file
  write.csv(
    data.frame(Message = "No significant genes found for enrichment analysis"),
    file = file.path(output_dir, "enrichment_results.csv"),
    row.names = FALSE
  )
  
  quit(status = 0)
}

# Separate up and down regulated genes
up_genes <- rownames(sig_genes[sig_genes$log2FoldChange > 0, ])
down_genes <- rownames(sig_genes[sig_genes$log2FoldChange < 0, ])

cat("Upregulated genes:", length(up_genes), "\n")
cat("Downregulated genes:", length(down_genes), "\n")

# ===== Parse GFF for Gene Annotations =====
cat("\nParsing GFF annotation file...\n")

# Read GFF file
gff_lines <- readLines(gff_file)
gff_lines <- gff_lines[!grepl("^#", gff_lines)]  # Remove comments
gff_lines <- gff_lines[gff_lines != ""]  # Remove empty lines

# Parse GFF to extract gene information
parse_gff_attributes <- function(attr_string) {
  attrs <- strsplit(attr_string, ";")[[1]]
  result <- list()
  for (attr in attrs) {
    if (grepl("=", attr)) {
      kv <- strsplit(attr, "=")[[1]]
      if (length(kv) == 2) {
        result[[kv[1]]] <- kv[2]
      }
    }
  }
  return(result)
}

# Extract gene annotations
gene_info <- list()
for (line in gff_lines) {
  fields <- strsplit(line, "\t")[[1]]
  if (length(fields) >= 9 && fields[3] == "gene") {
    attrs <- parse_gff_attributes(fields[9])
    gene_name <- attrs$Name
    if (!is.null(gene_name)) {
      gene_info[[gene_name]] <- list(
        product = ifelse(!is.null(attrs$product), attrs$product, "Unknown"),
        locus_tag = ifelse(!is.null(attrs$locus_tag), attrs$locus_tag, ""),
        gene_biotype = ifelse(!is.null(attrs$gene_biotype), attrs$gene_biotype, "")
      )
    }
  }
}

cat("Genes annotated in GFF:", length(gene_info), "\n")

# ===== Functional Categorization =====
cat("\nPerforming functional categorization...\n")

# Define functional categories based on common bacterial gene naming conventions
# This is a simplified approach for Staphylococcus aureus
functional_categories <- list(
  "DNA Replication/Repair" = c("dna", "gyr", "rec", "pol", "lig", "mut", "uvr", "xer"),
  "Transcription" = c("rpo", "sig", "nusA", "rho"),
  "Translation/Ribosome" = c("rps", "rpl", "rpm", "rrf", "tuf", "fus", "inf", "prf"),
  "tRNA/Aminoacyl-tRNA" = c("^[a-z]{3}S$", "ala", "arg", "asn", "asp", "cys", "gln", "glu", "gly", "his", "ile", "leu", "lys", "met", "phe", "pro", "ser", "thr", "trp", "tyr", "val"),
  "Cell Wall/Membrane" = c("mur", "pbp", "fem", "cap", "tag", "dlt"),
  "Cell Division" = c("fts", "min", "div", "sep"),
  "Metabolism" = c("pyk", "pfk", "eno", "gap", "ldh", "ack", "pta", "pfl"),
  "Amino Acid Metabolism" = c("hut", "aro", "trp", "phe", "tyr", "his", "leu", "ile", "val", "met", "cys", "gly", "ser", "thr", "ala", "glu", "gln", "asp", "asn", "arg", "lys", "pro"),
  "Nucleotide Metabolism" = c("pur", "pyr", "nrd", "adk", "gmk", "ndk"),
  "Energy/Electron Transport" = c("atp", "ndh", "cyd", "qox", "men"),
  "Stress Response" = c("clp", "dna[KJ]", "gro", "htr", "csp", "ctc", "sig[BH]"),
  "Virulence/Pathogenesis" = c("sae", "agr", "sar", "hla", "hlb", "hld", "spa", "coa", "fnb", "clfA", "sdr", "ica"),
  "Transport" = c("opp", "dpp", "abc", "pts", "pst", "mod", "nik", "znu", "mnt"),
  "Regulatory" = c("^[a-z]{3}R$", "ccpA", "codY", "fur", "perR", "sar", "rot"),
  "Unknown/Hypothetical" = c("^SAUSA300")
)

# Categorize genes
categorize_gene <- function(gene_name) {
  gene_lower <- tolower(gene_name)
  for (cat_name in names(functional_categories)) {
    patterns <- functional_categories[[cat_name]]
    for (pattern in patterns) {
      if (grepl(tolower(pattern), gene_lower, ignore.case = TRUE)) {
        return(cat_name)
      }
    }
  }
  return("Other")
}

# Categorize all significant genes
sig_gene_names <- rownames(sig_genes)
gene_categories <- sapply(sig_gene_names, categorize_gene)

# Create category summary
category_counts <- table(gene_categories)
category_df <- data.frame(
  Category = names(category_counts),
  Count = as.numeric(category_counts),
  stringsAsFactors = FALSE
)
category_df <- category_df[order(-category_df$Count), ]

cat("\nFunctional Category Summary:\n")
print(category_df)

# ===== Create Detailed Gene Report =====
cat("\nCreating detailed gene report...\n")

gene_report <- data.frame(
  Gene = sig_gene_names,
  log2FoldChange = sig_genes$log2FoldChange,
  padj = sig_genes$padj,
  Regulation = ifelse(sig_genes$log2FoldChange > 0, "Up", "Down"),
  Category = gene_categories,
  stringsAsFactors = FALSE
)

# Add annotation info if available
gene_report$Product <- sapply(gene_report$Gene, function(g) {
  if (g %in% names(gene_info)) {
    return(gene_info[[g]]$product)
  }
  return("Unknown")
})

gene_report <- gene_report[order(gene_report$padj), ]

# Save gene report
write.csv(gene_report, 
          file = file.path(output_dir, "enrichment_gene_report.csv"),
          row.names = FALSE)

# ===== Category Enrichment Analysis =====
cat("\nPerforming category enrichment analysis...\n")

# Get background gene set (all genes in count matrix)
all_results <- read.csv(file.path(deseq_dir, "deseq2_results.csv"), row.names = 1)
background_genes <- rownames(all_results)
background_categories <- sapply(background_genes, categorize_gene)
background_counts <- table(background_categories)

# Perform hypergeometric test for enrichment
enrichment_results <- data.frame(
  Category = character(),
  DE_Count = numeric(),
  Background_Count = numeric(),
  DE_Total = numeric(),
  Background_Total = numeric(),
  Fold_Enrichment = numeric(),
  P_Value = numeric(),
  stringsAsFactors = FALSE
)

de_total <- length(sig_gene_names)
bg_total <- length(background_genes)

for (cat in unique(c(names(category_counts), names(background_counts)))) {
  de_count <- ifelse(cat %in% names(category_counts), category_counts[cat], 0)
  bg_count <- ifelse(cat %in% names(background_counts), background_counts[cat], 0)
  
  if (bg_count > 0) {
    # Hypergeometric test
    # phyper(q, m, n, k)
    # q = DE genes in category - 1
    # m = background genes in category
    # n = background genes NOT in category
    # k = total DE genes
    p_val <- phyper(de_count - 1, bg_count, bg_total - bg_count, de_total, lower.tail = FALSE)
    
    # Fold enrichment
    expected <- (bg_count / bg_total) * de_total
    fold_enrich <- ifelse(expected > 0, de_count / expected, 0)
    
    enrichment_results <- rbind(enrichment_results, data.frame(
      Category = cat,
      DE_Count = de_count,
      Background_Count = bg_count,
      DE_Total = de_total,
      Background_Total = bg_total,
      Fold_Enrichment = round(fold_enrich, 2),
      P_Value = p_val
    ))
  }
}

# Adjust p-values
enrichment_results$FDR <- p.adjust(enrichment_results$P_Value, method = "BH")
enrichment_results <- enrichment_results[order(enrichment_results$P_Value), ]

# Save enrichment results
write.csv(enrichment_results,
          file = file.path(output_dir, "enrichment_results.csv"),
          row.names = FALSE)

cat("\nEnrichment Results:\n")
print(enrichment_results[enrichment_results$P_Value < 0.1, ])

# ===== Visualization =====
cat("\nCreating enrichment visualizations...\n")

# Bar plot of category counts
category_plot_data <- category_df[category_df$Count > 0, ]
category_plot_data$Category <- factor(category_plot_data$Category, 
                                        levels = category_plot_data$Category)

p1 <- ggplot(category_plot_data, aes(x = reorder(Category, Count), y = Count)) +
  geom_bar(stat = "identity", fill = "#4575b4") +
  coord_flip() +
  labs(
    title = "Functional Categories of DE Genes",
    x = "Functional Category",
    y = "Number of Genes"
  ) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

ggsave(file.path(output_dir, "category_barplot.png"), p1, width = 10, height = 8, dpi = 300)
ggsave(file.path(output_dir, "category_barplot.pdf"), p1, width = 10, height = 8)

# Enrichment dot plot
sig_enrichment <- enrichment_results[enrichment_results$P_Value < 0.2 & enrichment_results$DE_Count > 0, ]

if (nrow(sig_enrichment) > 0) {
  p2 <- ggplot(sig_enrichment, aes(x = Fold_Enrichment, y = reorder(Category, -P_Value))) +
    geom_point(aes(size = DE_Count, color = -log10(P_Value))) +
    scale_color_gradient(low = "blue", high = "red") +
    labs(
      title = "Functional Enrichment Analysis",
      x = "Fold Enrichment",
      y = "Category",
      size = "Gene Count",
      color = "-log10(p-value)"
    ) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
  
  ggsave(file.path(output_dir, "enrichment_dotplot.png"), p2, width = 10, height = 8, dpi = 300)
  ggsave(file.path(output_dir, "enrichment_dotplot.pdf"), p2, width = 10, height = 8)
}

# Up vs Down regulation by category
regulation_by_cat <- aggregate(
  log2FoldChange ~ Category,
  data = gene_report,
  FUN = function(x) c(Up = sum(x > 0), Down = sum(x < 0))
)

# Create stacked bar plot
up_down_data <- data.frame(
  Category = rep(gene_report$Category, 2),
  Regulation = rep(gene_report$Regulation, 2)
)
up_down_summary <- as.data.frame(table(gene_report$Category, gene_report$Regulation))
colnames(up_down_summary) <- c("Category", "Regulation", "Count")

p3 <- ggplot(up_down_summary, aes(x = reorder(Category, -Count), y = Count, fill = Regulation)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  scale_fill_manual(values = c("Down" = "#4575b4", "Up" = "#d73027")) +
  labs(
    title = "Up/Down Regulation by Functional Category",
    x = "Category",
    y = "Number of Genes"
  ) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

ggsave(file.path(output_dir, "regulation_by_category.png"), p3, width = 10, height = 8, dpi = 300)
ggsave(file.path(output_dir, "regulation_by_category.pdf"), p3, width = 10, height = 8)

cat("\n=== Enrichment Analysis Complete ===\n")
cat("Results saved to:", output_dir, "\n")
cat("Files generated:\n")
cat("  - enrichment_gene_report.csv\n")
cat("  - enrichment_results.csv\n")
cat("  - category_barplot.png/pdf\n")
cat("  - enrichment_dotplot.png/pdf\n")
cat("  - regulation_by_category.png/pdf\n")
