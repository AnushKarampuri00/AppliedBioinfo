# Makefile for BAM alignment, BigWig generation, Count Matrix, and Differential Expression Analysis
# Usage: make -f Makefile.mk usage
# Author: Anush karampuri
# Date: 2025-12-07 (Updated: 2025-12-14)
# Description: This Makefile downloads RNA-Seq data, aligns reads to a reference genome,
# generates BigWig files, creates a count matrix, performs differential expression analysis,
# and produces visualizations and enrichment analysis.

# THIS ANNOTATED MAKEFILE USES 'Name' AS GENE IDENTIFIER FOR COUNT MATRIX GENERATION

# NO CHANGE REQUIRED BELOW THIS LINE UNLESS CUSTOMIZING
 
SHELL := bash
.ONESHELL:
.SHELLFLAGS := -eu -o pipefail -c
.DELETE_ON_ERROR:
MAKEFLAGS += --warn-undefined-variables --no-builtin-rules

# Default accession
ACC ?= NC_007793.1

# Parameters
SRR_NUMBER := $(SRR)
SAMPLE := $(SAMPLE)
ACCESSION := $(ACC)
DESIRED_COVERAGE := $(COVERAGE)

# Directories
READ_DIR := reads
REF_DIR := refs
BAM_DIR := bam
BIGWIG_DIR := bigwig
FASTQC_DIR := fastqc
TEMP_DIR := temp
COUNTS_DIR := counts
DESEQ_DIR := deseq2_results
PLOTS_DIR := plots
ENRICHMENT_DIR := enrichment

# File paths
REF_GENOME := $(REF_DIR)/$(ACCESSION).fa
GFF_FILE := $(REF_DIR)/$(ACCESSION).gff
R1 := $(READ_DIR)/$(SAMPLE)_1.fastq
R2 := $(READ_DIR)/$(SAMPLE)_2.fastq
BAM_FILE := $(BAM_DIR)/$(SAMPLE).sorted.bam
BIGWIG_FILE := $(BIGWIG_DIR)/$(SAMPLE).bw

# R Scripts directory
SCRIPTS_DIR := scripts

# Create directories
$(shell mkdir -p $(READ_DIR) $(REF_DIR) $(BAM_DIR) $(BIGWIG_DIR) $(FASTQC_DIR) $(TEMP_DIR) $(COUNTS_DIR) $(DESEQ_DIR) $(PLOTS_DIR) $(ENRICHMENT_DIR) $(SCRIPTS_DIR))

# Usage instructions
usage:
	@echo "============================================================"
	@echo "       RNA-Seq Pipeline with Differential Expression        "
	@echo "============================================================"
	@echo ""
	@echo "STEP 1: Download genome, annotation, and index (run once)"
	@echo "  make -f Makefile.mk genome annotation index"
	@echo ""
	@echo "STEP 2: Process all samples in parallel"
	@echo "  cat design.csv | parallel --jobs 3 --colsep , --header : \\"
	@echo "  make -f Makefile.mk process_sample SRR={SRR} SAMPLE={name} COVERAGE={coverage}"
	@echo ""
	@echo "STEP 3: Generate count matrix (after all samples processed)"
	@echo "  make -f Makefile.mk count_matrix"
	@echo ""
	@echo "STEP 4: Run differential expression analysis (requires stats env)"
	@echo "  conda activate stats"
	@echo "  make -f Makefile.mk differential_expression"
	@echo ""
	@echo "STEP 5: Generate visualizations"
	@echo "  make -f Makefile.mk visualizations"
	@echo ""
	@echo "STEP 6: Run functional enrichment analysis"
	@echo "  make -f Makefile.mk enrichment"
	@echo ""
	@echo "Or run all DE analysis steps at once:"
	@echo "  make -f Makefile.mk full_de_analysis"
	@echo ""
	@echo "============================================================"

# Download reference genome
genome:
	@echo "Fetching reference genome..."
	efetch -db nucleotide -id $(ACCESSION) -format fasta > $(REF_GENOME)

# Download GFF annotation
annotation:
	@echo "Fetching GFF annotation..."
	wget -q -O $(GFF_FILE) "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=gff3&id=$(ACCESSION)"

# Index reference genome
index: $(REF_GENOME)
	@echo "Indexing reference..."
	bwa index $(REF_GENOME)
	samtools faidx $(REF_GENOME)

# Process single sample (called by parallel)
process_sample: calculate_coverage download_reads fastqc align stats bigwig
	@echo "✅ Done: $(SAMPLE)"

# Calculate required reads
calculate_coverage:
	@bash -c '\
	GENOME_SIZE=$$(efetch -db nucleotide -id $(ACCESSION) -format fasta | grep -v ">" | tr -d "\\n" | wc -c); \
	REQUIRED_READS=$$(( ($$GENOME_SIZE * $(DESIRED_COVERAGE)) / 150 )); \
	echo $$REQUIRED_READS > $(TEMP_DIR)/$(SAMPLE)_reads.txt'

# Download reads
download_reads:
	@REQUIRED_READS=$$(cat $(TEMP_DIR)/$(SAMPLE)_reads.txt); \
	fastq-dump -X $$REQUIRED_READS --outdir $(READ_DIR) --split-files $(SRR_NUMBER)
	@mv $(READ_DIR)/$(SRR_NUMBER)_1.fastq $(R1) || true
	@mv $(READ_DIR)/$(SRR_NUMBER)_2.fastq $(R2) || true

# FASTQC
fastqc:
	@if [ -f "$(R1)" ]; then fastqc $(R1) -o $(FASTQC_DIR); fi
	@if [ -f "$(R2)" ]; then fastqc $(R2) -o $(FASTQC_DIR); fi

# Align reads
align: $(REF_GENOME)
	@if [ -f "$(R1)" ] && [ -f "$(R2)" ]; then \
		bwa mem $(REF_GENOME) $(R1) $(R2) | samtools view -bS - | samtools sort -o $(BAM_FILE); \
	elif [ -f "$(R1)" ]; then \
		bwa mem $(REF_GENOME) $(R1) | samtools view -bS - | samtools sort -o $(BAM_FILE); \
	fi
	samtools index $(BAM_FILE)

# Alignment stats
stats: $(BAM_FILE)
	samtools flagstat $(BAM_FILE) > $(BAM_DIR)/$(SAMPLE)_stats.txt

# BigWig
bigwig: $(BAM_FILE)
	bedtools genomecov -ibam $(BAM_FILE) -split -bg | sort -k1,1 -k2,2n > $(BIGWIG_DIR)/$(SAMPLE).bedgraph
	bedGraphToBigWig $(BIGWIG_DIR)/$(SAMPLE).bedgraph $(REF_GENOME).fai $(BIGWIG_FILE)

# Count matrix (run after all samples processed)
count_matrix:
	@echo "Generating count matrix..."
	featureCounts -p -a $(GFF_FILE) -o $(COUNTS_DIR)/counts_raw.txt -t gene -g Name -F GFF $(BAM_DIR)/*.sorted.bam
	@cut -f1,7- $(COUNTS_DIR)/counts_raw.txt | tail -n +2 | sed 's|$(BAM_DIR)/||g' | sed 's|.sorted.bam||g' > $(COUNTS_DIR)/count_matrix.txt
	@echo "✅ Count matrix: $(COUNTS_DIR)/count_matrix.txt"

# ============================================================
#          DIFFERENTIAL EXPRESSION ANALYSIS TARGETS
# ============================================================
# Note: These targets require the 'stats' conda environment
# Run: conda activate stats before executing these targets

# Run DESeq2 differential expression analysis
differential_expression: $(COUNTS_DIR)/count_matrix.txt design.csv
	@echo "============================================================"
	@echo "Running DESeq2 Differential Expression Analysis..."
	@echo "============================================================"
	@echo ""
	@if [ ! -f "$(SCRIPTS_DIR)/deseq2_analysis.R" ]; then \
		echo "Error: R script not found at $(SCRIPTS_DIR)/deseq2_analysis.R"; \
		echo "Please ensure the R scripts are in the scripts/ directory."; \
		exit 1; \
	fi
	Rscript $(SCRIPTS_DIR)/deseq2_analysis.R $(COUNTS_DIR)/count_matrix.txt design.csv $(DESEQ_DIR)
	@echo ""
	@echo "✅ DESeq2 analysis complete!"
	@echo "   Results saved to: $(DESEQ_DIR)/"

# Generate visualizations (PCA, heatmaps, volcano plot)
visualizations: $(DESEQ_DIR)/dds_object.rds
	@echo "============================================================"
	@echo "Generating Visualizations..."
	@echo "============================================================"
	@echo ""
	@if [ ! -f "$(SCRIPTS_DIR)/visualizations.R" ]; then \
		echo "Error: R script not found at $(SCRIPTS_DIR)/visualizations.R"; \
		exit 1; \
	fi
	Rscript $(SCRIPTS_DIR)/visualizations.R $(DESEQ_DIR) $(PLOTS_DIR)
	@echo ""
	@echo "✅ Visualizations complete!"
	@echo "   Plots saved to: $(PLOTS_DIR)/"

# Run functional enrichment analysis
enrichment: $(DESEQ_DIR)/significant_genes.csv $(GFF_FILE)
	@echo "============================================================"
	@echo "Running Functional Enrichment Analysis..."
	@echo "============================================================"
	@echo ""
	@if [ ! -f "$(SCRIPTS_DIR)/enrichment_analysis.R" ]; then \
		echo "Error: R script not found at $(SCRIPTS_DIR)/enrichment_analysis.R"; \
		exit 1; \
	fi
	Rscript $(SCRIPTS_DIR)/enrichment_analysis.R $(DESEQ_DIR) $(GFF_FILE) $(ENRICHMENT_DIR)
	@echo ""
	@echo "✅ Enrichment analysis complete!"
	@echo "   Results saved to: $(ENRICHMENT_DIR)/"

# Full differential expression analysis pipeline
full_de_analysis: differential_expression visualizations enrichment
	@echo ""
	@echo "============================================================"
	@echo "       FULL DIFFERENTIAL EXPRESSION ANALYSIS COMPLETE       "
	@echo "============================================================"
	@echo ""
	@echo "Results Summary:"
	@echo "  • DESeq2 results:    $(DESEQ_DIR)/"
	@echo "  • Visualizations:    $(PLOTS_DIR)/"
	@echo "  • Enrichment:        $(ENRICHMENT_DIR)/"
	@echo ""
	@echo "Key output files:"
	@echo "  • $(DESEQ_DIR)/deseq2_results.csv"
	@echo "  • $(DESEQ_DIR)/significant_genes.csv"
	@echo "  • $(PLOTS_DIR)/pca_plot.png"
	@echo "  • $(PLOTS_DIR)/gene_expression_heatmap.png"
	@echo "  • $(PLOTS_DIR)/volcano_plot.png"
	@echo "  • $(ENRICHMENT_DIR)/enrichment_results.csv"
	@echo ""

# Generate summary report
report:
	@echo "============================================================"
	@echo "               ANALYSIS SUMMARY REPORT                      "
	@echo "============================================================"
	@echo ""
	@if [ -f "$(DESEQ_DIR)/deseq2_results.csv" ]; then \
		echo "DESeq2 Results:"; \
		echo "  Total genes analyzed: $$(tail -n +2 $(DESEQ_DIR)/deseq2_results.csv | wc -l)"; \
		echo "  Significant genes (padj < 0.05): $$(tail -n +2 $(DESEQ_DIR)/significant_genes.csv 2>/dev/null | wc -l || echo 0)"; \
	else \
		echo "DESeq2 results not found. Run 'make differential_expression' first."; \
	fi
	@echo ""
	@if [ -d "$(PLOTS_DIR)" ]; then \
		echo "Visualizations generated:"; \
		ls -1 $(PLOTS_DIR)/*.png 2>/dev/null || echo "  No plots found"; \
	fi
	@echo ""
	@if [ -d "$(ENRICHMENT_DIR)" ]; then \
		echo "Enrichment files:"; \
		ls -1 $(ENRICHMENT_DIR)/*.csv 2>/dev/null || echo "  No enrichment results found"; \
	fi
	@echo ""

# Clean all generated files
clean:
	rm -rf $(READ_DIR) $(REF_DIR) $(BAM_DIR) $(BIGWIG_DIR) $(FASTQC_DIR) $(TEMP_DIR) $(COUNTS_DIR)

# Clean only DE analysis files (keep alignment/counts)
clean_de:
	rm -rf $(DESEQ_DIR) $(PLOTS_DIR) $(ENRICHMENT_DIR)

# Clean everything including DE analysis
clean_all: clean clean_de
	@echo "All generated files removed."

.PHONY: usage genome annotation index process_sample calculate_coverage download_reads \
        fastqc align stats bigwig count_matrix differential_expression visualizations \
        enrichment full_de_analysis report clean clean_de clean_all

# End of Makefile