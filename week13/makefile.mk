# Makefile for BAM alignment, BigWig generation, and Count Matrix
# Usage: make -f makefile.mk usage
# Author: Anush karampuri
# Date: 2025-12-07
# Description: This Makefile downloads RNA-Seq data, aligns reads to a reference genome,
# generates BigWig files, and creates a count matrix from multiple samples.

# THIS ANNOTATED MAKEFILE USES 'Name' AS GENE IDENTIFIER FOR COUNT MATRIX GENERATION

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

# File paths
REF_GENOME := $(REF_DIR)/$(ACCESSION).fa
GFF_FILE := $(REF_DIR)/$(ACCESSION).gff
R1 := $(READ_DIR)/$(SAMPLE)_1.fastq
R2 := $(READ_DIR)/$(SAMPLE)_2.fastq
BAM_FILE := $(BAM_DIR)/$(SAMPLE).sorted.bam
BIGWIG_FILE := $(BIGWIG_DIR)/$(SAMPLE).bw

# Create directories
$(shell mkdir -p $(READ_DIR) $(REF_DIR) $(BAM_DIR) $(BIGWIG_DIR) $(FASTQC_DIR) $(TEMP_DIR) $(COUNTS_DIR))

# Usage instructions
usage:
	@echo "RNA-Seq Pipeline with Count Matrix Generation"
	@echo ""
	@echo "Step 1: Download genome, annotation, and index (run once)"
	@echo "  make -f makefile.mk genome annotation index"
	@echo ""
	@echo "Step 2: Process all samples in parallel"
	@echo "  cat design.csv | parallel --jobs 3 --colsep , --header : \\"
	@echo "  make -f makefile.mk process_sample SRR={SRR} SAMPLE={name} COVERAGE={coverage}"
	@echo ""
	@echo "Step 3: Generate count matrix"
	@echo "  make -f makefile.mk count_matrix"

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
process_sample: calculate_coverage download_reads fastqc align stats bigwig count_matrix
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



# Clean
clean:
	rm -rf $(READ_DIR) $(REF_DIR) $(BAM_DIR) $(BIGWIG_DIR) $(FASTQC_DIR) $(TEMP_DIR) $(COUNTS_DIR)

.PHONY: usage genome annotation index process_sample calculate_coverage download_reads fastqc align stats bigwig count_matrix clean

# End of Makefile
