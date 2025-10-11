# This makefile automates the process of counting the reads required for a desired coverage
# before downloading, aligning, and generating a BigWig file.
# It also handles both single-end and paired-end reads during alignment.	
# Usage: make -f <filename.mk> SRR=ERR15403387 ACC=NC_007793.1 COVERAGE=15

# --- NO EDIT BELOW THIS LINE ---

SHELL := bash
.ONESHELL:
.SHELLFLAGS := -eu -o pipefail -c
.DELETE_ON_ERROR:
MAKEFLAGS += --warn-undefined-variables --no-builtin-rules

# Check required parameters
ifndef SRR
  $(error Please provide SRR=<SRR_NUMBER>)
endif

ifndef ACC
  $(error Please provide ACC=<ACCESSION_NUMBER>)
endif

ifndef COVERAGE
  $(error Please provide COVERAGE=<DESIRED_COVERAGE> e.g., COVERAGE=10)
endif

# Parameters
SRR_NUMBER := $(SRR)
ACCESSION := $(ACC)
DESIRED_COVERAGE := $(COVERAGE)

# Directories and files
READ_DIR := reads
REF_DIR := refs
BAM_DIR := bam
BIGWIG_DIR := bigwig
TEMP_DIR := temp

REF_GENOME := $(REF_DIR)/$(ACCESSION).fa
R1 := $(READ_DIR)/$(SRR_NUMBER)_1.fastq
R2 := $(READ_DIR)/$(SRR_NUMBER)_2.fastq
BAM_FILE := $(BAM_DIR)/$(SRR_NUMBER).sorted.bam
BIGWIG_FILE := $(BIGWIG_DIR)/$(SRR_NUMBER).bw

# Create directories
$(shell mkdir -p $(READ_DIR) $(REF_DIR) $(BAM_DIR) $(BIGWIG_DIR) $(TEMP_DIR))

# Main pipeline
all: calculate_coverage download_ref genome index align stats bigwig
	@echo "✅ Pipeline completed for $(SRR_NUMBER) with $(DESIRED_COVERAGE)X coverage"

# 1. Calculate required reads based on coverage
calculate_coverage:
	@echo "Calculating reads needed for $(DESIRED_COVERAGE)X coverage..."
	@mkdir -p $(TEMP_DIR)
	
	@bash -c '\
	echo "📏 Fetching genome size for $(ACCESSION)..."; \
	GENOME_SIZE=$$(efetch -db nucleotide -id $(ACCESSION) -format fasta | grep -v ">" | tr -d "\\n" | wc -c); \
	echo "Genome size: $$GENOME_SIZE bp"; \
	\
	echo "Fetching read length from SRA..."; \
	fastq-dump -X 1 --outdir $(TEMP_DIR) --split-files $(SRR_NUMBER) 2>/dev/null || true; \
	if [ -f "$(TEMP_DIR)/$(SRR_NUMBER).fastq" ]; then \
		READ_LENGTH=$$(head -n 4 $(TEMP_DIR)/$(SRR_NUMBER).fastq | tail -n 1 | wc -c); \
		READ_LENGTH=$$((READ_LENGTH - 1)); \
		echo "Read length: $$READ_LENGTH bp"; \
	elif [ -f "$(TEMP_DIR)/$(SRR_NUMBER)_1.fastq" ]; then \
		READ_LENGTH=$$(head -n 4 $(TEMP_DIR)/$(SRR_NUMBER)_1.fastq | tail -n 1 | wc -c); \
		READ_LENGTH=$$((READ_LENGTH - 1)); \
		echo "Read length: $$READ_LENGTH bp"; \
	else \
		echo "Could not determine read length, using default 150 bp"; \
		READ_LENGTH=150; \
	fi; \
	\
	echo "🧪 Calculating required reads..."; \
	REQUIRED_READS=$$(( ($$GENOME_SIZE * $(DESIRED_COVERAGE)) / $$READ_LENGTH )); \
	echo "Required reads for $(DESIRED_COVERAGE)X coverage: $$REQUIRED_READS"; \
	\
	echo "$$REQUIRED_READS" > $(TEMP_DIR)/required_reads.txt; \
	echo "Calculation complete: $$REQUIRED_READS reads needed"; \
	\
	rm -f $(TEMP_DIR)/*.fastq; \
	'

# 2. Download reads based on calculated requirement
download_ref:
	@echo "Downloading reads..."
	@REQUIRED_READS=$$(cat $(TEMP_DIR)/required_reads.txt); \
	echo "Downloading $$REQUIRED_READS reads for $(DESIRED_COVERAGE)X coverage..."; \
	mkdir -p $(dir ${R1}); \
	fastq-dump -X $$REQUIRED_READS --outdir $(dir ${R1}) --split-files ${SRR_NUMBER}
	@echo "Checking downloaded files..."
	@ls -la $(READ_DIR)/*

# 3. Download reference genome
genome:
	@echo "Fetching reference genome..."
	efetch -db nucleotide -id $(ACCESSION) -format fasta > $(REF_GENOME)
	@echo "Reference saved: $(REF_GENOME)"

# 4. Index reference genome
index: $(REF_GENOME)
	@echo "Indexing reference..."
	bwa index $(REF_GENOME)
	samtools faidx $(REF_GENOME)

# 5. Align reads
align: $(REF_GENOME)
	@echo "Aligning reads..."
	@if [ -f "$(R1)" ] && [ -f "$(R2)" ]; then \
		echo "Aligning paired-end reads: $(R1) and $(R2)"; \
		bwa mem $(REF_GENOME) $(R1) $(R2) | \
		samtools view -bS - | \
		samtools sort -o $(BAM_FILE); \
	elif [ -f "$(R1)" ]; then \
		echo "Aligning single-end reads: $(R1)"; \
		bwa mem $(REF_GENOME) $(R1) | \
		samtools view -bS - | \
		samtools sort -o $(BAM_FILE); \
	else \
		echo "Available files:"; \
		ls -la $(READ_DIR)/*; \
		echo "No suitable read files found!"; \
		exit 1; \
	fi
	samtools index $(BAM_FILE)
	@echo "Alignment done: $(BAM_FILE)"

# 6. Alignment statistics
stats: $(BAM_FILE)
	@echo "Alignment statistics:"
	samtools flagstat $(BAM_FILE)

# 7. Generate BigWig file
bigwig: $(BAM_FILE)
	@echo "Creating BigWig..."
	bedtools genomecov -ibam $(BAM_FILE) -split -bg | \
	sort -k1,1 -k2,2n > $(BIGWIG_DIR)/$(SRR_NUMBER).bedgraph
	bedGraphToBigWig $(BIGWIG_DIR)/$(SRR_NUMBER).bedgraph $(REF_GENOME).fai $(BIGWIG_FILE)
	@echo "BigWig created: $(BIGWIG_FILE)"

# Show coverage information
coverage_info:
	@echo "Coverage Information:"
	@if [ -f "$(TEMP_DIR)/required_reads.txt" ]; then \
		REQUIRED_READS=$$(cat $(TEMP_DIR)/required_reads.txt); \
		echo "Calculated reads needed: $$REQUIRED_READS"; \
	else \
		echo "Run 'make calculate_coverage' first"; \
	fi

# Clean up
clean:
	rm -rf $(READ_DIR) $(REF_DIR) $(BAM_DIR) $(BIGWIG_DIR) $(TEMP_DIR)
	@echo "🧹 Cleaned all files"

.PHONY: all calculate_coverage download genome index align stats bigwig coverage_info clean

# End of Makefile!
