# Makefile for BAM alignment and BigWig generation using SRR and sample names
# Usage example:
# make -f makefile.mk SRR=SRR123456 SAMPLE=Sample1 ACC=NC_007793.1 COVERAGE=15

SHELL := bash
.ONESHELL:
.SHELLFLAGS := -eu -o pipefail -c
.DELETE_ON_ERROR:
MAKEFLAGS += --warn-undefined-variables --no-builtin-rules

# --- Required parameters ---
ifndef SRR
  $(error Please provide SRR=<SRR_NUMBER>)
endif

ifndef SAMPLE
  $(error Please provide SAMPLE=<SAMPLE_NAME>)
endif

ifndef ACC
  $(error Please provide ACC=<ACCESSION_NUMBER>)
endif

ifndef COVERAGE
  $(error Please provide COVERAGE=<DESIRED_COVERAGE>)
endif

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

# File paths
REF_GENOME := $(REF_DIR)/$(ACCESSION).fa
R1 := $(READ_DIR)/$(SAMPLE)_1.fastq
R2 := $(READ_DIR)/$(SAMPLE)_2.fastq
BAM_FILE := $(BAM_DIR)/$(SAMPLE).sorted.bam
BIGWIG_FILE := $(BIGWIG_DIR)/$(SAMPLE).bw

# Create directories
$(shell mkdir -p $(READ_DIR) $(REF_DIR) $(BAM_DIR) $(BIGWIG_DIR) $(FASTQC_DIR) $(TEMP_DIR))

# Main pipeline
all: calculate_coverage download_reads fastqc genome index align stats bigwig
	@echo "✅ Pipeline completed for $(SAMPLE) ($(SRR_NUMBER)) with $(DESIRED_COVERAGE)X coverage"
	@echo "$(SAMPLE)" >> $(TEMP_DIR)/processed_samples.txt
	@if [ "$$(wc -l < $(TEMP_DIR)/processed_samples.txt)" -eq "$$(wc -l < design.csv)" ]; then \
		echo "🎯 All samples processed:"; \
		cat $(TEMP_DIR)/processed_samples.txt; \
	fi

# 1. Calculate required reads
calculate_coverage:
	@echo "Calculating reads needed for $(DESIRED_COVERAGE)X coverage..."
	@mkdir -p $(TEMP_DIR)
	@bash -c '\
	GENOME_SIZE=$$(efetch -db nucleotide -id $(ACCESSION) -format fasta | grep -v ">" | tr -d "\\n" | wc -c); \
	echo "Genome size: $$GENOME_SIZE bp"; \
	READ_LENGTH=150; \
	echo "Using read length: $$READ_LENGTH"; \
	REQUIRED_READS=$$(( ($$GENOME_SIZE * $(DESIRED_COVERAGE)) / $$READ_LENGTH )); \
	echo "Required reads: $$REQUIRED_READS"; \
	echo $$REQUIRED_READS > $(TEMP_DIR)/required_reads.txt; \
	'

# 2. Download reads
download_reads:
	@echo "Downloading reads for $(SAMPLE)..."
	@REQUIRED_READS=$$(cat $(TEMP_DIR)/required_reads.txt); \
	mkdir -p $(READ_DIR); \
	fastq-dump -X $$REQUIRED_READS --outdir $(READ_DIR) --split-files $(SRR_NUMBER)
	@mv $(READ_DIR)/$(SRR_NUMBER)_1.fastq $(R1) || true
	@mv $(READ_DIR)/$(SRR_NUMBER)_2.fastq $(R2) || true
	@ls -lh $(READ_DIR)

# 3. FASTQC quality reports
fastqc:
	@echo "Running FastQC for $(SAMPLE)..."
	@if [ -f "$(R1)" ]; then fastqc $(R1) -o $(FASTQC_DIR); fi
	@if [ -f "$(R2)" ]; then fastqc $(R2) -o $(FASTQC_DIR); fi
	@echo "FASTQC reports generated in $(FASTQC_DIR)"

# 4. Download reference genome
genome:
	@echo "Fetching reference genome..."
	efetch -db nucleotide -id $(ACCESSION) -format fasta > $(REF_GENOME)
	@echo "Reference saved: $(REF_GENOME)"

# 5. Index reference genome
index: $(REF_GENOME)
	@echo "Indexing reference..."
	bwa index $(REF_GENOME)
	samtools faidx $(REF_GENOME)

# 6. Align reads
align: $(REF_GENOME)
	@echo "Aligning reads for $(SAMPLE)..."
	@if [ -f "$(R1)" ] && [ -f "$(R2)" ]; then \
		echo "Paired-end alignment"; \
		bwa mem $(REF_GENOME) $(R1) $(R2) | samtools view -bS - | samtools sort -o $(BAM_FILE); \
	elif [ -f "$(R1)" ]; then \
		echo "Single-end alignment"; \
		bwa mem $(REF_GENOME) $(R1) | samtools view -bS - | samtools sort -o $(BAM_FILE); \
	else \
		echo "No FASTQ files found for $(SAMPLE)"; \
		exit 1; \
	fi
	samtools index $(BAM_FILE)
	@echo "Alignment done: $(BAM_FILE)"

# 7. Alignment statistics
stats: $(BAM_FILE)
	@echo "Alignment statistics for $(SAMPLE):"
	samtools flagstat $(BAM_FILE) > $(BAM_DIR)/$(SAMPLE)_alignment_stats.txt
	@echo "Stats saved to $(BAM_DIR)/$(SAMPLE)_alignment_stats.txt"

# 8. Generate BigWig file
bigwig: $(BAM_FILE)
	@echo "Creating BigWig for $(SAMPLE)..."
	bedtools genomecov -ibam $(BAM_FILE) -split -bg | sort -k1,1 -k2,2n > $(BIGWIG_DIR)/$(SAMPLE).bedgraph
	bedGraphToBigWig $(BIGWIG_DIR)/$(SAMPLE).bedgraph $(REF_GENOME).fai $(BIGWIG_FILE)
	@echo "BigWig created: $(BIGWIG_FILE)"

# Clean up
clean:
	rm -rf $(READ_DIR) $(REF_DIR) $(BAM_DIR) $(BIGWIG_DIR) $(FASTQC_DIR) $(TEMP_DIR)
	@echo "🧹 Cleaned all files"

.PHONY: all calculate_coverage download_reads fastqc genome index align stats bigwig clean
