SHELL := bash
.ONESHELL:
.SHELLFLAGS := -eu -o pipefail -c
.DELETE_ON_ERROR:
MAKEFLAGS += --warn-undefined-variables --no-builtin-rules

# Parameters
ACCESSION := NC_007793.1
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

# Create directories
$(shell mkdir -p $(READ_DIR) $(REF_DIR) $(BAM_DIR) $(BIGWIG_DIR) $(FASTQC_DIR) $(TEMP_DIR))

# Master target
all: genome index calculate_coverage download_reads fastqc align stats bigwig
	@echo "Finished processing individual samples using "all" target."


# Download reference genome (this will be done only once)
genome:
	@echo "Fetching reference genome..."
	@if [ ! -f "$(REF_GENOME)" ]; then \
		efetch -db nucleotide -id $(ACCESSION) -format fasta > $(REF_GENOME); \
		echo "Reference saved: $(REF_GENOME)"; \
	else \
		echo "Reference genome already exists. Skipping download."; \
	fi

# Index reference genome (only once)
index: $(REF_GENOME)
	@echo "Indexing reference..."
	@if [ ! -f "$(REF_GENOME).bwt" ]; then \
		bwa index $(REF_GENOME); \
		samtools faidx $(REF_GENOME); \
		echo "Indexing done."; \
	else \
		echo "Index files already exist. Skipping indexing."; \
	fi

# Per-sample workflow
process_sample: calculate_coverage download_reads fastqc align stats bigwig
	@echo "Pipeline completed for $(SAMPLE) ($(SRR)) with $(COVERAGE)X coverage"
	@echo "$(SAMPLE)" >> $(TEMP_DIR)/processed_samples.txt
	@if [ "$$(wc -l < $(TEMP_DIR)/processed_samples.txt)" -eq "$$(wc -l < design.csv)" ]; then \
		echo "All samples processed:"; \
		cat $(TEMP_DIR)/processed_samples.txt; \
	fi

# Calculate coverage
calculate_coverage:
	@echo "Calculating reads needed for $(COVERAGE)X coverage..."
	@mkdir -p $(TEMP_DIR)
	@bash -c '\
	GENOME_SIZE=$$(grep -v ">" $(REF_GENOME) | tr -d "\\n" | wc -c); \
	echo "Genome size: $$GENOME_SIZE bp"; \
	READ_LENGTH=150; \
	REQUIRED_READS=$$(( ($$GENOME_SIZE * $(COVERAGE)) / $$READ_LENGTH )); \
	echo "Required reads: $$REQUIRED_READS"; \
	echo $$REQUIRED_READS > $(TEMP_DIR)/$(SAMPLE)_required_reads.txt; \
	'

# Download reads
download_reads:
	@echo "Downloading reads for $(SAMPLE)..."
	@REQUIRED_READS=$$(cat $(TEMP_DIR)/$(SAMPLE)_required_reads.txt); \
	mkdir -p $(READ_DIR); \
	fastq-dump -X $$REQUIRED_READS --outdir $(READ_DIR) --split-files $(SRR)
	@mv $(READ_DIR)/$(SRR)_1.fastq $(READ_DIR)/$(SAMPLE)_1.fastq || true
	@mv $(READ_DIR)/$(SRR)_2.fastq $(READ_DIR)/$(SAMPLE)_2.fastq || true
	@ls -lh $(READ_DIR)

# FASTQC
fastqc:
	@echo "Running FastQC for $(SAMPLE)..."
	@if [ -f "$(READ_DIR)/$(SAMPLE)_1.fastq" ]; then fastqc $(READ_DIR)/$(SAMPLE)_1.fastq -o $(FASTQC_DIR); fi
	@if [ -f "$(READ_DIR)/$(SAMPLE)_2.fastq" ]; then fastqc $(READ_DIR)/$(SAMPLE)_2.fastq -o $(FASTQC_DIR); fi

# Align
align:
	@echo "Aligning reads for $(SAMPLE)..."
	@if [ -f "$(READ_DIR)/$(SAMPLE)_1.fastq" ] && [ -f "$(READ_DIR)/$(SAMPLE)_2.fastq" ]; then \
		bwa mem $(REF_GENOME) $(READ_DIR)/$(SAMPLE)_1.fastq $(READ_DIR)/$(SAMPLE)_2.fastq | samtools view -bS - | samtools sort -o $(BAM_DIR)/$(SAMPLE).sorted.bam; \
	echo "Your FASTQ data was paired end sequenced"; \
	elif [ -f "$(READ_DIR)/$(SAMPLE)_1.fastq" ]; then \
		bwa mem $(REF_GENOME) $(READ_DIR)/$(SAMPLE)_1.fastq | samtools view -bS - | samtools sort -o $(BAM_DIR)/$(SAMPLE).sorted.bam; \
	echo "Your FASTQ data was single end sequenced"; \
	else \
		echo "No FASTQ files found for $(SAMPLE)"; \
		exit 1; \
	fi
	samtools index $(BAM_DIR)/$(SAMPLE).sorted.bam

# Stats
stats:
	samtools flagstat $(BAM_DIR)/$(SAMPLE).sorted.bam > $(BAM_DIR)/$(SAMPLE)_alignment_stats.txt
	@echo "Stats saved to $(BAM_DIR)/$(SAMPLE)_alignment_stats.txt"

# BigWig
bigwig:
	@echo "Creating BigWig for $(SAMPLE)..."
	bedtools genomecov -ibam $(BAM_DIR)/$(SAMPLE).sorted.bam -split -bg | sort -k1,1 -k2,2n > $(BIGWIG_DIR)/$(SAMPLE).bedgraph
	bedGraphToBigWig $(BIGWIG_DIR)/$(SAMPLE).bedgraph $(REF_GENOME).fai $(BIGWIG_DIR)/$(SAMPLE).bw
	@echo "BigWig created: $(BIGWIG_DIR)/$(SAMPLE).bw"

# Clean
clean:
	rm -rf $(READ_DIR) $(REF_DIR) $(BAM_DIR) $(BIGWIG_DIR) $(FASTQC_DIR) $(TEMP_DIR)
	@echo "🧹 Cleaned all files"

.PHONY: all genome index parallel_run sample_pipeline calculate_coverage download_reads fastqc align stats bigwig clean
