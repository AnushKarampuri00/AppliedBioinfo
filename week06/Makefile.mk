# Set the shell the commands run in.
SHELL = bash

# Execute all commands in a single shell.
.ONESHELL:

# Run the shell with strict error checking.
.SHELLFLAGS = -eu -o pipefail -c

# Delete target files if the command fails.
.DELETE_ON_ERROR:

# Warn if a variable is not defined.
MAKEFLAGS += --warn-undefined-variables

# Disable built-in rules.
MAKEFLAGS += --no-builtin-rules

# Parameters
SRR_NUMBER = SRR21835901
REQUIRED_READS = 137931
ACC = NC_007793.1
READ_DIR = reads
REF_GENOME = refs/$(ACC).fa
BAM_FILE = $(SRR_NUMBER).sorted.bam

#----NO CHANGES BELOW THIS LINE----

# Create directories if they don't exist
$(shell mkdir -p $(READ_DIR) refs)

# Main pipeline target
all: download genome index align check
	@echo "Pipeline completed successfully!"

# 1. Download sequencing reads
download:
	@echo "Downloading reads for $(SRR_NUMBER)..."
	fastq-dump -X $(REQUIRED_READS) --split-files --outdir $(READ_DIR) $(SRR_NUMBER)
	@echo "Finished downloading $(SRR_NUMBER)"

# 2. Get reference genome
genome:
	@echo "Fetching reference genome $(ACC)..."
	efetch -db nucleotide -id $(ACC) -format fasta > $(REF_GENOME)
	@echo "Reference genome downloaded: $(ACC)"

# 3. Index the genome
index: $(REF_GENOME)
	@echo "Indexing genome..."
	bwa index $(REF_GENOME)
	@echo "Genome indexed successfully"

# 4. Align reads and generate sorted BAM file
align: $(REF_GENOME) $(READ_DIR)/$(SRR_NUMBER)_1.fastq $(READ_DIR)/$(SRR_NUMBER)_2.fastq
	@echo "Starting alignment for $(SRR_NUMBER)..."
	bwa mem $(REF_GENOME) $(READ_DIR)/$(SRR_NUMBER)_1.fastq $(READ_DIR)/$(SRR_NUMBER)_2.fastq | \
		samtools view -bS - | \
		samtools sort -o $(BAM_FILE)
	samtools index $(BAM_FILE)
	@echo "Alignment complete: $(BAM_FILE) and index created"

# 5. Check alignment statistics
check: $(BAM_FILE)
	@echo "Alignment statistics:"
	samtools flagstat $(BAM_FILE)


.PHONY: all download genome index align check clean 
