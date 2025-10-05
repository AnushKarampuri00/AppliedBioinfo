#
#
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

## Genbank accession for the reference genome 
ACC=NC_007793.1
# SRR number
SRR=SRR21835896

# Reference genome
REF=refs/${ACC}.fa

# Read 1
R1=reads/${SRR}_1.fastq

# Read 2
R2=reads/${SRR}_2.fastq

# BAM file
BAM=bam/${SRR}.bam

# How many reads to download
N=137931

# Setting useful defaults.
SHELL = bash
.ONESHELL:
.SHELLFLAGS = -eu -o pipefail -c
.DELETE_ON_ERROR:
MAKEFLAGS += --warn-undefined-variables
MAKEFLAGS += --no-builtin-rules

# Prints the usage message
usage:
	@echo "#"
	@echo "# Usage: make [all|refs|fastq|index|align|clean]"    
	@echo "#"
# Obtain the reference genome
refs:
	mkdir -p $(dir ${REF})
    # Download the reference genome
	bio fetch ${ACC} --format fasta > ${REF}

    # Show some information about the genome
	seqkit stats ${REF}

# Download a subset of reads from SRA
# Remove the -X flag to get all data.

fastq:
#	 Make the reads directory
	mkdir -p $(dir ${R1})
    # Download the reads
	fastq-dump -X ${N} --outdir $(dir ${R1}) --split-files ${SRR}
    # Show some information about the reads
	seqkit stats ${R1} ${R2}

# Index the reference genome
index:
	bwa index ${REF}

# Align the reads and convert to BAM. Use 4 threads
# Works for paired-end reads. Modify for single-end reads.
align:
    # Make the BAM directory
	mkdir -p $(dir ${BAM})

    # Align the reads
	bwa mem -t 4 ${REF} ${R1} ${R2} | samtools sort  > ${BAM}

    # Index the BAM file
	samtools index ${BAM}

# Generate alignment statistics
stats:
	samtools flagstat ${BAM}
	samtools coverage ${BAM}

# Clean up generated files
clean:
	rm -rf ${REF} ${R1} ${R2} ${BAM} ${BAM}.bai

# Create necessary directories
all: refs fastq index align stats

# Create necessary directories
.PHONY: all refs fastq index align clean stats


######