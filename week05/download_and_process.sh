#!/bin/bash
set -uexo pipefail

# Step 1: Get the SRR numbers from the publication
# No changes above this line

echo "=== STEP 1: Getting SRR numbers ==="

# This is the code to get the list of SRR numbers and its corresponding experiment data
esearch -db gds -query "GSE78711" | elink -target sra | efetch -format runinfo > runinfo1.csv

echo "SRA run information saved to runinfo.csv"

# see what's in the file
echo "First few lines of runinfo1.csv:"
head -n 3 runinfo1.csv

# Extract just the SRR numbers
echo "SRR numbers found:"
cat runinfo1.csv | cut -f 1 -d ',' | grep "SRR"

# Make a new file with SRR list 
cut -f1 -d',' runinfo1.csv | grep "SRR" > SRR_list.txt

echo "=== Step 1 complete! ==="

# STEP 2: Downloading the subset for 10X coverage
echo "=== STEP 2: Downloading Subset for 10x Coverage ==="

GENOME_SIZE=10389
COVERAGE=10
REQUIRED_BASES=$((GENOME_SIZE * COVERAGE))

echo "Genome size: $GENOME_SIZE bp"
echo "Required for ${COVERAGE}x coverage: $REQUIRED_BASES bases"

SRR_NUMBER=$(cut -f1 -d',' runinfo.csv | grep "SRR" | head -1 || true)

if [ -z "$SRR_NUMBER" ]; then
    echo "Error: No SRR number found!"
    exit 1
fi
echo "Selected SRR number: $SRR_NUMBER"

# Calculate required reads
AVG_READ_LEN=150
REQUIRED_READS=$((REQUIRED_BASES / AVG_READ_LEN))
if [ "$REQUIRED_READS" -le 0 ]; then
    echo "Error: calculated 0 reads, check genome size or avg read length"
    exit 1
fi

echo "Assuming avg read length = $AVG_READ_LEN bp"
echo "Reads required for ${COVERAGE}x coverage: $REQUIRED_READS"

# Run fasterq-dump
if ! command -v fasterq-dump &> /dev/null; then
    echo "Error: fasterq-dump not installed. Install with: conda install -c bioconda sra-tools"
    exit 1
fi

echo "Downloading ~${COVERAGE}x coverage from $SRR_NUMBER..."
fasterq-dump "$SRR_NUMBER" \
    --split-files \
    --progress \
    --row-limit $REQUIRED_READS

echo "=== Download Verification ==="
ls -lh ${SRR_NUMBER}*.fastq || { echo "Download failed"; exit 1; }



