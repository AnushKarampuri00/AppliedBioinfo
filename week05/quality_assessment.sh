#!/bin/bash

set -uexo pipefail

# Step 4: Run FASTQC for quality assessment
# No changes above this line

echo "=== STEP 4: FASTQC Quality Assessment ==="

# Find the FASTQ file
FASTQ_FILE=$(ls *.fastq 2>/dev/null | head -1)

if [ -z "$FASTQ_FILE" ]; then
    echo "Error: No FASTQ files found!"
    exit 1
fi

echo "Running FASTQC on: $FASTQ_FILE"

# Create directory for results
mkdir -p fastqc_results



# Run FASTQC
fastqc "$FASTQ_FILE" --outdir=fastqc_results --threads=2

echo "FASTQC analysis complete!"

# List the output files
echo "Output files:"
ls -la fastqc_results/


echo "=== Step 4 complete ==="


#!/bin/bash
set -uexo pipefail

# Step 5: Evaluate FASTQC report and summarize findings

echo "=== STEP 5: FASTQC Report Evaluation ==="

# Find the FASTQC text report
FASTQC_DIR=$(find fastqc_results -name "*_fastqc" -type d | head -1)
SUMMARY_FILE="${FASTQC_DIR}/summary.txt"
DATA_FILE="${FASTQC_DIR}/fastqc_data.txt"

if [ ! -f "$SUMMARY_FILE" ]; then
    echo "Error: FASTQC summary not found. Run step4 first."
    exit 1
fi

echo "FASTQC Summary for: $(basename $FASTQC_DIR)"
echo "=========================================="

# Display the summary
cat "$SUMMARY_FILE"

echo ""
echo "=== KEY METRICS ==="

# Extract important metrics from the data file
if [ -f "$DATA_FILE" ]; then
    echo "Basic Statistics:"
    grep -A 10 ">>Basic Statistics" "$DATA_FILE" | grep -v ">>" | head -10
    
    echo ""
    echo "Per Base Sequence Quality:"
    grep -A 2 ">>Per base sequence quality" "$DATA_FILE" | tail -1
    
    echo ""
    echo "Sequence Length Distribution:"
    grep -A 2 ">>Sequence Length Distribution" "$DATA_FILE" | tail -1
    
    echo ""
    echo "GC Content:"
    grep -A 2 ">>Per sequence GC content" "$DATA_FILE" | tail -1
fi

echo ""
echo "=== INTERPRETATION GUIDE ==="
echo "PASS: Good quality"
echo "WARN: Potential issues - check carefully"
echo "FAIL: Serious problems - data may be unusable"

echo "=== Step 5 complete ==="
# No changes below this line
