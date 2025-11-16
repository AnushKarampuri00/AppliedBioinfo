#!/bin/bash

# Set the chromosome and region of interest
CHR=chr18
REGION=chr18:51,030,213-51,085,042

# Create directory structure
mkdir -p refs bam reports

# ===== REFERENCE GENOME =====
echo "Downloading reference genome..."
REF_URL=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/references/GRCh38/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta.gz
REF_FULL=refs/GRCh38.fasta.gz
REF_SUBSET=refs/${CHR}.fasta

# Download full reference
curl -L ${REF_URL} -o ${REF_FULL}

# Extract only the chromosome of interest
echo "Extracting ${CHR} from reference..."
samtools faidx ${REF_FULL} ${CHR} > ${REF_SUBSET}

# Index the subset reference
samtools faidx ${REF_SUBSET}

echo "Reference genome prepared: ${REF_SUBSET}"

# ===== BAM FILES DOWNLOAD =====

echo "=========================================="
echo "Processing Illumina Control (Normal)..."
BAM_URL="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_somatic/HG008/Liss_lab/BCM_Illumina-WGS_20240313/HG008-N-D_Illumina_169x_GRCh38-GIABv3.bam"
BAM_OUTPUT="bam/HG008-N-D_Illumina_${CHR}.bam"
echo "Downloading region ${REGION}..."
samtools view -b "${BAM_URL}" ${REGION} > "${BAM_OUTPUT}"
echo "Indexing..."
samtools index "${BAM_OUTPUT}"
echo "Statistics:"
samtools flagstat "${BAM_OUTPUT}"
echo ""

echo "=========================================="
echo "Processing Illumina Cancer (Tumor)..."
BAM_URL="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_somatic/HG008/Liss_lab/BCM_Illumina-WGS_20240313/HG008-T_Illumina_195x_GRCh38-GIABv3.bam"
BAM_OUTPUT="bam/HG008-T_Illumina_${CHR}.bam"
echo "Downloading region ${REGION}..."
samtools view -b "${BAM_URL}" ${REGION} > "${BAM_OUTPUT}"
echo "Indexing..."
samtools index "${BAM_OUTPUT}"
echo "Statistics:"
samtools flagstat "${BAM_OUTPUT}"
echo ""

echo "=========================================="
echo "Processing Element Control (Normal)..."
BAM_URL="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_somatic/HG008/Liss_lab/Element-AVITI-20241216/HG008-N-D_Element-StdInsert_77x_GRCh38-GIABv3.bam"
BAM_OUTPUT="bam/HG008-N-D_Element_${CHR}.bam"
echo "Downloading region ${REGION}..."
samtools view -b "${BAM_URL}" ${REGION} > "${BAM_OUTPUT}"
echo "Indexing..."
samtools index "${BAM_OUTPUT}"
echo "Statistics:"
samtools flagstat "${BAM_OUTPUT}"
echo ""

echo "=========================================="
echo "Processing Element Cancer (Tumor)..."
BAM_URL="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_somatic/HG008/Liss_lab/Element-AVITI-20241216/HG008-T_Element-StdInsert_111x_GRCh38-GIABv3.bam"
BAM_OUTPUT="bam/HG008-T_Element_${CHR}.bam"
echo "Downloading region ${REGION}..."
samtools view -b "${BAM_URL}" ${REGION} > "${BAM_OUTPUT}"
echo "Indexing..."
samtools index "${BAM_OUTPUT}"
echo "Statistics:"
samtools flagstat "${BAM_OUTPUT}"
echo ""

echo "=========================================="
echo "Processing PacBio Control (Normal)..."
BAM_URL="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_somatic/HG008/Liss_lab/PacBio_Onso_20240415/HG008-N-D_Pacbio-onso_48x_GRCh38-GIABv3.bam"
BAM_OUTPUT="bam/HG008-N-D_PacBio_${CHR}.bam"
echo "Downloading region ${REGION}..."
samtools view -b "${BAM_URL}" ${REGION} > "${BAM_OUTPUT}"
echo "Indexing..."
samtools index "${BAM_OUTPUT}"
echo "Statistics:"
samtools flagstat "${BAM_OUTPUT}"
echo ""

echo "=========================================="
echo "Processing PacBio Cancer (Tumor)..."
BAM_URL="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_somatic/HG008/Liss_lab/PacBio_Onso_20240415/HG008-T_Pacbio-onso_136x_GRCh38-GIABv3.bam"
BAM_OUTPUT="bam/HG008-T_PacBio_${CHR}.bam"
echo "Downloading region ${REGION}..."
samtools view -b "${BAM_URL}" ${REGION} > "${BAM_OUTPUT}"
echo "Indexing..."
samtools index "${BAM_OUTPUT}"
echo "Statistics:"
samtools flagstat "${BAM_OUTPUT}"
echo ""

echo "=========================================="
echo "Processing Oxford Nanopore Cancer (Tumor)..."
BAM_URL="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_somatic/HG008/Liss_lab/Northeastern-ONT-UL-20241216/HG008-T_GRCh38_GIABv3_ONT-UL-R10.4.1-dorado_0.8.1_sup.5mC_5hmC_54x_20241216.bam"
BAM_OUTPUT="bam/HG008-T_ONT_${CHR}.bam"
echo "Downloading region ${REGION}..."
samtools view -b "${BAM_URL}" ${REGION} > "${BAM_OUTPUT}"
echo "Indexing..."
samtools index "${BAM_OUTPUT}"
echo "Statistics:"
samtools flagstat "${BAM_OUTPUT}"
echo ""

echo "=========================================="
echo "Download complete!"
echo ""
echo "All BAM files:"
ls -lh bam/*.bam
echo ""

# ===== ALIGNMENT COMPARISON ANALYSIS =====
echo "=========================================="
echo "DETAILED ALIGNMENT COMPARISON ANALYSIS"
echo "=========================================="
echo ""

REPORT_FILE="reports/alignment_comparison_report.txt"
echo "Alignment Comparison Report" > ${REPORT_FILE}
echo "Region: ${REGION}" >> ${REPORT_FILE}
echo "Date: $(date)" >> ${REPORT_FILE}
echo "========================================" >> ${REPORT_FILE}
echo "" >> ${REPORT_FILE}

# Function to analyze each BAM file
analyze_bam() {
    local bam_file=$1
    local sample_name=$(basename ${bam_file} .bam)
    
    echo "----------------------------------------"
    echo "Sample: ${sample_name}"
    echo "----------------------------------------"
    
    echo "" >> ${REPORT_FILE}
    echo "========== ${sample_name} ==========" >> ${REPORT_FILE}
    
    # 1. Basic alignment statistics
    echo "1. Basic Alignment Statistics:" | tee -a ${REPORT_FILE}
    samtools flagstat ${bam_file} | tee -a ${REPORT_FILE}
    echo "" >> ${REPORT_FILE}
    
    # 2. Read depth (coverage)
    echo "2. Read Depth Analysis:" | tee -a ${REPORT_FILE}
    samtools depth ${bam_file} > reports/${sample_name}_depth.txt
    
    # Calculate average, min, max depth
    awk '{sum+=$3; if(NR==1){min=$3; max=$3} if($3<min){min=$3} if($3>max){max=$3}} END {print "Average Depth: " sum/NR "\nMin Depth: " min "\nMax Depth: " max}' reports/${sample_name}_depth.txt | tee -a ${REPORT_FILE}
    echo "" >> ${REPORT_FILE}
    
    # 3. Mapping quality distribution
    echo "3. Mapping Quality Distribution:" | tee -a ${REPORT_FILE}
    samtools view ${bam_file} | awk '{print $5}' | sort -n | uniq -c | sort -rn | head -20 | tee -a ${REPORT_FILE}
    echo "" >> ${REPORT_FILE}
    
    # 4. Coverage uniformity (coefficient of variation)
    echo "4. Coverage Uniformity:" | tee -a ${REPORT_FILE}
    awk '{sum+=$3; sumsq+=$3*$3} END {mean=sum/NR; variance=sumsq/NR - mean*mean; sd=sqrt(variance); cv=sd/mean*100; print "Mean Coverage: " mean "\nStd Dev: " sd "\nCoefficient of Variation: " cv "%"}' reports/${sample_name}_depth.txt | tee -a ${REPORT_FILE}
    echo "" >> ${REPORT_FILE}
    
    # 5. Insert size statistics (for paired-end reads)
    echo "5. Insert Size Statistics:" | tee -a ${REPORT_FILE}
    samtools stats ${bam_file} | grep "^SN" | grep "insert size" | tee -a ${REPORT_FILE}
    echo "" >> ${REPORT_FILE}
    
    # 6. Read length distribution
    echo "6. Read Length Distribution:" | tee -a ${REPORT_FILE}
    samtools view ${bam_file} | awk '{print length($10)}' | sort -n | uniq -c | sort -rn | head -10 | tee -a ${REPORT_FILE}
    echo "" >> ${REPORT_FILE}
    
    echo ""
}

# Analyze all BAM files
for bamfile in bam/*.bam; do
    analyze_bam ${bamfile}
done

# ===== COMPARATIVE SUMMARY =====
echo "=========================================="
echo "COMPARATIVE SUMMARY"
echo "=========================================="

SUMMARY_FILE="reports/comparative_summary.txt"
echo "Comparative Summary Across All Samples" > ${SUMMARY_FILE}
echo "Region: ${REGION}" >> ${SUMMARY_FILE}
echo "Date: $(date)" >> ${SUMMARY_FILE}
echo "" >> ${SUMMARY_FILE}

# Create a table comparing key metrics
printf "%-30s %-15s %-15s %-15s %-15s\n" "Sample" "Total Reads" "Mapped %" "Avg Depth" "Avg MapQ" | tee -a ${SUMMARY_FILE}
printf "%-30s %-15s %-15s %-15s %-15s\n" "==============================" "===============" "===============" "===============" "===============" | tee -a ${SUMMARY_FILE}

for bamfile in bam/*.bam; do
    sample_name=$(basename ${bamfile} .bam)
    
    # Extract metrics
    total_reads=$(samtools view -c ${bamfile})
    mapped_pct=$(samtools flagstat ${bamfile} | grep "mapped (" | head -1 | awk '{print $5}' | tr -d '(%)' | tr -d ':')
    avg_depth=$(awk '{sum+=$3} END {print sum/NR}' reports/${sample_name}_depth.txt | xargs printf "%.2f")
    avg_mapq=$(samtools view ${bamfile} | awk '{sum+=$5; n++} END {print sum/n}' | xargs printf "%.2f")
    
    printf "%-30s %-15s %-15s %-15s %-15s\n" "${sample_name}" "${total_reads}" "${mapped_pct}" "${avg_depth}" "${avg_mapq}" | tee -a ${SUMMARY_FILE}
done

echo "" | tee -a ${SUMMARY_FILE}
echo "=========================================="
echo "Analysis complete!"
echo ""
echo "Reports saved to:"
echo "  - ${REPORT_FILE}"
echo "  - ${SUMMARY_FILE}"
echo "  - Individual depth files in reports/ directory"
echo ""

