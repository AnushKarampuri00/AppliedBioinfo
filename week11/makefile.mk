# USAGE : cat design.csv | parallel --jobs 3 --colsep , --header : --eta --bar --verbose \
# make -f makefile_2.mk process_sample SRR={SRR} SAMPLE={name} COVERAGE={coverage}

SHELL := bash
.ONESHELL:
.SHELLFLAGS := -eu -o pipefail -c
.DELETE_ON_ERROR:
MAKEFLAGS += --warn-undefined-variables --no-builtin-rules

# Parameters
ACCESSION := NC_007793.1
CHROMOSOME := Chromosome
DESIRED_COVERAGE := $(COVERAGE)
SNPEFF_DB := Staphylococcus_aureus_subsp_aureus_usa300_fpr3757_gca_000013465

# Directories
READ_DIR := reads
REF_DIR := refs
BAM_DIR := bam
VCF_DIR := vcf
BIGWIG_DIR := bigwig
FASTQC_DIR := fastqc
TEMP_DIR := temp
ANNOTATION_DIR := vcf_annotations

# Merged VCF Output files
MERGED_VCF := $(VCF_DIR)/merged_multisample.vcf.gz

# File paths
REF_GENOME := $(REF_DIR)/$(ACCESSION).fa

# Create directories
$(shell mkdir -p $(READ_DIR) $(REF_DIR) $(BAM_DIR) $(VCF_DIR) $(BIGWIG_DIR) $(FASTQC_DIR) $(TEMP_DIR) $(ANNOTATION_DIR))

# Master target (for single sample)
all: genome index process_sample
	@echo "✅ Full pipeline completed for sample $(SAMPLE) ($(SRR)) with $(COVERAGE)X coverage"

# Download reference genome (once)
genome:
	@echo "📥 Fetching reference genome..."
	@if [ ! -f "$(REF_GENOME)" ]; then \
		efetch -db nucleotide -id $(ACCESSION) -format fasta > $(REF_GENOME); \
		echo "Reference saved: $(REF_GENOME)"; \
	else \
		echo "Reference genome already exists. Skipping download."; \
	fi

# Index reference genome (once)
index: $(REF_GENOME)
	@echo "🧩 Indexing reference..."
	@if [ ! -f "$(REF_GENOME).bwt" ]; then \
		bwa index $(REF_GENOME); \
		samtools faidx $(REF_GENOME); \
		echo "Indexing done."; \
	else \
		echo "Index files already exist. Skipping indexing."; \
	fi

# Per-sample workflow (removed merge_vcfs from here)
process_sample: calculate_coverage download_reads fastqc align stats vcf bigwig vcf_annotation
	@echo "✅ Pipeline completed for $(SAMPLE) ($(SRR)) with $(COVERAGE)X coverage"
	@echo "$(SAMPLE)" >> $(TEMP_DIR)/processed_samples.txt

# Calculate coverage
calculate_coverage:
	@echo "🧮 Calculating reads needed for $(COVERAGE)X coverage..."
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
	@echo "📥 Downloading reads for $(SAMPLE)..."
	@REQUIRED_READS=$$(cat $(TEMP_DIR)/$(SAMPLE)_required_reads.txt); \
	mkdir -p $(READ_DIR); \
	fastq-dump -X $$REQUIRED_READS --outdir $(READ_DIR) --split-files $(SRR)
	@mv $(READ_DIR)/$(SRR)_1.fastq $(READ_DIR)/$(SAMPLE)_1.fastq || true
	@mv $(READ_DIR)/$(SRR)_2.fastq $(READ_DIR)/$(SAMPLE)_2.fastq || true
	@ls -lh $(READ_DIR)

# FASTQC
fastqc:
	@echo "🔬 Running FastQC for $(SAMPLE)..."
	@if [ -f "$(READ_DIR)/$(SAMPLE)_1.fastq" ]; then fastqc $(READ_DIR)/$(SAMPLE)_1.fastq -o $(FASTQC_DIR); fi
	@if [ -f "$(READ_DIR)/$(SAMPLE)_2.fastq" ]; then fastqc $(READ_DIR)/$(SAMPLE)_2.fastq -o $(FASTQC_DIR); fi

# Align reads
align:
	@echo "🧬 Aligning reads for $(SAMPLE)..."
	@if [ -f "$(READ_DIR)/$(SAMPLE)_1.fastq" ] && [ -f "$(READ_DIR)/$(SAMPLE)_2.fastq" ]; then \
		bwa mem $(REF_GENOME) $(READ_DIR)/$(SAMPLE)_1.fastq $(READ_DIR)/$(SAMPLE)_2.fastq | samtools view -bS - | samtools sort -o $(BAM_DIR)/$(SAMPLE).sorted.bam; \
		echo "Your FASTQ data was paired-end sequenced."; \
	elif [ -f "$(READ_DIR)/$(SAMPLE)_1.fastq" ]; then \
		bwa mem $(REF_GENOME) $(READ_DIR)/$(SAMPLE)_1.fastq | samtools view -bS - | samtools sort -o $(BAM_DIR)/$(SAMPLE).sorted.bam; \
		echo "Your FASTQ data was single-end sequenced."; \
	else \
		echo "No FASTQ files found for $(SAMPLE)"; \
		exit 1; \
	fi
	samtools index $(BAM_DIR)/$(SAMPLE).sorted.bam

# Alignment statistics
stats:
	samtools flagstat $(BAM_DIR)/$(SAMPLE).sorted.bam > $(BAM_DIR)/$(SAMPLE)_alignment_stats.txt
	@echo "📊 Stats saved to $(BAM_DIR)/$(SAMPLE)_alignment_stats.txt"

# Variant calling
vcf:
	@echo "🧫 Calling variants for $(SAMPLE)..."
	@if [ -f "$(BAM_DIR)/$(SAMPLE).sorted.bam" ]; then \
		bcftools mpileup -Ou -f $(REF_GENOME) $(BAM_DIR)/$(SAMPLE).sorted.bam | \
		bcftools call -mv -Oz -o $(VCF_DIR)/$(SAMPLE).vcf.gz; \
		bcftools index $(VCF_DIR)/$(SAMPLE).vcf.gz; \
		echo "✅ VCF generated: $(VCF_DIR)/$(SAMPLE).vcf.gz"; \
	else \
		echo "❌ BAM file not found for $(SAMPLE). Skipping VCF generation."; \
	fi

# VCF Annotation for individual samples
vcf_annotation:
	@echo "🧬 Annotating individual VCF for $(SAMPLE)..."
	@mkdir -p $(ANNOTATION_DIR)/individual/$(SAMPLE)
	@if [ -f "$(VCF_DIR)/$(SAMPLE).vcf.gz" ]; then \
		echo "Processing $(SAMPLE).vcf.gz..."; \
		gunzip -c $(VCF_DIR)/$(SAMPLE).vcf.gz | sed 's/^$(ACCESSION)/$(CHROMOSOME)/g' > $(TEMP_DIR)/$(SAMPLE)_fixed.vcf; \
		snpEff -c snpEff.config $(SNPEFF_DB) $(TEMP_DIR)/$(SAMPLE)_fixed.vcf > $(ANNOTATION_DIR)/individual/$(SAMPLE)/$(SAMPLE)_annotated.vcf 2>&1; \
		if [ -f "snpEff_summary.html" ]; then mv snpEff_summary.html $(ANNOTATION_DIR)/individual/$(SAMPLE)/$(SAMPLE)_snpEff_summary.html; fi; \
		if [ -f "snpEff_genes.txt" ]; then mv snpEff_genes.txt $(ANNOTATION_DIR)/individual/$(SAMPLE)/$(SAMPLE)_snpEff_genes.txt; fi; \
		rm $(TEMP_DIR)/$(SAMPLE)_fixed.vcf; \
		echo "✅ Annotation completed for $(SAMPLE)"; \
		echo "   - Annotated VCF: $(ANNOTATION_DIR)/individual/$(SAMPLE)/$(SAMPLE)_annotated.vcf"; \
		echo "   - Summary HTML: $(ANNOTATION_DIR)/individual/$(SAMPLE)/$(SAMPLE)_snpEff_summary.html"; \
		echo "   - Genes TXT: $(ANNOTATION_DIR)/individual/$(SAMPLE)/$(SAMPLE)_snpEff_genes.txt"; \
	else \
		echo "❌ VCF file not found for $(SAMPLE). Skipping annotation."; \
	fi

# BigWig generation
bigwig:
	@echo "📈 Creating BigWig for $(SAMPLE)..."
	bedtools genomecov -ibam $(BAM_DIR)/$(SAMPLE).sorted.bam -split -bg | sort -k1,1 -k2,2n > $(BIGWIG_DIR)/$(SAMPLE).bedgraph
	bedGraphToBigWig $(BIGWIG_DIR)/$(SAMPLE).bedgraph $(REF_GENOME).fai $(BIGWIG_DIR)/$(SAMPLE).bw
	@echo "BigWig created: $(BIGWIG_DIR)/$(SAMPLE).bw"

# Merge VCFs (run this AFTER all samples are processed)
merge_vcfs:
	@set -euo pipefail; \
	echo "🔍 Searching for individual sample .vcf.gz files in '$(VCF_DIR)'..."; \
	VCF_FILES=$$(find "$(VCF_DIR)" -maxdepth 1 -type f -name "*.vcf.gz" ! -name "merged_*.vcf.gz" | sort); \
	if [ -z "$$VCF_FILES" ]; then \
		echo "❌ No individual VCF files found in $(VCF_DIR)/"; \
		exit 1; \
	fi; \
	VCF_COUNT=$$(echo "$$VCF_FILES" | wc -l); \
	echo "✅ Found $$VCF_COUNT VCF file(s):"; \
	echo "$$VCF_FILES" | sed 's/^/   • /'; \
	if [ "$$VCF_COUNT" -eq 1 ]; then \
		echo "⚠️  Only one VCF file found. Skipping merge (need at least 2 files to merge)."; \
		echo "   If you want to annotate this single file, it's already annotated in $(ANNOTATION_DIR)/individual/"; \
		exit 0; \
	fi; \
	echo "🧩 Checking and indexing VCF files if necessary..."; \
	for vcf in $$VCF_FILES; do \
		if [ ! -f "$${vcf}.csi" ]; then \
			echo "   Indexing $$vcf..."; \
			bcftools index --csi "$$vcf"; \
		else \
			echo "   $$vcf already indexed."; \
		fi; \
	done; \
	echo "🧬 Merging $$VCF_COUNT VCFs into a single multi-sample file..."; \
	bcftools merge $$VCF_FILES -O z -o "$(MERGED_VCF)"; \
	echo "📦 Indexing merged file with CSI..."; \
	bcftools index --csi "$(MERGED_VCF)"; \
	echo "✅ Multi-sample VCF generated successfully:"; \
	echo "   → $(MERGED_VCF)"; \
	echo "   → $(MERGED_VCF).csi"

# Annotate merged VCF (run this AFTER merge_vcfs)
annotate_merged:
	@echo "🧬 Annotating merged multi-sample VCF..."
	@mkdir -p $(ANNOTATION_DIR)/merged
	@if [ ! -f "$(MERGED_VCF)" ]; then \
		echo "❌ Merged VCF file not found: $(MERGED_VCF)"; \
		echo "   Please run 'make merge_vcfs' first."; \
		exit 1; \
	fi
	@echo "Processing merged_multisample.vcf.gz..."
	@gunzip -c $(MERGED_VCF) | sed 's/^$(ACCESSION)/$(CHROMOSOME)/g' > $(TEMP_DIR)/merged_multisample_fixed.vcf
	@snpEff -c snpEff.config $(SNPEFF_DB) $(TEMP_DIR)/merged_multisample_fixed.vcf > $(ANNOTATION_DIR)/merged/merged_multisample_annotated.vcf 2>&1
	@if [ -f "snpEff_summary.html" ]; then mv snpEff_summary.html $(ANNOTATION_DIR)/merged/merged_multisample_snpEff_summary.html; fi
	@if [ -f "snpEff_genes.txt" ]; then mv snpEff_genes.txt $(ANNOTATION_DIR)/merged/merged_multisample_snpEff_genes.txt; fi
	@rm $(TEMP_DIR)/merged_multisample_fixed.vcf
	@echo "✅ Merged VCF annotation completed"
	@echo "   - Annotated VCF: $(ANNOTATION_DIR)/merged/merged_multisample_annotated.vcf"
	@echo "   - Summary HTML: $(ANNOTATION_DIR)/merged/merged_multisample_snpEff_summary.html"
	@echo "   - Genes TXT: $(ANNOTATION_DIR)/merged/merged_multisample_snpEff_genes.txt"

# Complete workflow: merge and annotate (run after all samples)
finalize: merge_vcfs annotate_merged
	@echo "✅ All VCF files merged and annotated successfully!"
	@echo ""
	@echo "📂 Output structure:"
	@echo "   Individual annotations: $(ANNOTATION_DIR)/individual/<sample_name>/"
	@echo "   Merged VCF: $(MERGED_VCF)"
	@echo "   Merged annotation: $(ANNOTATION_DIR)/merged/"

# Clean up
clean:
	rm -rf $(READ_DIR) $(REF_DIR) $(BAM_DIR) $(VCF_DIR) $(BIGWIG_DIR) $(FASTQC_DIR) $(TEMP_DIR) $(ANNOTATION_DIR)
	@echo "🧹 Cleaned all files"

.PHONY: all genome index process_sample calculate_coverage download_reads fastqc align stats vcf vcf_annotation bigwig merge_vcfs annotate_merged finalize clean

