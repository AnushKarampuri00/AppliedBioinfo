# ================================
# Makefile for merging VCF files
# ================================

# Variables
VCF_DIR := vcf
OUTPUT := $(VCF_DIR)/merged_multisample.vcf.gz

# Default target
all: merge_vcfs

# Step 1: Merge all VCFs and index the merged file
merge_vcfs:
	@set -euo pipefail; \
	echo "🔍 Searching for .vcf.gz files in '$(VCF_DIR)'..."; \
	VCF_FILES=$$(find "$(VCF_DIR)" -type f -name "*.vcf.gz" | sort); \
	if [ -z "$$VCF_FILES" ]; then \
		echo "❌ No VCF files found in $(VCF_DIR)/"; \
		exit 1; \
	fi; \
	echo "✅ Found the following VCF files:"; \
	echo "$$VCF_FILES" | sed 's/^/   • /'; \
	echo "🧩 Checking and indexing VCF files if necessary..."; \
	for vcf in $$VCF_FILES; do \
		if [ ! -f "$${vcf}.csi" ]; then \
			echo "   Indexing $$vcf..."; \
			bcftools index --csi "$$vcf"; \
		else \
			echo "   $$vcf already indexed."; \
		fi; \
	done; \
	echo "🧬 Merging VCFs into a single multi-sample file..."; \
	bcftools merge $$VCF_FILES -O z -o "$(OUTPUT)"; \
	echo "📦 Indexing merged file with CSI..."; \
	bcftools index --csi "$(OUTPUT)"; \
	echo "✅ Multi-sample VCF generated successfully:"; \
	echo "   → $(OUTPUT)"; \
	echo "   → $(OUTPUT).csi"

# Step 2: Clean up merged files if needed
clean:
	@echo "🧹 Removing merged files..."
	rm -f $(OUTPUT) $(OUTPUT).csi
	@echo "✅ Clean complete."
	
.PHONY: all merge_vcfs clean