# Variant effect evaluation - Week11

This Makefile is an extension of the previous assignment by including a new target named "vcf_annotation". This makefile can perform reads alignment, variant calling and variant annotation.

## Requirement
1. Installing the required package
* I used snpeff package for this assignment.

```bash
# install SnpEff using conda
conda install -c bioconda snpeff
```

2. Download the reference database for the organism
* I downloaded the database for Staphylococcus aureus subsp. aureus USA300_FPR3757 based on my accession ID : NC_007793
  
```bash

# check available SnpEff databases for Staphylococcus aureus USA300
snpEff databases | grep -i "aureus" | grep -i "usa300"

# copy the snpEff.config file to the current directory
# This step will enable us to download the files required for vcf annotation in the current directory.
# Otherwise it downloads in the default path where the snpEFF was downloaded
# This step just makes things easy to locate, you can skip this and check in the default path as well for download.

cp /opt/anaconda3/envs/bioinformatics/share/snpeff-5.3.0a-1/snpEff.config ./snpEff.config

# to download the SnpEff database for Staphylococcus aureus USA300_FPR3757
snpEff download Staphylococcus_aureus_subsp_aureus_usa300_fpr3757_gca_000013465

# list the files in the SnpEff data directory for Staphylococcus aureus USA300_FPR3757

ls -l ./data/Staphylococcus_aureus_subsp_aureus_usa300_fpr3757_gca_000013465/
# You should find these files within the created directory - sequence.Chromosome.bin, sequence.bin, snpEffectPredictor.bin


```


## Usage

To work with single sample

```bash
make -f makefile.mk all  SRR=SRR35862149 SAMPLE=S1 COVERAGE=17
```

To work on multiple samples

```bash
cat design.csv | parallel --jobs 3 --colsep , --header : --eta --bar --verbose \ make -f makefile.mk process_sample SRR={SRR} SAMPLE={name} COVERAGE={coverage}
```

## Variants identified




