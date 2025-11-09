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

Primary task is to download the genome and index it

```bash
make -f makefile.mk genome index
```

To work with a single sample 

```bash
make -f makefile.mk process_sample SRR=SRR35862149 SAMPLE=S1 COVERAGE=17
```

To work on multiple samples

```bash
cat design.csv | parallel --jobs 3 --colsep , --header : --eta --bar --verbose \
  make -f makefile.mk process_sample SRR={SRR} SAMPLE={name} COVERAGE={coverage}


# To merge the VCF files - if you work with multiple samples this step is needed

make -f makefile.mk finalize
```


Description of tasks from the makefile

| Task                | Description                                                                 |
|----------------------|------------------------------------------------------------------------------|
| calculate_coverage   | Determines how many reads to download based on genome size and desired coverage |
| download_reads       | Downloads only the required number of reads from SRA                        |
| fastqc               | Generates quality control reports for the downloaded reads                  |
| genome               | Fetches the reference genome sequence                                       |
| index                | Creates index files for the reference genome                                |
| process_sample       | Complete pipeline for individual sample                                     |
| align                | Aligns reads to the reference and creates sorted BAM files                  |
| stats                | Generates alignment statistics and metrics                                  |
| bigwig               | Creates BigWig coverage tracks from aligned reads                           |
| vcf                  | Calling variants                                                            |
| merge_vcfs           | merge the generated vcf files                                               |
| vcf_annotation       | performs annotation using snpeff tool                                       |
| finalize             | This performs merging vcf and annotating the merged vcf file                |


## Variants identified


**Variant 1: Position 2,865,224 T>A**

*  **Basic Change:** T to A substitution at position 2,865,224.

* **Genotype:** Homozygous Alternate (1/1) in all samples except Sample8 (missing).

* **Impact:** Missense Variant (MODERATE) — results in an amino acid change from Glutamic acid (Glu) to Aspartic acid (Asp) at position 52 (p.Glu52Asp) in the coding region.

**Key Affected Genes & Region Context:**

* **Primary gene:** ENSB:2zQP9zbD3cehAHU - it is a protein-coding gene

**Additional context:** Multiple downstream and upstream variants noted across neighboring genes such as bceB_3, cspLA, and rsmG, mostly MODIFIER effects.

* **Quality:** High (225.417); DP=104, with 100+ reads supporting the alternate allele.

**Variant 2: Position 612,319 C>T**

* **Basic Change:** C to T substitution at position 612,319.

* **Genotype:** Homozygous Alternate (1/1) in all samples except Sample8 (missing).

* **Impact:** Synonymous Variant (LOW) — does not alter the amino acid sequence (p.Phe391Phe) in the sdrC gene.

**Key Affected Genes & Region Context:**

* **Primary gene:** sdrC — synonymous change in the coding region.

* **Nearby genes:** Upstream variants in dck, dgk, and sdrD; downstream variants in tadA_1, ywpJ_1, and azo1, all MODIFIER effects.

* **Quality:** High (225.417); DP=129, well-supported by multiple reads.

**Variant 3: Position 612,864 T>C**

* **Basic Change:** T to C substitution at position 612,864.

* **Genotype:** Homozygous Alternate (1/1) across all eight samples.

* **Impact:** Missense Variant (MODERATE) — amino acid substitution Valine (Val) → Alanine (Ala) at position 573 (p.Val573Ala) within sdrC.

* **Key Affected Genes & Region Context:**

* **Primary gene:** sdrC — moderate effect within same gene as Variant 2.

* **Additional context:** Upstream and downstream effects in genes like dck, dgk, and tadA_1 (all MODIFIER).

* **Quality:** High (225.417); DP=82 with strong alternate allele support.



## HTML Summary file

![HTML Summary file](images/summary_1.png)



![HTML Summary file](images/summary_2.png)



![HTML Summary file](images/summary_3.png)





