# Sequencing data download and quality assessment

## Requeirements
* Please install the required dependencies before performing the data retreival and quality assessment.

``` bash

# Create a new conda environment for bioinformatics
conda create -n bioinfo python=3.9 -y

# Activate the environment
conda activate bioinfo

# Install all required packages
conda install -c bioconda -y \
    sra-tools \
    fastqc \
    entrez-direct \
    multiqc \
    seqtk
```

## 1. Downloading Data

* Step-1: Retreiving SRR numbers from the accession number
* Step-2: Downloading the FASTQ file from the SRR numbers obtained
Please run the given script - "data_download.sh" script in your working directory to retreive the data corresponding to a SRR number from the list.

```bash
bash download_and_process.sh
```

* **Note:** I calculated the required number of reads by dividing the total bases needed for 10× coverage by the average read length (150 bp). Then I used fasterq-dump with --row-limit to download only that subset instead of the full dataset.


## 2. Quality assessment

Please run the given script - "quality_assessment.sh" to perfrom basic statistics calculations. The outputs we obtain from the script are given below

```bash
bash quality_assessment.sh
```
My output:

```
TOTAL_READS=4152
TOTAL_BASES=377832
AVG_READ_LENGTH=91
GENOME_SIZE=10389
ACTUAL_COVERAGE=36.36
Target coverage: 10x


PASS    Basic Statistics        SRR32692327_3.fastq
PASS    Per base sequence quality       SRR32692327_3.fastq
PASS    Per tile sequence quality       SRR32692327_3.fastq
PASS    Per sequence quality scores     SRR32692327_3.fastq
WARN    Per base sequence content       SRR32692327_3.fastq
PASS    Per sequence GC content SRR32692327_3.fastq
PASS    Per base N content      SRR32692327_3.fastq
PASS    Sequence Length Distribution    SRR32692327_3.fastq
FAIL    Sequence Duplication Levels     SRR32692327_3.fastq
PASS    Overrepresented sequences       SRR32692327_3.fastq
PASS    Adapter Content SRR32692327_3.fastq


#Measure        Value
Filename        SRR32692327_3.fastq
File type       Conventional base calls
Encoding        Sanger / Illumina 1.9
Total Sequences 14764041
Total Bases     1.3 Gbp
Sequences flagged as poor quality       0
Sequence length 91
%GC     46

```

**Interpretation:** Based on my calculations, ~692 reads (~10× coverage) were needed for a genome of 10,389 bp at 10× depth. After downloading a subset using fasterq-dump, I obtained 4,152 reads totaling ~378 Kbp, which corresponds to ~36× coverage—exceeding the target. Quality assessment with FastQC showed generally high-quality reads (PASS in most metrics), with only minor biases in base composition and a failure in duplication levels, likely due to the subset size and technical redundancy. Overall, the data is of sufficient quality and depth for downstream genomic analyses.

``` calculation:

# Calculation:
REQUIRED_BASES = GENOME_SIZE × COVERAGE
               = 10,389 × 10
               = 103,890 bases
AVG_READ_LEN = 150 bp
REQUIRED_READS = REQUIRED_BASES ÷ AVG_READ_LEN
               = 103,890 ÷ 150
               = 692 reads (rounded down)
```

## 3. Comparing sequencing platforms:

* I have identified another SRA ID - **SRR34913726**. This is from oxford nanopore sequncing method.

The quality metrics calculated are:

```bash output

TOTAL_READS=60
TOTAL_BASES=101065
AVG_READ_LENGTH=1684.42
GENOME_SIZE=10389
ACTUAL_COVERAGE=9.72

```

Breif comparision
* The Illumina dataset (SRR32692327) contained thousands of short reads (~91 bp) with high coverage (~36×) and mostly passed all quality metrics, showing only minor warnings in per-base sequence content and some duplication. In contrast, the Oxford Nanopore dataset (SRR34913726) had very few long reads (~1.68 kb) with slightly lower coverage (~9.7×), and several metrics, including per-sequence quality, per-base content, and GC content, showed warnings or errors. Illumina reads are highly uniform and accurate, ideal for detecting small variants, while Nanopore reads are much longer, capturing structural information but with higher error rates and overrepresented sequences. Overall, the two technologies complement each other: Illumina provides depth and precision, whereas Nanopore provides long-range context.
