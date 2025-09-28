# Sequencing data download and quality assessment

## Requeirements
* Please install the required dependencies before performing the data retreival and quality assessment.

``` bash

conda create -n sra_env -c bioconda -c conda-forge \
    entrez-direct \
    sra-tools \
    fastqc \
    coreutils

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

I have identified another SRA number - SRR34913726. This is from oxford nanopore sequncing method.




