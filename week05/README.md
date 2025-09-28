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
bash data_download.sh
```
* **Note:** I calculated the required number of reads by dividing the total bases needed for 10× coverage by the average read length (150 bp). Then I used fasterq-dump with --row-limit to download only that subset instead of the full dataset.

## 2. 

ststistics output

```
TOTAL_READS=4152
TOTAL_BASES=377832
AVG_READ_LENGTH=91
GENOME_SIZE=10389
ACTUAL_COVERAGE=36.36
Target coverage: 10x
```



