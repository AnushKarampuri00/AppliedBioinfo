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

## 2. 
