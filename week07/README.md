# Make file - Week 07 - Bigwig file generation

This readme file contains the information needed to run the makefile generated. The script written is able to It automatically calculates the exact number of reads needed to achieve your desired sequencing coverage, then downloads and processes the data through alignment and BigWig file generation.

## Basic usage

```bash

make -f bound_script2.mk all SRR=SRRXXXXXXX ACC=NC_XXXXXX.1 COVERAGE=15

```

## Example

```bash

make -f bound_script2.mk all SRR=ERR15403387 ACC=NC_007793.1 COVERAGE=15

```

## Required Parameters

Parameter	Description	Example
SRR	SRA Run accession number	SRR=ERR15403387
ACC	NCBI Reference genome accession	ACC=NC_007793.1
COVERAGE	Desired sequencing coverage	COVERAGE=15

## Output Files

Output directories:

```bash

project/
├── reads/           # Downloaded FASTQ files
├── refs/           # Reference genome files
├── bam/            # Aligned BAM files + indexes
└── bigwig/         # BigWig coverage files

```


