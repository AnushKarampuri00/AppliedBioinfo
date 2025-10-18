# Assignment week 8 - Using a Design file and parallel command

This Makefile automates the process of downloading sequencing data, aligning it to a reference genome, generating fastqc reports and BAM, BigWig files . It's written for processing multiple samples from SRA datasets with controlled coverage depth using parallel command.

## Makefile tasks


| Task                | Description                                                                 |
|----------------------|------------------------------------------------------------------------------|
| calculate_coverage   | Determines how many reads to download based on genome size and desired coverage |
| download_reads       | Downloads only the required number of reads from SRA                        |
| fastqc               | Generates quality control reports for the downloaded reads                  |
| genome               | Fetches the reference genome sequence                                       |
| index                | Creates index files for the reference genome                                |
| align                | Aligns reads to the reference and creates sorted BAM files                  |
| stats                | Generates alignment statistics and metrics                                  |
| bigwig               | Creates BigWig coverage tracks from aligned reads                           |



## Usage
Basic Usage
The makefile is written to be run using GNU parallel with a design CSV file:

```bash
parallel --jobs 1 --colsep , --header : \
  make -f makefile.mk SRR={run_accession} SAMPLE={sample_name} ACC={ACC} COVERAGE={COVERAGE} < design.csv
```

your design.csv file should have columns like:

```text

run_accession,sample_name,ACC,COVERAGE
SRR123456,Sample1,NC_007793.1,15
SRR789012,Sample2,NC_007793.1,20

```

## Dry Run

To see what commands would be executed without actually running them:

```bash
parallel --dry-run --jobs 1 --colsep , --header : \
  make -f makefile.mk SRR={run_accession} SAMPLE={sample_name} ACC={ACC} COVERAGE={COVERAGE} < design.csv
```

## Manual Single Sample Processing

You can also process a single sample manually:

```bash
make -f makefile2.mk SRR=SRR123456 SAMPLE=MySample ACC=NC_007793.1 COVERAGE=15
```

## Output Files:

| Output File        | Location | Description                                                |
|--------------------|-----------|------------------------------------------------------------|
| Raw reads          | reads/    | Downloaded FASTQ files (forward and reverse if paired-end) |
| Quality reports    | fastqc/   | HTML and ZIP files with read quality metrics              |
| Reference genome   | refs/     | FASTA file of the reference sequence                      |
| Aligned reads      | bam/      | Sorted BAM file and its index                             |
| Alignment stats    | bam/      | Text file with alignment statistics                       |
| Coverage track     | bigwig/   | BigWig file for visualization in genome browsers           |
| Temp files         | temp/     | Intermediate files and processing logs                    |


## Requirements

Make sure you have these tools installed and in your PATH:

* bwa (for alignment)
* samtools (for BAM processing)
* bedtools (for coverage calculations)
* bedGraphToBigWig (from UCSC tools)
* fastq-dump (from SRA Toolkit)
* efetch (from NCBI E-utilities)
* fastqc (for quality control)
* GNU parallel


## Cleanup

To remove all generated files and start fresh:

```bash
make -f makefile.mk clean
```




