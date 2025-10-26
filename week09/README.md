# Week09 - Revised automation code

This readme file provides the code using parallel command for automating the sequence alignment for a set of samples.

## Usage

```bash

# To check the targets and the usage of this make file

make -f makefile.mk usage

# To download the genome and index it.
# This step only needs to be run once since the same genome is used for all samples.

make -f makefile.mk genome index

# To run the remaining important tasks other than downloading the genome and indexing, the below input is used.
# This runs "process_sample", which performs counting reads, aligning, BAM and BigWig file generation (except reference download and indexing tasks) using parallel command to iterate over the samples present in the design.csv file.
 
cat design.csv | parallel --jobs 3 --colsep , --header : \ make -f makefile.mk process_sample SRR={SRR} SAMPLE={name} COVERAGE={coverage}

# Preview the commands without executing them (dry run).
# Use this to verify which commands would be executed before running the full pipeline.

cat design.csv | parallel --jobs 3 --colsep , --header : \ make --dry-run -f makefile.mk process_sample SRR={SRR} SAMPLE={name} COVERAGE={coverage}


```

## Design file
The design file contains three columns:
1. SRR numbers
2. Sample Names
3. Coverage required (I included this to have the code calculate the required reads for different coverages for all the samples)



## Output

This makefile generates the below given directories and saves the files accordingly in them.

```
  |-temp
  |-bam
  |-reads
  |-refs
  |-fastqc
  |-bigwig
```

It took me around 20 seconds to finish running this code
```bash

[main] Version: 0.7.19-r1273
[main] CMD: bwa mem refs/NC_007793.1.fa reads/Sample8_1.fastq reads/Sample8_2.fastq
[main] Real time: 20.306 sec; CPU: 20.451 sec

```

## To run single sample

```bash

# use "all" for running a single sample, this way all the tasks from reference downloading till bigwig file generation will run sequentially.

make -f makefile.mk all SRR=SRR35862149 SAMPLE=S1 COVERAGE=15

```






