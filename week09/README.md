# Week09 - Revised automation code

This readme file provides the code using parallel command for automating the sequence alignment for a set of samples.

## Usage

```bash
# To download the genome and index it.
# This step only needs to be run once since the same genome is used for all samples.

make -f makefile.mk genome index

# To run the remaining important tasks other than downloading the genome and indexing, the below input is used.
# This runs "process_sample", which performs counting reads, aligning, BAM and BigWig file generation using parallel command to iterate over the samples present in the design.csv file.
 
cat design1.csv | parallel --jobs 3 --colsep , --header : \ make -f week9_trial2.mk process_sample SRR={SRR} SAMPLE={name} COVERAGE={coverage}

# Preview the commands without executing them (dry run).
# Use this to verify which commands would be executed before running the full pipeline.

cat design1.csv | parallel --jobs 3 --colsep , --header : \ make --dry-run -f week9_trial2.mk process_sample SRR={SRR} SAMPLE={name} COVERAGE={coverage}


```
