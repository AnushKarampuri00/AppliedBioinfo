# Week 1 Assignment

This document explains the steps I followed for the assignment.

---

## 1. Create a new directory for this week's assignment
I created a directory named `week01` inside my cloned repository.

```bash
mkdir week01
cd week01
```

## 2. Add a README.md file to the directory

I created this README.md file using:
```
touch README.md
```

## 3. Check samtools version in the bioinfo environment

After activating the bioinfo environment, I checked the version of samtools:

```
conda activate bioinfo
samtools --version
```

* Sam tools version iam using is:
```
samtools 1.22.1
Using htslib 1.22.1
Copyright (C) 2025 Genome Research Ltd.
```


## 4. Commands to create a nested directory structure

I created nested directories (dir1/dir2/dir3) as follows:
- **Week01(Assignment01)/** → Each assignment has its own folder.  
- **Data/** → Contains input datasets used for the assignment.  
- **Results/** → Output files, plots, and reports generated from the scripts.

```
mkdir -p data/raw data results/figures results/tables
```

## 5. Commands that create files in different directories

I created some files at different levels of the nested structure:

```
echo "sample1 ACGTACGT" > data/raw/sample1.txt
echo "sample2 GGGGAAAA" > data/raw/sample2.txt
touch results/tables/summary.tsv
```

## 6. Access files using relative and absolute paths

Relative path :
```
cat data/raw/sample1.txt
```
Absolute path :
```
/home/anushkarampuri/projects/AppliedBioinfo_Week1/week01/data/raw/sample1.txt
```

## 7. Commit and push changes to the repository
```
git add .
git commit -m "Week1 assignment with README and directories"
git push origin main
```

## 8. Repository Link:

https://github.com/AnushKarampuri00/AppliedBioinfo




 
 





