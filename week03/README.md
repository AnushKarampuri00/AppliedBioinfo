# Visualizing Genomic Data

## Download Data
* I downloaded the **Acetonema longum DSM 6540** genome FASTA and GFF annotation files from Ensembl Bacteria


```bash

#FASTA Files
curl -O https://ftp.ensemblgenomes.ebi.ac.uk/pub/bacteria/current/fasta/bacteria_1_collection/acetonema_longum_dsm_6540_gca_000219125/dna/Acetonema_longum_dsm_6540_gca_000219125.ASM21912v1_.dna.toplevel.fa.gz
curl -O https://ftp.ensemblgenomes.ebi.ac.uk/pub/bacteria/current/gff3/bacteria_1_collection/acetonema_longum_dsm_6540_gca_000219125/Acetonema_longum_dsm_6540_gca_000219125.ASM21912v1.62.gff3.gz

# Annotation files
gunzip Acetonema_longum_dsm_6540_gca_000219125.ASM21912v1_.dna.toplevel.fa.gz
gunzip Acetonema_longum_dsm_6540_gca_000219125.ASM21912v1.62.gff3.gz 

```

## 1. How big is the genome, and how many features of each type does the GFF file contain?

```bash

grep -v ">" Acetonema_longum_dsm_6540_gca_000219125.ASM21912v1_.dna.toplevel.fa | wc -m

```

My Output:

```
The Genome is 4395200 basepair long

```
**Feature Counts in GFF:**

``` bash
cut -f3 Acetonema_longum_dsm_6540_gca_000219125.ASM21912v1.62.gff3 | sort | uniq -c
```

 My Output:
 ```
 4047 CDS
   4107 exon
   4047 gene
   4047 mRNA
     60 ncRNA
     60 ncRNA_gene
    296 region
```


## 2. From your GFF file, separate the intervals of type "gene" or "transcript" into a different file.

```bash
grep "^#" Acetonema_longum_dsm_6540_gca_000219125.ASM21912v1.62.gff3 > simplified.gff3

grep -w -E "gene|transcript" Acetonema_longum_dsm_6540_gca_000219125.ASM21912v1.62.gff3 >> simplified.gf
f3

```


# IGV Genome Visualization

![Use IGV to visualize your genome and the annotations relative to the genome.](images/Acetpnema_longum1.png)



