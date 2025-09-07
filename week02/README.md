# GFF3 File Analysis using UNIX Commands

This repository demonstrates the analysis of the **Amphilophus citrinellus** (Midas cichlid) genome GFF3 file using **UNIX command-line tools**.

---

## Download the GFF3 File

The GFF3 file used for this analysis is **not included** due to its large size. You can download it directly from Ensembl using:

```bash
wget ftp://ftp.ensembl.org/pub/current_gff3/Amphilophus_citrinellus/Amphilophus_citrinellus.Midas_v5.115.gff3.gz
```

or

```bash
curl -O ftp://ftp.ensembl.org/pub/current_gff3/Amphilophus_citrinellus/Amphilophus_citrinellus.Midas_v5.115.gff3.gz
```


## 1. Tell us about the organism
_Amphilophus citrinellus_ is a large cichlid fish native to the San Juan River and surrounding watersheds in Costa Rica and Nicaragua. In genomics and aquaculture contexts, it is commonly referred to as the **Midas cichlid**. The species has a diploid chromosome number of 2n = 48, and its genome exhibits a high degree of gene duplication and structural complexity typical of cichlid fishes. This genome assembly provides a framework for studying gene regulation, transcript diversity, and evolutionary dynamics within the Cichlidae family.

## 2. How many sequence regions (chromosomes) does the file contain? 
* In GFF3, sequence regions are in the header lines starting with ##sequence-region.

```
zgrep "##sequence-region" Amphilophus_citrinellus.Midas_v5.115.gff3.gz | wc -l
```

**My output:**
```
The file contains 6637 sequence regions.
```
**Does that match with the expectation for this organism?**

* No — the expected chromosome number for **Amphilophus citrinellus is 48**, but the GFF3 lists **6,637 sequence regions**. This is because the assembly is fragmented into many scaffolds/contigs rather than full chromosomes. Such fragmentation is common in genome assemblies and reflects technical limitations, not the actual biological chromosome count.

## 3. How many features does the file contain?

```
zgrep -v "^#" Amphilophus_citrinellus.Midas_v5.115.gff3.gz | wc -l

```
**My Output:**
```
Total number of features: 631,012

```
* Note: 631,012 features include all annotated elements such as genes, mRNAs, exons, CDS, regulatory regions, and other feature types.

## 4. How many genes are listed for this organism?

```
zgrep -v "^#" Amphilophus_citrinellus.Midas_v5.115.gff3.gz | awk '$3=="gene"' | wc -l

```

**My Output:**
```
23696
```

## 5. Is there a feature type that you may have not heard about before? What is the feature and how is it defined? (If there is no such feature, pick a common feature.)

```
zgrep -v "^#" Amphilophus_citrinellus.Midas_v5.115.gff3.gz | cut -f3 | sort | uniq

```

**My output:**

```
biological_region
CDS
exon
five_prime_UTR
gene
J_gene_segment
lnc_RNA
miRNA
mRNA
ncRNA_gene
pseudogene
pseudogenic_transcript
region
rRNA
scRNA
snoRNA
snRNA
three_prime_UTR
transcript
V_gene_segment
```

I have also checked how many times each feature appears

```
zgrep -v "^#" Amphilophus_citrinellus.Midas_v5.115.gff3.gz | cut -f3 | sort | uniq -c | sort -nr

```

**My output:**

```
 279905 exon
 279011 CDS
  31744 mRNA
  23696 gene
   8332 biological_region
   6637 region
    782 ncRNA_gene
    324 rRNA
    199 snoRNA
    180 miRNA
     50 pseudogenic_transcript
     50 pseudogene
     48 snRNA
     19 transcript
     15 V_gene_segment
      9 scRNA
      6 J_gene_segment
      3 lnc_RNA
      1 three_prime_UTR
      1 five_prime_UTR
```


* I found the **pseudogenic transcript** feature new and interesting. A pseudogenic transcript is an RNA sequence produced from a pseudogene, which is a non-functional copy of a gene. Unlike regular genes, pseudogenes generally do not code for functional proteins, but their transcripts can still have regulatory roles or provide insight into genome evolution.


## 6. What are the top-ten most annotated feature types (column 3) across the genome?

```
zgrep -v "^#" Amphilophus_citrinellus.Midas_v5.115.gff3.gz | cut -f3 | sort | uniq -c | sort -nr | head -n 10

```

**Top-ten most annotated feature types across the genome are:**

```
 279905 exon
 279011 CDS
  31744 mRNA
  23696 gene
   8332 biological_region
   6637 region
    782 ncRNA_gene
    324 rRNA
    199 snoRNA
    180 miRNA
```

## 7.Having analyzed this GFF file, does it seem like a complete and well-annotated organism?

* No — while the file includes all major features such as gene, mRNA, exon, and CDS, it also contains a very large number of additional sequence regions (6,637 scaffolds) and many non-coding or regulatory features. This suggests the genome is highly annotated but fragmented, which is common for draft assemblies. Overall, it is reasonably well-annotated, but the large number of scaffolds indicates it may not be fully complete at the chromosome level.

## 8.My insights:

* The GFF3 contains a large number of scaffolds (6,637) compared to the expected 48 chromosomes, reflecting a fragmented assembly.
* There is a rich variety of feature types, including coding (gene, mRNA, CDS, exon) and non-coding or regulatory features (pseudogenic_transcript, miRNA, UTRs).
* The annotation seems comprehensive, allowing analysis of gene structure, transcript diversity, and potential regulatory elements.
* Features like pseudogenic transcripts were new and interesting, highlighting genome complexity and evolution.


* All commands used to generate this analysis are provided in the file:

```
GFF_file_analysis.sh
```






