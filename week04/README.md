# Retrieving the genomic Data
I selected the Zika virus-related study (group 2), I have selected the RefSeq assembly GCF_000882815.3 from Zika virus to retrieve the genome and annotation data for analysis.

## Download the genome FASTA and GFF files

``` bash
# Create a directory for the genome data
mkdir -p ~/genome_data
cd ~/genome_data

# Download the genome FASTA file using curl
curl -O https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/882/815/GCF_000882815.3_ViralProj36615/GCF_000882815.3_ViralProj36615_cds_from_genomic.fna.gz

# Download the GFF annotation file using curl
curl -O https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/882/815/GCF_000882815.3_ViralProj36615/GCF_000882815.3_ViralProj36615_genomic.gff.gz


# Unzip the downloaded files
gunzip *.gz

```

## Length of the genome


```
grep -v ">" GCF_000882815.3_ViralProj36615_cds_from_genomic.fna | wc -c

```
Output:

```
The genome is 10389 basepairs long
```

## How many features of each type exist in the GFF file

```
grep -v "#" GCF_000882815.3_ViralProj36615_genomic.gff | cut -f3 | sort | uniq -c | sort -nr
```

Output

```
14 mature_protein_region_of_CDS
   1 three_prime_UTR
   1 region
   1 gene
   1 five_prime_UTR
   1 CDS

```

## Identifying the longest gene

```
awk '$3 == "gene" {print $9, $5 - $4}' GCF_000882815.3_ViralProj36615_genomic.gff | sort -k2,2nr | head -1

```

Output:

```
ID=gene-ZIKV_gp1;Dbxref=GeneID:7751225;Name=POLY;gbkey=Gene;gene=POLY;gene_biotype=protein_coding;locus_tag=ZIKV_gp1 10259

```

**Explaination:**

Gene name : POLY 

NCBI Gene ID: 7751225

**Function:** POLY is a poly protein consisting three structural proteins (capsid, pre-membrane and envelope)


* The POLY gene of Zika virus encodes a single long polyprotein that is later cleaved into individual viral proteins. These include the structural proteins (C, prM, and E), which form the viral particle, and the non-structural proteins (NS1–NS5), which are essential for viral replication and immune evasion. This gene represents the entire coding capacity of the Zika virus genome.

## Identifying another gene:

```
awk -F"\t" '$3=="gene" {print $5-$4 "\t" $9}' GCF_000882815.3_ViralProj36615_genomic.gff | sort -n | head -1

```

Output:
```
10259	ID=gene-ZIKV_gp1;Dbxref=GeneID:7751225;Name=POLY;gbkey=Gene;gene=POLY;gene_biotype=protein_coding;locus_tag=ZIKV_gp1
```
* The Zika virus genome contains a single polyprotein gene (POLY), which is cleaved into structural and non-structural proteins. Since this is the only feature labeled as gene in the GFF file, no additional genes can be identified; other proteins are annotated as CDS features rather than separate genes.


## Analysing the intragenomic space and packing nature of the genome


**1. Finding the number of genes:**

```
awk '$3=="gene"' GCF_000882815.3_ViralProj36615_genomic.gff | wc -l
```

Output:

```
1
```

**2. Measure the genomic size - total number of basepairs**

```
grep -v ">" GCF_000882815.3_ViralProj36615_cds_from_genomic.fna | tr -d '\n' | wc -c

```

Output:

```
10260
```

**3. total length of all genes - how many basepairs does all the genes cover**
```
awk '$3=="gene" {sum += $5-$4} END {print sum}' GCF_000882815.3_ViralProj36615_genomic.gff

```

Output:

```
10259
```

**4. To find the largest and average intragenic spaces among the genes**

```
# for largest intragenic spaces
awk '$3=="gene" {print $1,$4,$5}' GCF_000882815.3_ViralProj36615_genomic.gff | \
  sort -k1,1 -k2,2n | \
  awk '{if(prev_chr==$1){dist=$2-prev_end; if(dist>0) print dist} prev_chr=$1; prev_end=$3}' | \
  sort -nr | head

# for average intragenic spaces
awk '$3=="gene" {print $1,$4,$5}' GCF_000882815.3_ViralProj36615_genomic.gff | \
  sort -k1,1 -k2,2n | \
  awk '{if(prev_chr==$1){dist=$2-prev_end; if(dist>0){print dist; sum+=dist; count++}} prev_chr=$1; prev_end=$3} END {print "Average intergenic distance:", sum/count}'

```

**Interpretation:** 

* The Zika virus genome contains only one gene feature (POLY). Therefore, there are no consecutive genes to calculate intergenic distances, and no output is produced by the distance calculation command.


## An estimation on how much genome is covered through IGV

Simplify the actual GFF file

```
awk -F'\t' '$3=="CDS"' GCF_000882815.3_ViralProj36615_genomic.gff > CDS_only.gff
```

**Visualise it in IGV along with the genome file**

![CDS coverage on chromosome 11](images/intragenic_spaces.png)
Figure 1: Visualizing the genome of Zika virus in IGV.

* **Note:** The GFF file for Zika virus contains only a single polyprotein gene (POLY) and its corresponding CDS spanning almost the entire genome. As a result, there are no additional gene or CDS features to display. This is why the annotation track appears empty in IGV at a typical zoom level.


## Alternative genome builds and research questions

Looking beyond GRCh38, other genome references can open up new perspectives.
* The Zika virus genome **(NC_012532.1)** could be used to study viral proteins and how they interact with host responses.
* The older human build, **GRCh37**, might reveal differences in annotations or alternative isoforms linked to cell cycle control.
* Meanwhile, a primate genome such as rhesus macaque **(GCF_003339765.1)** could help compare infection patterns across species and point to host-specific resistance factors.



