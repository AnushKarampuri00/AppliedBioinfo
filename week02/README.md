# GFF File Analysis using UNIX commands

## 1. Tell us about the organism
* Amphilophus citrinellus is a large cichlid fish endemic to the San Juan River and adjacent watersheds in Costa Rica and Nicaragua. In the aquarium trade _A. citrinellus_ is often sold under the trade name of **Midas cichlid**. The Midas cichlid (_Amphilophus citrinellus_) has 48 chromosomes. This means its diploid number is 2n = 48.

## 2. How many sequence regions (chromosomes) does the file contain? Does that match with the expectation for this organism?
* In GFF3, sequence regions are in the header lines starting with ##sequence-region.

"""

zgrep "##sequence-region" Amphilophus_citrinellus.Midas_v5.115.gff3.gz | wc -l

"""




