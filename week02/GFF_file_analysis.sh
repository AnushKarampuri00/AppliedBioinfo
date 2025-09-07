wget ftp://ftp.ensembl.org/pub/current_gff3/Amphilophus_citrinellus/Amphilophus_citrinellus.Midas_v5.115.gff3.gz
zgrep "##sequence-region" Amphilophus_citrinellus.Midas_v5.115.gff3.gz | wc -l
zgrep -v "^#" Amphilophus_citrinellus.Midas_v5.115.gff3.gz | wc -l
zgrep -v "^#" Amphilophus_citrinellus.Midas_v5.115.gff3.gz | awk '$3=="gene"' | wc -l
zgrep -v "^#" Amphilophus_citrinellus.Midas_v5.115.gff3.gz | cut -f3 | sort | uniq
zgrep -v "^#" Amphilophus_citrinellus.Midas_v5.115.gff3.gz | cut -f3 | sort | uniq -c | sort -nr
zgrep -v "^#" Amphilophus_citrinellus.Midas_v5.115.gff3.gz | cut -f3 | sort | uniq -c | sort -nr | head -n 10
