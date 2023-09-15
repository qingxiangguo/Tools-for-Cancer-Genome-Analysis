# Bedtools
Bedtools allows one to intersect, merge, count, complement, and shuffle genomic intervals from multiple files in widely-used genomic file formats such as BAM, BED, GFF/GTF, VCF. While each individual tool is designed to do a relatively simple task (e.g., intersect two interval files), quite sophisticated analyses can be conducted by combining multiple bedtools operations on the UNIX command line.

```
wget https://github.com/arq5x/bedtools2/releases/download/v2.29.1/bedtools-2.29.1.tar.gz
tar -zxvf bedtools-2.29.1.tar.gz
cd bedtools2
make
```

### Sort bed file

```bash
sortBed -g ~/qgn1237/1_my_database/GRCh38_p13/GRCh38.p13.genome.fa.fai -i HACk.random.bed  > test
# The chromosomes have been sorted alphabetically by default. You have to specify -g, honestly this is a stupid design.
```
