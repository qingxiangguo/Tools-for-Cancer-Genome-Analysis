# The installation and usage of SVDSS

## 1. About
A SV caller for long-read single-molecular sequencing data.

## 2. Usage

```bash
# Create the genome index first
SVDSS index --fasta ./GRCh38.p13.genome.fa --index ./GRCh38.p13.genome.fmd --threads 6

SVDSS smooth -t 6 --workdir ./ --bam ../minimap2/PC3_ONT_bulk_all.bam --reference  /projects/b1171/qgn1237/1_my_database/GRCh38_p13/GRCh38.p13.genome.fa

samtools index smoothed.selective.bam smoothed.selective.bai

SVDSS search --index ~/qgn1237/1_my_database/GRCh38_p13/GRCh38.p13.genome.fmd -t 6 --bam ./smoothed.selective.bam

SVDSS assemble --workdir ./ --batches 3 -t 6

SVDSS call -t 6 --workdir ./ --bam ./smoothed.selective.bam --reference ~/qgn1237/1_my_database/GRCh38_p13/GRCh38.p13.genome.fa
```
