# The installation and usage of SVABA

## 1. About

SvABA is a method for detecting structural variants in sequencing data using genome-wide local assembly

## 2. Usage

### 2.1 Assume that you've already index and sort the BAM file first, you also have a BWM indexed genome

```bash
svaba run -t in.bam -p 8 -a your_tag -G /projects/b1171/qgn1237/1_my_database/GRCh38_p13/GRCh38.p13.genome.fa
```

### 2.2 Output files

SRR11951439_sort.var.vcf
