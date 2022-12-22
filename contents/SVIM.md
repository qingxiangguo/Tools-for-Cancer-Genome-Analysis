# The installation and usage of SVIM

## 1. About

SVIM is able to detect and classify the following six classes of structural variation: deletions, insertions, inversions, tandem duplications, interspersed duplications and translocations

## 2. Installation and Usage

### mamba

```bash
mamba install -c bioconda svim
```

### 2.1 Assume that you've already index and sort the BAM file first

```bash
svim alignment /home/qgn1237/qgn1237/4_single_cell_SV_chimera/1_smooth_seq_95_sc_K562_SMRT/SRR11951439/svim /projects/b1171/qgn1237/4_single_cell_SV_chimera/1_smooth_seq_95_sc_K562_SMRT/SRR11951439/SRR11951439_sort.bam ~/qgn1237/1_my_database/GRCh38_p13/GRCh38.p13.genome.fa
```

### 2.3 Output files

SRR11951439_sort.var.vcf
