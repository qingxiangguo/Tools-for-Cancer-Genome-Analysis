# The installation and usage of PBSV

## 1. About

pbsv calls insertions, deletions, inversions, duplications, and translocations. Both single-sample calling and joint (multi-sample) calling are provided. pbsv is most effective for:

insertions 20 bp to 10 kb  
deletions 20 bp to 100 kb  
inversions 200 bp to 10 kb  
duplications 20 bp to 10 kb  
translocations between different chromosomes or further than 100kb apart on a single chromosome  

## 2. Installation and Usage

### mamba

```bash
mamba install -c bioconda PBSV=2.8.0
```

### 2.1 Discover signatures

It is highly recommended to provide one tandem repeat annotation .bed file of your reference to pbsv discover via --tandem-repeats. This increases sensitivity and recall.

```bash
pbsv discover --tandem-repeats ~/qgn1237/1_my_database/GRCh38_p13/tandem_repeats/human_GRCh38_no_alt_analysis_set.trf.bed ./SRR11951439_sort.bam ./SRR11951439_sort.svsig.gz
```

### 2.2 Call structural variants and assign genotypes

```bash
pbsv call --ccs -j 6 ~/qgn1237/1_my_database/GRCh38_p13/GRCh38.p13.genome.fa /home/qgn1237/qgn1237/4_single_cell_SV_chimera/1_smooth_seq_95_sc_K562_SMRT/SRR11951439/SRR11951439_sort.svsig.gz /home/qgn1237/qgn1237/4_single_cell_SV_chimera/1_smooth_seq_95_sc_K562_SMRT/SRR11951439/SRR11951439_sort.var.vcf
```

-j means tthreads

### 2.3 Output files

SRR11951439_sort.var.vcf
