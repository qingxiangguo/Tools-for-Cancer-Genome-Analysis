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

### 2.1 Reads mapping

```bash
minimap2 -ax map-hifi --MD -t 16 -Y -R '@RG\tID:SRR11951494\tPL:pacbio\tLB:library\tSM:SRR11951494' /home/qgn1237/qgn1237/1_my_database/GRCh38_p13/minimap2_index/GRCh38.p13.genome.mmi /projects/b1171/qgn1237/2_raw_data/smooth_seq_95_sc_K562_SMRT/SRR11951494/SRR11951494.fastq | samtools sort -@ 16 -m 2G -O BAM -o PC3.bam && samtools index PC3.bam PC3.bai
```

### 2.2 Discover signatures

It is highly recommended to provide one tandem repeat annotation .bed file of your reference to pbsv discover via --tandem-repeats. This increases sensitivity and recall.

```bash
pbsv discover --tandem-repeats ~/qgn1237/1_my_database/GRCh38_p13/tandem_repeats/human_GRCh38_no_alt_analysis_set.trf.bed ./SRR11951439_sort.bam ./SRR11951439_sort.svsig.gz
```

### 2.3 Call structural variants and assign genotypes

For pacbio

```bash
pbsv call --ccs -j 6 ~/qgn1237/1_my_database/GRCh38_p13/GRCh38.p13.genome.fa /home/qgn1237/qgn1237/4_single_cell_SV_chimera/1_smooth_seq_95_sc_K562_SMRT/SRR11951439/SRR11951439_sort.svsig.gz /home/qgn1237/qgn1237/4_single_cell_SV_chimera/1_smooth_seq_95_sc_K562_SMRT/SRR11951439/SRR11951439_sort.var.vcf
```

For nanopore

```bash
pbsv call --ccs -j 6 ~/qgn1237/1_my_database/GRCh38_p13/GRCh38.p13.genome.fa /home/qgn1237/qgn1237/4_single_cell_SV_chimera/1_smooth_seq_95_sc_K562_SMRT/SRR11951439/SRR11951439_sort.svsig.gz /home/qgn1237/qgn1237/4_single_cell_SV_chimera/1_smooth_seq_95_sc_K562_SMRT/SRR11951439/SRR11951439_sort.var.vcf
```

-j means threads

### 2.4 Output files

SRR11951439_sort.var.vcf
