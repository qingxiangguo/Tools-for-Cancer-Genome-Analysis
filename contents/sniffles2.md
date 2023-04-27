# The installation and usage of cuteSV

# 1. About

Sniffles2 accurately detect SVs on germline, somatic and population-level for PacBio and Oxford Nanopore read data.

# 2. Installation and Usage

## Installation

```bash
mamba install -c bioconda sniffles=2.0
```

## Read mapping

```
minimap2 -ax map-hifi --MD -t 16 -Y -R '@RG\tID:SRR11951494\tPL:pacbio\tLB:library\tSM:SRR11951494' /home/qgn1237/qgn1237/1_my_database/GRCh38_p13/minimap2_index/GRCh38.p13.genome.mmi /projects/b1171/qgn1237/2_raw_data/smooth_seq_95_sc_K562_SMRT/SRR11951494/SRR11951494.fastq | samtools sort -@ 16 -m 2G -O BAM -o PC3.bam && samtools index PC3.bam PC3.bai
```

## SV calling

Specify tandem repeat annotations (for improved call accuracy), reference (for DEL sequences) and non-germline mode for detecting rare SVs

```bash
sniffles --input /projects/b1171/qgn1237/5_impact_of_depth_sv_detection/2_cell_line_SV/SKBR3_CLR/5X_depth/SKBR3_CLR_5X.bam --vcf out.vcf --tandem-repeats /projects/b1171/qgn1237/1_my_database/GRCh38_p13/sniffles2_compatible_tandem_repeat/human_GRCh38_no_alt_analysis_set.trf.bed --reference /projects/b1171/qgn1237/1_my_database/GRCh38_p13/GRCh38.p13.genome.fa --minsupport auto --minsvlen 50
```