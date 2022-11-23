# The installation and usage of Minimap2
# 1. About
the use-case, dataset and the running machine.
# 2. Installation and Usage
## mamba
```
mamba install -c bioconda minimap2
```

## 2.1 Build genome index, so you don't need to build this everytime
The mmi file is the index
```
minimap2 -d /home/qgn1237/qgn1237/1_my_database/GRCh38_p13/minimap2_index/GRCh38.p13.genome.mmi /projects/b1171/qgn1237/1_my_database/GRCh38_p13/GRCh38.p13.genome.fa -t 6
```

## 2.2 Running mapping job with minimap2 against the same index

You must use the parameter -Y, use soft clipping for supplementary alignments, or it will not be compatible with PBSV.

You also need a -R parameter, and add read group information, or else it will not be compatible with PBSV.

```
minimap2 -ax map-hifi -t 24 -Y -R '@RG\tID:SRR11951494\tPL:pacbio\tLB:library\tSM:SRR11951494' /home/qgn1237/qgn1237/1_my_database/GRCh38_p13/minimap2_index/GRCh38.p13.genome.mmi /projects/b1171/qgn1237/2_raw_data/smooth_seq_95_sc_K562_SMRT/SRR11951494/SRR11951494.fastq > /home/qgn1237/qgn1237/4_single_cell_SV_chimera/1_smooth_seq_95_sc_K562_SMRT/SRR11951494/SRR11951494.sam
```
-a means outputing in SAM format; -x means a preset depends on your task aim.
