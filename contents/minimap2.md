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

## 2.2 Mapping DNA reads - Running mapping job with minimap2 against the same index

You must use the parameter -Y, use soft clipping for supplementary alignments, or it will not be compatible with PBSV.

You also need a -R parameter, and add read group information, or else it will not be compatible with PBSV.

```
minimap2 -ax map-hifi --MD -t 16 -Y -R '@RG\tID:SRR11951494\tPL:pacbio\tLB:library\tSM:SRR11951494' /home/qgn1237/qgn1237/1_my_database/GRCh38_p13/minimap2_index/GRCh38.p13.genome.mmi /projects/b1171/qgn1237/2_raw_data/smooth_seq_95_sc_K562_SMRT/SRR11951494/SRR11951494.fastq | samtools sort -@ 16 -m 2G -O BAM -o PC3.bam && samtools index PC3.bam PC3.bai
```

```
minimap2 -ax map-ont --MD -t 8 -Y -R '@RG\tID:XXX\tPL:ont\tLB:library\tSM:XXX' /home/qgn1237/qgn1237/1_my_database/GRCh38_p13/minimap2_index/GRCh38.p13.genome.mmi ../read.fastq | samtools sort -@ 8 -m 2G -O BAM -o mapped.bam && samtools index mapped.bam mapped.bai
```

-a means outputing in SAM format; -x means a preset depends on your task aim.


## 2.3 Mapping RNA reads - It's different from DNA, you need allowing junction.

```
minimap2 -ax splice -k14 --MD -t 8 -Y -R '@RG\tID:PC310cells\tPL:ont\tLB:library\tSM:PC310cells' /home/qgn1237/qgn1237/1_my_database/GRCh38_p13/minimap2_index/GRCh38.p13.genome.mmi ../dorado.fastq | samtools sort -@ 8 -m 2G -O BAM -o VCaP_PCLC_input.bam && samtools index VCaP_PCLC_input.bam  VCaP_PCLC_input.bam.bai
```

```bash
# PacBio data
./MINIMAP2_steps_generator.py \
    --fastq input.fastq \
    --reference ref.fa \
    --sample-name HG002 \
    --data-type pacbio

# ONT DNA data
./MINIMAP2_steps_generator.py \
    --fastq input.fastq \
    --reference ref.fa \
    --sample-name HG002 \
    --data-type ont

# RNA splicing data
./MINIMAP2_steps_generator.py \
    --fastq input.fastq \
    --reference ref.fa \
    --sample-name PC310cells \
    --data-type rna
```