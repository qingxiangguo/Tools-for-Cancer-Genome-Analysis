# The installation and usage of cuteSV

# 1. About

CuteSV is a reference alignment based SV caller. 

# 2. Installation and Usage

## Installation

```bash
mamba install -c conda-forge biopython # Install BioPython first
mamba install -c bioconda cutesv
```

## Read mapping

```bash
minimap2 -ax map-hifi --MD -t 16 -Y -R '@RG\tID:SRR11951494\tPL:pacbio\tLB:library\tSM:SRR11951494' /home/qgn1237/qgn1237/1_my_database/GRCh38_p13/minimap2_index/GRCh38.p13.genome.mmi /projects/b1171/qgn1237/2_raw_data/smooth_seq_95_sc_K562_SMRT/SRR11951494/SRR11951494.fastq | samtools sort -@ 16 -m 2G -O BAM -o PC3.bam && samtools index PC3.bam PC3.bai
```

## SV calling

The -s parameter should changed based on Depth

For PacBio CCS(HIFI) data:

```bash
cuteSV --threads 16 --genotype --report_readid -l 50 -L 5000000 -r 1000 -q 20 -s 3 --max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5 sorted.bam ~/qgn1237/1_my_database/GRCh38_p13/GRCh38.p13.genome.fa cutesv.vcf ./
```

For ONT data:

```bash
cuteSV --threads 16 --genotype --report_readid -l 50 -L 5000000 -r 1000 -q 20 -s 3 --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3 ../minimap2/PC3_ONT_bulk_all.bam  ~/qgn1237/1_my_database/GRCh38_p13/GRCh38.p13.genome.fa cutesv.vcf ./
```

* min_sv_length is a minimum length threshold a call must have, which is set to 30 nucleotides by default.  
* max_sv_length is a maximum length threshold above which we will not consider calls, set to 10000 by default.  
* min_read_length is a minimum threshold for read lengths. This is used to reject the shortest mapping sequences and is by default defined as 1000 nucleotides.  
* read_support is a minimum threshold for the number of required supporting reads. The default value requires that three of more reads cover an SV.  
* min_read_mapping_quality was defined in the configuration form at the top of the analysis section (a value of 20 by default) and is a minimum threshold used to reject reads with lower quality mapping scores.  

```bash
# short cut
# For ONT data
./CuteSV_steps_generator.py \
    --bam input.bam \
    --reference ref.fa \
    --data-type ont

# For PacBio data
./CuteSV_steps_generator.py \
    --bam input.bam \
    --reference ref.fa \
    --data-type pacbio
```