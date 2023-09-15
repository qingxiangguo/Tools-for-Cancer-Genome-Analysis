# The installation and usage of Guppy

## 1. About

## 2. Installation

```bash
 wget https://cdn.oxfordnanoportal.com/software/analysis/dorado-0.3.1-linux-x64.tar.gz
 tar zxf dorado-0.3.1-linux-x64.tar.gz
```

## 3. Usage

## Run duplex basecalling in Dorado

```bash
#!/bin/bash -l
#SBATCH --account=b1171
#SBATCH --partition=b1171
#SBATCH --gres=gpu:a100:2
# set the number of nodes
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --mem=100G
# set max wallclock time
#SBATCH --time=70:00:00

# run the application
cd $SLURM_SUBMIT_DIR

/home/qgn1237/2_software/dorado-0.3.1-linux-x64/bin/dorado duplex /home/qgn1237/2_software/dorado-0.3.1-linux-x64/bin/dna_r10.4.1_e8.2_400bps_sup@v4.2.0 /projects/b1171/qgn1237/2_raw_data/20230616_PC3_bulk_genome_ONT/pc3_dna_bulk_1/pc3_bulk/20230612_1911_MC-114785_FAW84522_f68d5359/pod5 > /home/qgn1237/qgn1237/2_raw_data/20230616_PC3_bulk_genome_ONT/pc3_dna_bulk_1/pc3_bulk/20230612_1911_MC-114785_FAW84522_f68d5359/20230630_dorado_duplex_basecalling/duplex_also_simplex_dorado.bam
```

Then you'll get a single bam file, now try to separate them into simplex fastq reads and duplex fastq reads.

```bash
samtools view duplex_also_simplex_dorado.bam -d dx:0 | samtools fastq > dorado.simplex.fastq

samtools view duplex_also_simplex_dorado.bam -d dx:1 | samtools fastq > dorado.duplex.fastq
```