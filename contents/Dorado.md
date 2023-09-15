# The installation and usage of Guppy

## 1. About

## 2. Installation

```bash
 wget https://cdn.oxfordnanoportal.com/software/analysis/dorado-0.3.1-linux-x64.tar.gz
 tar zxf dorado-0.3.1-linux-x64.tar.gz
```

For optimal performance, Dorado requires POD5 file input. Please convert your .fast5 files before basecalling.
https://github.com/nanoporetech/pod5-file-format

## 3. Usage

## Run simplex basecalling in Dorado

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

/home/qgn1237/2_software/dorado-0.3.1-linux-x64/bin/dorado basecaller sample.pod5 > simplex_dorado.bam
```

Then you'll get a single bam file, now try to get fastq reads from it.

```bash
samtools view duplex_also_simplex_dorado.bam -d dx:0 | samtools fastq > dorado.simplex.fastq
```

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

/home/qgn1237/2_software/dorado-0.3.1-linux-x64/bin/dorado duplex sample.pod5 > duplex_also_simplex_dorado.bam
```

Then you'll get a single bam file, now try to separate them into simplex fastq reads and duplex fastq reads.

```bash
samtools view duplex_also_simplex_dorado.bam -d dx:0 | samtools fastq > dorado.simplex.fastq

samtools view duplex_also_simplex_dorado.bam -d dx:1 | samtools fastq > dorado.duplex.fastq
```