# Dorado Basecalling Guide

This guide offers a comprehensive overview of how to set up and run Dorado for ONT sequencing data, specifically focusing on the basecalling step.

## Table of Contents
1. About Dorado
2. Installation
3. Usage
   - Simplex Basecalling
   - Duplex Basecalling
4. Channel-wise Data Splitting with POD5
5. Running Dorado with SLURM
6. (Optional) Merging BAM Files

---

## 1. About Dorado

(Detailed introduction about Dorado, its advantages, applications, etc.)

---

## 2. Installation

```bash
 wget https://cdn.oxfordnanoportal.com/software/analysis/dorado-0.3.1-linux-x64.tar.gz
 tar zxf dorado-0.3.1-linux-x64.tar.gz
```

**Note:** For optimal performance, Dorado requires POD5 file input. Convert your `.fast5` files before basecalling. [POD5 File Format](https://github.com/nanoporetech/pod5-file-format)

---

## 3. Usage

### Simplex Basecalling

For running simplex basecalling in Dorado:

```bash
#!/bin/bash -l
#SBATCH --account=b1171
#SBATCH --partition=b1171
#SBATCH --gres=gpu:a100:2
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --mem=100G
#SBATCH --time=70:00:00

cd $SLURM_SUBMIT_DIR

/home/qgn1237/2_software/dorado-0.3.1-linux-x64/bin/dorado basecaller /home/qgn1237/2_software/dorado-0.3.4-linux-x64/bin/dna_r10.4.1_e8.2_400bps_sup@v4.2.0 sample.pod5 > simplex_dorado.bam
```

After obtaining the BAM file, extract fastq reads:

```bash
samtools view duplex_also_simplex_dorado.bam -d dx:0 | samtools fastq > dorado.simplex.fastq
```

### Duplex Basecalling

For running duplex basecalling:

```bash
#!/bin/bash -l
#SBATCH --account=b1171
#SBATCH --partition=b1171
#SBATCH --gres=gpu:a100:2
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --mem=100G
#SBATCH --time=70:00:00

cd $SLURM_SUBMIT_DIR

/home/qgn1237/2_software/dorado-0.3.1-linux-x64/bin/dorado duplex /home/qgn1237/2_software/dorado-0.3.4-linux-x64/bin/dna_r10.4.1_e8.2_400bps_sup@v4.2.0 sample.pod5 > duplex_also_simplex_dorado.bam -t 48
```

After this step, you'll have a BAM file. Separate it into simplex and duplex fastq reads:

```bash
samtools view duplex_also_simplex_dorado.bam -d dx:0 | samtools fastq > dorado.simplex.fastq

samtools view duplex_also_simplex_dorado.bam -d dx:1 | samtools fastq > dorado.duplex.fastq
```

```bash
# Or combine all to fastq
samtools view duplex_also_simplex_dorado.bam | samtools fastq > dorado.simplex.fastq
```

---

## 4. Channel-wise Data Splitting with POD5

Firstly, generate the `summary.tsv`:

```bash
pod5 view <path_to_your_pod5_file> --include "read_id, channel" --output summary.tsv
```

Then, split the data based on channels:

```bash
pod5 subset <path_to_your_pod5_file> --summary summary.tsv --columns channel --output <output_directory_path>
```

---

## 5. Running Dorado with SLURM

Create a new SLURM script:

```bash
#!/bin/bash -l
#SBATCH --account=b1171
#SBATCH --partition=b1171
#SBATCH --gres=gpu:a100:2
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --mem=100G
#SBATCH --time=70:00:00

cd $SLURM_SUBMIT_DIR

for file in <path_to_split_pod5_files>/*.pod5; do
    /home/qgn1237/2_software/dorado-0.3.4-linux-x64/bin/dorado duplex <model_file_path> $file > "${file%.pod5}_dorado.bam"
done
```

Save the script and submit:

```bash
sbatch <name_of_script.sh>
```

---

## 6. (Optional) Merging BAM Files

Merge all BAM files using `samtools`:

```bash
samtools merge merged.bam <path_to_individual_bam_files>/*_dorado.bam
```

---

**End of Guide**

Ensure all software dependencies are installed and paths are correctly specified before running the commands. Always validate the results with the appropriate tools or visualization methods.