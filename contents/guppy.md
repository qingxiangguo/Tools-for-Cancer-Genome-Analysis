# The installation and usage of Guppy

## 1. About

Read fast5/pod5 and give you fastq.

## 2. Installation

wget https://cdn.oxfordnanoportal.com/software/analysis/ont-guppy_6.5.7_linux64.tar.gz


## 3. Usage 

### First round of simplex basecalling

Select config file according to your flow cell, kit and modified base requirement in data folder:

sup: Super-accurate basecalling. hac: High accuracy basecalling. These are the configurations that will be selected when a kit and flow cell are specified on the command-line instead of a specific config file. fast: Fast basecalling. sketch: Sketch basecalling. This is primarily for use with adaptive sampling on the MinION Mk1C device to minimise latency.

For sqk LSK114/R10.4.1 and 5mc and 5hmc (can't do 6ma by Guppy now)

dna_r10.4.1_e8.2_400bps_5khz_modbases_5hmc_5mc_cg_hac_mk1c.cfg


```bash
fish_add_path /home/qgn1237/2_software/ont-guppy/bin
```

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

/home/qgn1237/2_software/ont-guppy/bin/guppy_basecaller --disable_pings --input_path /projects/b1171/qgn1237/2_raw_data/20230616_PC3_bulk_genome_ONT/pc3_dna_bulk_1/pc3_bulk/20230612_1911_MC-114785_FAW84522_f68d5359/pod5  --recursive  --gpu_runners_per_device 50 --num_callers 48 --save_path /home/qgn1237/qgn1237/2_raw_data/20230616_PC3_bulk_genome_ONT/pc3_dna_bulk_1/pc3_bulk/20230612_1911_MC-114785_FAW84522_f68d5359/20230627_2nd_duplex_basecalling/1st_simplex_basecalling --config /home/qgn1237/2_software/ont-guppy/data/dna_r10.4.1_e8.2_400bps_5khz_modbases_5hmc_5mc_cg_fast_mk1c.cfg --device "cuda:0 cuda:1" --do_read_splitting
```

```
- `--gpu_runners_per_device`: Sets the number of tasks run in parallel per GPU.
- `--num_callers`: Sets the number of parallel basecalling tasks.
- `--ntasks`: Sets the total number of CPU cores used for the job.
```

In your case, with `--gpu_runners_per_device=8` and `--num_callers=16` across 2 GPUs, you have balanced the load. Each GPU supports 8 tasks, and you have total 16 tasks running in parallel, (less than) aligning with your CPU cores specified by `--ntasks=24`.

The methylation information in the FASTQ file produced by Guppy, if a suitable configuration file is used (like dna_r9.4.1_450bps_modbases_dam-dcm-cpg_hac.cfg), can be utilized by certain downstream analysis tools specialized in DNA methylation. These tools, including Megalodon and DeepSignal, can read the methylation probability scores from the Guppy output FASTQ file and perform further analysis, such as pinpointing the precise location and type of methylation.

Utilizing the methylation scores from the Guppy-produced FASTQ can significantly enhance the speed and efficiency of these methylation analysis tools, as the scores are precomputed during the basecalling process by Guppy, preventing the need for duplicate computations in downstream analyses.

However, note that not all methylation analysis tools can directly use the methylation information from the Guppy-produced FASTQ file. Certain tools like Nanopolish depend on the raw ionic current data stored in the FAST5 files and use their own algorithms to detect methylation. For such tools, even if you used a methylation-related configuration file while running Guppy, the methylation information from Guppy cannot be directly utilized and the methylation detection needs to be recomputed.

If you find error like "guppy_basecaller with gpu: error while loading shared libraries: libcuda.so.1", no worry, because you don't have cuda in your loca, it's on a GPU node.

If running, it will show something like,

```bash
ONT Guppy basecalling software version 6.5.7+ca6d6af, minimap2 version 2.24-r1122
config file:        /home/qgn1237/2_software/ont-guppy/data/dna_r10.4.1_e8.2_400bps_5khz_modbases_5hmc_5mc_cg_hac_mk1c.cfg
model file:         /home/qgn1237/2_software/ont-guppy/data/template_r10.4.1_e8.2_400bps_5khz_hac.jsn
input path:         /projects/b1171/qgn1237/2_raw_data/20230616_PC3_bulk_genome_ONT/pc3_dna_bulk_1/pc3_bulk/20230612_1911_MC-114785_FAW84522_f68d5359/pod5
save path:          /projects/b1171/qgn1237/2_raw_data/20230616_PC3_bulk_genome_ONT/pc3_dna_bulk_1/pc3_bulk/20230612_1911_MC-114785_FAW84522_f68d5359/basecalled_fastq
chunk size:         2000
chunks per runner:  48
minimum qscore:     9
records per file:   4000
num basecallers:    32
gpu device:         cuda:0 cuda:1
kernel path:
runners per device: 50

Use of this software is permitted solely under the terms of the end user license agreement (EULA).
By running, copying or accessing this software, you are demonstrating your acceptance of the EULA.
The EULA may be found in /home/qgn1237/2_software/ont-guppy/bin
Found 353 input read files to process.
Init time: 4437 ms

0%   10   20   30   40   50   60   70   80   90   100%
|----|----|----|----|----|----|----|----|----|----|
```

Merged all FASTQ files inside the "pass" folder of Guppy results with "cat" command and obtained single FASTQ.

My nanopore DNA pipeline is Guppy + PycoQC (analysis the sequencing_summary.txt) + Porechop (remove adaptors) + Nanoflit trimming (get clean reads) + mapping with minimap2 + Qualimap BAM file QC.


### Duplex basecalling with the help of duplex_tools

You need to do simplex calling and enable read splitting first.

Please refer to Duplex tools to get the pair_ids_filtered.txt

Then

Start duplex basecalling.

The input is the original .fast5/.pod5 file and the pair_ids_filtered.txt file generated above, as follows

```bash
#!/bin/bash -l
#SBATCH --account=b1171
#SBATCH --partition=b1171
#SBATCH --gres=gpu:a100:2
# set the number of nodes
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --mem=100G
#SBATCH --time=80:00:00

# run the application
cd $SLURM_SUBMIT_DIR

/home/qgn1237/2_software/ont-guppy/bin/guppy_basecaller_duplex --disable_pings --input_path /projects/b1171/qgn1237/2_raw_data/20230616_PC3_bulk_genome_ONT/pc3_dna_bulk_1/pc3_bulk/20230612_1911_MC-114785_FAW84522_f68d5359/pod5 --recursive --gpu_runners_per_device 50 --num_callers 48 --save_path /home/qgn1237/qgn1237/2_raw_data/20230616_PC3_bulk_genome_ONT/pc3_dna_bulk_1/pc3_bulk/20230612_1911_MC-114785_FAW84522_f68d5359/20230627_2nd_duplex_basecalling/2nd_duplex_basecalling/ --config /home/qgn1237/2_software/ont-guppy/data/dna_r10.4.1_e8.2_400bps_5khz_modbases_5hmc_5mc_cg_sup.cfg --device "cuda:0 cuda:1" --duplex_pairing_mode from_pair_list --duplex_pairing_file /projects/b1171/qgn1237/2_raw_data/20230616_PC3_bulk_genome_ONT/pc3_dna_bulk_1/pc3_bulk/20230612_1911_MC-114785_FAW84522_f68d5359/20230627_2nd_duplex_basecalling/duplex_output/pair_ids_filtered.txt --chunks_per_runner 1200
```

Note that you must specify the --chunks_per_runner 1200, or else there will be a length problem.

This will generate duplex fastq files.

### Combine the results of duplex and simplex fastqs 

### Overview

In nanopore sequencing, duplex sequencing offers an additional layer of accuracy in base calling, owing to the independent sequencing of both DNA strands. However, the data obtained from duplex sequencing is not readily combined with the simplex sequencing data. This guide will walk you through the steps to combine the FASTQ files from both types of sequencing, removing any redundant reads in the process.

### Prerequisites

You'll need:

- Your duplex and simplex FASTQ files
- [Seqkit](https://github.com/shenwei356/seqkit) installed

### Procedure

1. **Concatenate all Simplex FASTQ files**: Use the `cat` command to combine all the FASTQ files from the simplex sequencing. 

    ```bash
    cat simplex-guppy/*.fastq.gz > simplex.fastq.gz 
    ```

2. **Extract all Simplex run names**: The `seqkit` tool can be used to extract all read IDs from the combined FASTQ file.

    ```bash
    seqkit seq --name simplex.fastq.gz > simplex_ids_txt 
    ```

3. **Substitute Simplex reads with Duplex reads and combine**: Use the `seqkit grep` command to extract all the reads from the simplex file that are not present in the duplex file, and combine them with all reads from the duplex file.

    ```bash
    { sed 's/ /\n/' pair_ids_filtered.txt | \
    seqkit grep -v -f - simplex.fastq.gz ; \
    zcat duplex-guppy/../fastq_pass/*.fastq.gz ; } \
    | gzip - > combined.fastq.gz 
    ```

4. **Output**: You'll end up with a combined FASTQ file (`combined.fastq.gz`) that includes both simplex and duplex reads, with any redundant reads removed.
