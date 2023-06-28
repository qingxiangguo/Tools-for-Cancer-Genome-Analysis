# The installation and usage of Duplex tools

# 1. About

for duplex calling of guppy, needs the output of first round of guppy simplex calling.

# 2. Installation and Usage

## Installation


## Usage

```bash
duplex_tools pairs_from_summary /path/to/your/output/basecalled_fastq/sequencing_summary.txt /path/to/your/output/
```

Above step will generate pair_ids.txt.

```bash
duplex_tools filter_pairs /path/to/your/output/pair_ids.txt /path/to/your/simplex_fastq_output/basecalled_fastq
```

Get pair_ids_filtered.txt.

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
#SBATCH --chdir=/home/qgn1237/qgn1237/2_raw_data/20230616_PC3_bulk_genome_ONT/pc3_dna_bulk_1/pc3_bulk/20230612_1911_MC-114785_FAW84522_f68d5359
# set max wallclock time
#SBATCH --time=48:00:00

# run the application
cd $SLURM_SUBMIT_DIR

/home/qgn1237/2_software/ont-guppy/bin/guppy_basecaller_duplex --disable_pings --input_path simplex_pod5  --recursive  --gpu_runners_per_device 50 --num_callers 48 --save_path duplex_fastq --config /home/qgn1237/2_software/ont-guppy/data/dna_r10.4.1_e8.2_400bps_5khz_modbases_5hmc_5mc_cg_sup.cfg --device "cuda:0 cuda:1" --duplex_pairing_mode from_pair_list --duplex_pairing_file /path/to/your/output/pair_ids_filtered.txt