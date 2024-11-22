# The installation and usage of Lumpy
# 1. About
A probabilistic framework for short-read structural variant discovery
# 2. Installation and Usage

## mamba
```
mamba create -n lumpy_env python=2.7
mamba install -c bioconda lumpy-sv
```

## Run

```bash
#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics ## Required: (buyin, short, normal, long, gengpu, genhimem, etc)
#SBATCH --time=48:00:00 ## Required: How long will the job need to run (remember different partitions have restrictions on this parameter)
#SBATCH --nodes=1 ## how many computers/nodes do you need (no default)
#SBATCH --ntasks-per-node=1 ## how many cpus or processors do you need on per computer/node (default value 1)
#SBATCH --mem=40G ## how much RAM do you need per computer/node (this affects your FairShare score so be careful to not ask for more than you need))
#SBATCH --job-name=allen_sam1 ## When you run squeue -u  this is how you can identify the job
#SBATCH --output=output.log ## standard out and standard error goes to this file

# A regular comment in Bash
/home/qgn1237/2_software/mambaforge/bin/mamba init
source ~/.bashrc
mamba activate lumpy_env
lumpyexpress -B /projects/b1171/qgn1237/5_impact_of_depth_sv_detection/1_reads_mapping/bwa-mem2/SKBR3_NGS.bam -o SKBR3_NGS.vcf
```

```bash
python LUMPY_steps_generator.py --bam input.bam

python LUMPY_steps_generator.py --bam input.bam --output custom_output.vcf
```