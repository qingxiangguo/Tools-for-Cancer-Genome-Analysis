# fasterq-dump
The fasterq-dump tool extracts data in FASTQ- or FASTA-format from SRA-accessions.  

Download raw data using prefetch + fasterq-dump

```
prefetch SRR9736820 -O /home/qgn1237/qgn1237/2_raw_data/
```

```
#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics ## Required: (buyin, short, normal, long, gengpu, genhimem, etc)
#SBATCH --time=03:00:00 ## Required: How long will the job need to run (remember different partitions have restrictions on this parameter)
#SBATCH --nodes=1 ## how many computers/nodes do you need (no default)
#SBATCH --ntasks-per-node=1 ## how many cpus or processors do you need on per computer/node (default value 1)
#SBATCH --mem=1G ## how much RAM do you need per computer/node (this affects your FairShare score so be careful to not ask for more than you need))
#SBATCH --job-name=allen_download ## When you run squeue -u  this is how you can identify the job
#SBATCH --output=output.log ## standard out and standard error goes to this file

# A regular comment in Bash
fasterq-dump --split-files /home/qgn1237/qgn1237/2_raw_data/SRR9736820 -O /home/qgn1237/qgn1237/2_raw_data/SRR9736820 
```
