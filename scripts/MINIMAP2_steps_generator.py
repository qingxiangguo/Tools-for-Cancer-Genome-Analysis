#!/usr/bin/env python3
import argparse
import os
from pathlib import Path

def create_parser():
    parser = argparse.ArgumentParser(
        description='Generate SLURM script for Minimap2 pipeline')
    
    parser.add_argument('--fastq', required=True, type=str,
                      help='Input FASTQ file path')
    parser.add_argument('--reference', required=True, type=str,
                      help='Reference genome path')
    parser.add_argument('--sample-name', required=True, type=str,
                      help='Sample name (used in read group and output files)')
    parser.add_argument('--data-type', required=True, 
                      choices=['pacbio', 'ont', 'rna'],
                      help='Data type: pacbio, ont, or rna')
    parser.add_argument('--threads', type=int, default=8,
                      help='Number of threads (default: 8)')
    parser.add_argument('--mem', type=str, default='40G',
                      help='Memory allocation (default: 40G)')
    parser.add_argument('--mem-per-thread', type=str, default='2G',
                      help='Memory per thread for samtools sort (default: 2G)')
    parser.add_argument('--time', type=str, default='05:00:00',
                      help='Wall time limit (format: HH:MM:SS, default: 05:00:00)')
    
    return parser

def get_minimap2_params(data_type, sample_name):
    """Return Minimap2 parameters based on data type"""
    params = {
        'pacbio': {
            'mode': 'map-hifi',
            'platform': 'pacbio'
        },
        'ont': {
            'mode': 'map-ont',
            'platform': 'ONT'
        },
        'rna': {
            'mode': 'splice -k14',
            'platform': 'ont'
        }
    }
    mode = params[data_type]['mode']
    platform = params[data_type]['platform']
    read_group = f"'@RG\\tID:{sample_name}\\tPL:{platform}\\tLB:library\\tSM:{sample_name}'"
    return mode, read_group

def generate_slurm_script(args):
    # Get data type specific parameters
    mode, read_group = get_minimap2_params(args.data_type, args.sample_name)
    output_bam = f"{args.sample_name}.bam"
    
    script_content = f'''#!/bin/bash -l
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --nodes=1
#SBATCH --ntasks={args.threads}
#SBATCH --mem={args.mem}
#SBATCH --time={args.time}
#SBATCH --mail-type=ALL

cd $SLURM_SUBMIT_DIR

/home/qgn1237/2_software/mambaforge/bin/mamba init
source ~/.bashrc
mamba activate mamba666

# Run Minimap2 and pipe to samtools sort
minimap2 -ax {mode} \\
    --MD \\
    -t {args.threads} \\
    -Y \\
    -R {read_group} \\
    {os.path.abspath(args.reference)} \\
    {os.path.abspath(args.fastq)} | \\
    samtools sort \\
        -@ {args.threads} \\
        -m {args.mem_per_thread} \\
        -O BAM \\
        -o {output_bam}

if [ $? -eq 0 ]; then
    # Index the BAM file
    samtools index {output_bam} {output_bam}.bai
else
    echo "Minimap2 alignment or sorting failed"
    exit 1
fi

if [ $? -ne 0 ]; then
    echo "BAM indexing failed"
    exit 1
fi
'''
    
    with open('run_minimap2_steps.slurm', 'w') as f:
        f.write(script_content)
    
    print(f"Generated SLURM script: run_minimap2_steps.slurm")

def main():
    parser = create_parser()
    args = parser.parse_args()
    generate_slurm_script(args)

if __name__ == "__main__":
    main()
