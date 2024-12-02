#!/usr/bin/env python3
import argparse
import os
from pathlib import Path

def create_parser():
    parser = argparse.ArgumentParser(
        description='Generate SLURM script for Sniffles2 pipeline')
    
    parser.add_argument('--bam', required=True, type=str,
                      help='Input BAM file path')
    parser.add_argument('--reference', required=True, type=str,
                      help='Reference genome path')
    parser.add_argument('--tandem-repeats', required=True, type=str,
                      help='Tandem repeats BED file path')
    parser.add_argument('--threads', type=int, default=6,
                      help='Number of threads (default: 6)')
    parser.add_argument('--mem', type=str, default='50G',
                      help='Memory allocation (default: 50G)')
    parser.add_argument('--time', type=str, default='06:00:00',
                      help='Wall time limit (format: HH:MM:SS, default: 06:00:00)')
    parser.add_argument('--output', type=str, default='sniffles.vcf',
                      help='Output VCF filename (default: sniffles.vcf)')
    
    return parser

def generate_slurm_script(args):
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
mamba activate sniffles2

# Run Sniffles2
sniffles \\
    --input {os.path.abspath(args.bam)} \\
    --vcf {args.output} \\
    --tandem-repeats {os.path.abspath(args.tandem_repeats)} \\
    --reference {os.path.abspath(args.reference)} \\
    --minsupport auto \\
    --minsvlen 50 \\
    --threads {args.threads}

if [ $? -ne 0 ]; then
    echo "Sniffles2 analysis failed"
    exit 1
fi
'''
    
    with open('run_sniffles2_steps.slurm', 'w') as f:
        f.write(script_content)
    
    print(f"Generated SLURM script: run_sniffles2_steps.slurm")

def main():
    parser = create_parser()
    args = parser.parse_args()
    generate_slurm_script(args)

if __name__ == "__main__":
    main()
