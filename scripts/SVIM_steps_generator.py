#!/usr/bin/env python3
import argparse
import os
from pathlib import Path

def create_parser():
    parser = argparse.ArgumentParser(
        description='Generate SLURM script for SVIM pipeline')
    
    parser.add_argument('--bam', required=True, type=str,
                      help='Input BAM file path')
    parser.add_argument('--reference', required=True, type=str,
                      help='Reference genome path')
    parser.add_argument('--threads', type=int, default=1,
                      help='Number of threads (default: 1)')
    parser.add_argument('--mem', type=str, default='10G',
                      help='Memory allocation (default: 10G)')
    parser.add_argument('--time', type=str, default='08:00:00',
                      help='Wall time limit (format: HH:MM:SS, default: 08:00:00)')
    parser.add_argument('--quality', type=int, default=10,
                      help='Quality threshold for VCF filtering (default: 10)')
    
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
mamba activate mamba666

# Run SVIM
svim alignment ./ {os.path.abspath(args.bam)} {os.path.abspath(args.reference)}

# Filter VCF based on quality (wait for SVIM to complete)
if [ $? -eq 0 ]; then
    filter_vcf_based_on_quality.py variants.vcf {args.quality} > filtered_variant_q{args.quality}.vcf
else
    echo "SVIM analysis failed"
    exit 1
fi

# Remove temporary files and directories (wait for filtering to complete)
if [ $? -eq 0 ]; then
    rm -rf candidates signatures sv-* SVIM*
else
    echo "VCF filtering failed"
    exit 1
fi
'''
    
    with open('run_svim_steps.slurm', 'w') as f:
        f.write(script_content)
    
    print(f"Generated SLURM script: run_svim_steps.slurm")

def main():
    parser = create_parser()
    args = parser.parse_args()
    generate_slurm_script(args)

if __name__ == "__main__":
    main()
