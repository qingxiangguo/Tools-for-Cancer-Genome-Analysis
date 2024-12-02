#!/usr/bin/env python3
import argparse
import os
from pathlib import Path

def create_parser():
    parser = argparse.ArgumentParser(
        description='Generate SLURM script for SVABA pipeline')
    
    parser.add_argument('--bam', required=True, type=str,
                      help='Input BAM file path')
    parser.add_argument('--reference', required=True, type=str,
                      help='Reference genome path')
    parser.add_argument('--tag', required=True, type=str,
                      help='Analysis tag for output files')
    parser.add_argument('--threads', type=int, default=8,
                      help='Number of threads (default: 8)')
    parser.add_argument('--mem', type=str, default='50G',
                      help='Memory allocation (default: 50G)')
    parser.add_argument('--time', type=str, default='16:00:00',
                      help='Wall time limit (format: HH:MM:SS, default: 16:00:00)')
    
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

# Initialize mamba
/home/qgn1237/2_software/mambaforge/bin/mamba init
source ~/.bashrc
mamba activate mamba666

# Run SVABA
svaba run -t {os.path.abspath(args.bam)} \\
    -p {args.threads} \\
    -a {args.tag} \\
    -G {os.path.abspath(args.reference)}

if [ $? -eq 0 ]; then
    echo "SVABA analysis completed successfully"
else
    echo "SVABA analysis failed"
    exit 1
fi
'''
    
    with open('run_svaba_steps.slurm', 'w') as f:
        f.write(script_content)
    
    print(f"Generated SLURM script: run_svaba_steps.slurm")

def main():
    parser = create_parser()
    args = parser.parse_args()
    generate_slurm_script(args)

if __name__ == "__main__":
    main()
