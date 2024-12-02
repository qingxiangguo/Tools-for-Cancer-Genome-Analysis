#!/usr/bin/env python3
import argparse
import os
from pathlib import Path

def create_parser():
    parser = argparse.ArgumentParser(
        description='Generate SLURM script for LUMPY pipeline')
    
    parser.add_argument('--bam', required=True, type=str,
                      help='Input BAM file path')
    parser.add_argument('--output', type=str, default='lumpy.vcf',
                      help='Output VCF file name (default: lumpy.vcf)')
    parser.add_argument('--ntasks', type=int, default=1,
                      help='Number of tasks (default: 1)')
    parser.add_argument('--mem', type=str, default='50G',
                      help='Memory allocation (default: 50G)')
    parser.add_argument('--time', type=str, default='15:00:00',
                      help='Wall time limit (format: HH:MM:SS, default: 15:00:00)')
    
    return parser

def generate_slurm_script(args):
    script_content = f'''#!/bin/bash -l
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --nodes=1
#SBATCH --ntasks={args.ntasks}
#SBATCH --mem={args.mem}
#SBATCH --time={args.time}
#SBATCH --mail-type=ALL

cd $SLURM_SUBMIT_DIR

# Initialize mamba
/home/qgn1237/2_software/mambaforge/bin/mamba init
source ~/.bashrc
mamba activate lumpy_env

# Run LUMPY
lumpyexpress -B {os.path.abspath(args.bam)} -o {args.output}

if [ $? -eq 0 ]; then
    echo "LUMPY analysis completed successfully"
else
    echo "LUMPY analysis failed"
    exit 1
fi
'''
    
    with open('run_lumpy_steps.slurm', 'w') as f:
        f.write(script_content)
    
    print(f"Generated SLURM script: run_lumpy_steps.slurm")

def main():
    parser = create_parser()
    args = parser.parse_args()
    generate_slurm_script(args)

if __name__ == "__main__":
    main()
