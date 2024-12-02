#!/usr/bin/env python3
import argparse
import os
from pathlib import Path

def create_parser():
    parser = argparse.ArgumentParser(
        description='Generate SLURM script for Debreak pipeline')
    
    parser.add_argument('--bam', required=True, type=str,
                      help='Input BAM file path')
    parser.add_argument('--reference', required=True, type=str,
                      help='Reference genome path')
    parser.add_argument('--threads', type=int, default=8,
                      help='Number of threads (default: 8)')
    parser.add_argument('--min-support', type=int, default=2,
                      help='Minimum support reads (default: 2)')
    parser.add_argument('--mem', type=str, default='40G',
                      help='Memory allocation (default: 40G)')
    parser.add_argument('--time', type=str, default='05:00:00',
                      help='Wall time limit (format: HH:MM:SS, default: 05:00:00)')
    
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
mamba activate debreak

# Run Debreak
debreak \\
    --min_support {args.min_support} \\
    -t {args.threads} \\
    --bam {os.path.abspath(args.bam)} \\
    -o ./ \\
    --rescue_large_ins \\
    --rescue_dup \\
    --poa \\
    --ref {os.path.abspath(args.reference)}

# Check if Debreak completed successfully
if [ $? -eq 0 ]; then
    # Clean up temporary files and directories
    rm -rf debreak-* \\
          debreak_* \\
          map_depth/ \\
          sv_raw_calls/ \\
          deletion-merged \\
          duplication-merged \\
          insertion-merged \\
          inversion-merged \\
          translocation-merged
else
    echo "Debreak analysis failed"
    exit 1
fi
'''
    
    with open('run_debreak_steps.slurm', 'w') as f:
        f.write(script_content)
    
    print(f"Generated SLURM script: run_debreak_steps.slurm")

def main():
    parser = create_parser()
    args = parser.parse_args()
    generate_slurm_script(args)

if __name__ == "__main__":
    main()
