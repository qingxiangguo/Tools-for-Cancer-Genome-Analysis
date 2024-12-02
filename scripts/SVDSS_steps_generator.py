#!/usr/bin/env python3
import argparse
import os
from pathlib import Path

def create_parser():
    parser = argparse.ArgumentParser(
        description='Generate SLURM script for SVDSS pipeline')
    
    parser.add_argument('--bam', required=True, type=str,
                      help='Input BAM file path')
    parser.add_argument('--reference', required=True, type=str,
                      help='Reference genome path')
    parser.add_argument('--index', required=True, type=str,
                      help='Reference genome index path (.fmd file)')
    parser.add_argument('--threads', type=int, default=6,
                      help='Number of threads (default: 6)')
    parser.add_argument('--mem', type=str, default='50G',
                      help='Memory allocation (default: 50G)')
    parser.add_argument('--time', type=str, default='06:00:00',
                      help='Wall time limit (format: HH:MM:SS, default: 06:00:00)')
    
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
mamba activate SVDSS

# Step 1: SVDSS smooth
SVDSS smooth \\
    --threads {args.threads} \\
    --bam {os.path.abspath(args.bam)} \\
    --reference {os.path.abspath(args.reference)} > smoothed.bam

# Check if smooth completed successfully
if [ $? -eq 0 ]; then
    # Step 2: Generate BAM index for smoothed.bam
    samtools index smoothed.bam smoothed.bai
    # Check if indexing completed successfully
    if [ $? -eq 0 ]; then
        # Step 3: SVDSS search
        SVDSS search \\
            --threads {args.threads} \\
            --index {os.path.abspath(args.index)} \\
            --bam smoothed.bam > specifics.txt
    else
        echo "BAM indexing failed"
        exit 1
    fi
else
    echo "SVDSS smooth failed"
    exit 1
fi

# Check if search completed successfully
if [ $? -eq 0 ]; then
    # Step 4: SVDSS call
    SVDSS call \\
        --reference {os.path.abspath(args.reference)} \\
        --bam ./smoothed.bam \\
        --sfs specifics.txt > calls.vcf
else
    echo "SVDSS search failed"
    exit 1
fi

if [ $? -ne 0 ]; then
    echo "SVDSS call failed"
    exit 1
fi

# Clean up intermediate files
rm smoothed.bam smoothed.bai specifics.txt
'''
    
    with open('run_svdss_steps.slurm', 'w') as f:
        f.write(script_content)
    
    print(f"Generated SLURM script: run_svdss_steps.slurm")

def main():
    parser = create_parser()
    args = parser.parse_args()
    generate_slurm_script(args)

if __name__ == "__main__":
    main()
