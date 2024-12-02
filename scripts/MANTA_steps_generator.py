#!/usr/bin/env python3
import argparse
import os
from pathlib import Path

def create_parser():
    parser = argparse.ArgumentParser(
        description='Generate SLURM script for MANTA pipeline')
    
    parser.add_argument('--bam', required=True, type=str,
                      help='Input BAM file path')
    parser.add_argument('--reference', required=True, type=str,
                      help='Reference genome path')
    parser.add_argument('--work-dir', required=True, type=str,
                      help='Working directory for MANTA output')
    parser.add_argument('--threads', type=int, default=8,
                      help='Number of threads (default: 8)')
    parser.add_argument('--mem', type=str, default='80G',
                      help='Memory allocation (default: 80G)')
    parser.add_argument('--time', type=str, default='48:00:00',
                      help='Wall time limit (format: HH:MM:SS, default: 48:00:00)')
    
    return parser

def generate_slurm_script(args):
    # Get absolute path for work directory
    work_dir = os.path.abspath(args.work_dir)
    
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
mamba activate mamba_py2

# Step 1: Run configManta.py
configManta.py \\
    --tumorBam={os.path.abspath(args.bam)} \\
    --referenceFasta={os.path.abspath(args.reference)} \\
    --runDir={work_dir}

# Check if configManta.py was successful
if [ $? -eq 0 ]; then
    echo "configManta.py completed successfully"
    
    # Step 2: Run the workflow script
    {work_dir}/runWorkflow.py -j {args.threads}
    
    if [ $? -eq 0 ]; then
        echo "MANTA workflow completed successfully"
        
        # Step 3: Clean up and organize results
        cd {work_dir}
        
        # Move variant results to working directory
        mv results/variants/* ./
        
        # Clean up intermediate and unnecessary files
        rm -rf results/ workspace/
        rm runWorkflow.py* workf*
        
        # Decompress the VCF file
        gunzip tumorSV.vcf.gz
        
        echo "Results cleaned and organized successfully"
    else
        echo "MANTA workflow failed"
        exit 1
    fi
else
    echo "configManta.py failed"
    exit 1
fi
'''
    
    with open('run_manta_steps.slurm', 'w') as f:
        f.write(script_content)
    
    print(f"Generated SLURM script: run_manta_steps.slurm")

def main():
    parser = create_parser()
    args = parser.parse_args()
    
    # Create work directory if it doesn't exist
    work_dir = Path(args.work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)
    
    generate_slurm_script(args)

if __name__ == "__main__":
    main()
