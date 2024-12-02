#!/usr/bin/env python3
import argparse
import os
from pathlib import Path

def create_parser():
    parser = argparse.ArgumentParser(
        description='Generate SLURM script for CuteSV pipeline')
    
    parser.add_argument('--bam', required=True, type=str,
                      help='Input BAM file path')
    parser.add_argument('--reference', required=True, type=str,
                      help='Reference genome path')
    parser.add_argument('--data-type', required=True, choices=['pacbio', 'ont'],
                      help='Sequencing data type: pacbio (HIFI/CCS) or ont')
    parser.add_argument('--threads', type=int, default=16,
                      help='Number of threads (default: 16)')
    parser.add_argument('--mem', type=str, default='50G',
                      help='Memory allocation (default: 50G)')
    parser.add_argument('--time', type=str, default='03:00:00',
                      help='Wall time limit (format: HH:MM:SS, default: 03:00:00)')
    parser.add_argument('--output', type=str, default='cutesv.vcf',
                      help='Output VCF filename (default: cutesv.vcf)')
    
    return parser

def get_cutesv_params(data_type):
    """Return CuteSV parameters based on data type"""
    params = {
        'pacbio': {
            'max_cluster_bias_INS': 1000,
            'diff_ratio_merging_INS': 0.9,
            'max_cluster_bias_DEL': 1000,
            'diff_ratio_merging_DEL': 0.5
        },
        'ont': {
            'max_cluster_bias_INS': 100,
            'diff_ratio_merging_INS': 0.3,
            'max_cluster_bias_DEL': 100,
            'diff_ratio_merging_DEL': 0.3
        }
    }
    return params[data_type]

def generate_slurm_script(args):
    # Get data type specific parameters
    params = get_cutesv_params(args.data_type)
    
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

# Run CuteSV with {args.data_type.upper()} specific parameters
cuteSV \\
    --threads {args.threads} \\
    --genotype \\
    --report_readid \\
    -l 50 \\
    -L 5000000 \\
    -r 1000 \\
    -q 20 \\
    -s 3 \\
    --max_cluster_bias_INS {params['max_cluster_bias_INS']} \\
    --diff_ratio_merging_INS {params['diff_ratio_merging_INS']} \\
    --max_cluster_bias_DEL {params['max_cluster_bias_DEL']} \\
    --diff_ratio_merging_DEL {params['diff_ratio_merging_DEL']} \\
    {os.path.abspath(args.bam)} \\
    {os.path.abspath(args.reference)} \\
    {args.output} \\
    ./

if [ $? -ne 0 ]; then
    echo "CuteSV analysis failed"
    exit 1
fi
'''
    
    with open('run_cutesv_steps.slurm', 'w') as f:
        f.write(script_content)
    
    print(f"Generated SLURM script: run_cutesv_steps.slurm")

def main():
    parser = create_parser()
    args = parser.parse_args()
    generate_slurm_script(args)

if __name__ == "__main__":
    main()
