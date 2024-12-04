#!/usr/bin/env python3
import argparse
import os
from pathlib import Path

def create_parser():
    parser = argparse.ArgumentParser(
        description='Generate SLURM script for ONT Direct RNA004 analysis pipeline')
    
    # Required arguments
    parser.add_argument('--bam', required=True, type=str,
                      help='Input BAM file from Dorado basecalling')
    parser.add_argument('--reference', required=True, type=str,
                      help='Reference genome path (minimap2 index)')
    
    # Optional arguments with defaults
    parser.add_argument('--min_quality', type=int, default=8,
                      help='Minimum read quality score (default: 8)')
    parser.add_argument('--min_length', type=int, default=500,
                      help='Minimum read length in bp (default: 500)')
    parser.add_argument('--threads', type=int, default=8,
                      help='Number of threads (default: 8)')
    parser.add_argument('--mem', type=str, default='70G',
                      help='Memory allocation (default: 70G)')
    parser.add_argument('--time', type=str, default='48:00:00',
                      help='Wall time limit (format: HH:MM:SS, default: 48:00:00)')
    parser.add_argument('--sample_name', type=str, default='sample',
                      help='Sample name for output files (default: sample)')
    
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

# Step 1: Convert BAM to FASTQ
echo "Converting BAM to FASTQ..."
samtools view {os.path.abspath(args.bam)} | samtools fastq > dorado.fastq

# Check if conversion completed successfully
if [ $? -ne 0 ]; then
    echo "BAM to FASTQ conversion failed"
    exit 1
fi

# Step 2: Quality filtering with chopper
echo "Filtering reads with chopper..."
cat dorado.fastq | chopper -t {args.threads} \\
    -q {args.min_quality} \\
    -l {args.min_length} > clean_reads.fastq

if [ $? -ne 0 ]; then
    echo "Chopper filtering failed"
    exit 1
fi

# Step 3: Align with minimap2 and sort with samtools
echo "Aligning reads to reference genome..."
minimap2 -ax splice -uf -k14 --MD -t {args.threads} -Y \\
    -R '@RG\\tID:{args.sample_name}\\tPL:ont\\tLB:library\\tSM:{args.sample_name}' \\
    {os.path.abspath(args.reference)} clean_reads.fastq | \\
    samtools sort -@ {args.threads} -m 2G -O BAM \\
    -o {args.sample_name}_direct_RNA.bam

if [ $? -ne 0 ]; then
    echo "Alignment or sorting failed"
    exit 1
fi

# Step 4: Index BAM file
echo "Indexing BAM file..."
samtools index {args.sample_name}_direct_RNA.bam {args.sample_name}_direct_RNA.bam.bai

if [ $? -ne 0 ]; then
    echo "BAM indexing failed"
    exit 1
fi

# Step 5: Quality control with NanoPlot
echo "Running NanoPlot QC on cleaned reads..."
NanoPlot -t {args.threads} --fastq clean_reads.fastq --format pdf \\
    -o nanoplot_qc_results

if [ $? -ne 0 ]; then
    echo "NanoPlot QC failed"
    exit 1
fi

# Step 6: Quality control with Qualimap
echo "Running Qualimap QC on aligned BAM..."
qualimap bamqc \\
    -bam {args.sample_name}_direct_RNA.bam \\
    -c ./qualimap \\
    --java-mem-size={args.mem}

if [ $? -ne 0 ]; then
    echo "Qualimap QC failed"
    exit 1
fi

echo "Pipeline completed successfully!"
'''
    
    output_script = f'run_ont_rna_{args.sample_name}.slurm'
    with open(output_script, 'w') as f:
        f.write(script_content)
    
    print(f"Generated SLURM script: {output_script}")

def main():
    parser = create_parser()
    args = parser.parse_args()
    generate_slurm_script(args)

if __name__ == "__main__":
    main()
