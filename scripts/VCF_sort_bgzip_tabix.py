#!/usr/bin/env python3

import argparse
import os
import subprocess
import sys

def get_prefix(vcf_path):
    """
    Extract the prefix from input VCF file path
    Example: /path/to/test.vcf -> test
    """
    base = os.path.basename(vcf_path)  # Get filename from path
    prefix = os.path.splitext(base)[0]  # Remove extension
    return prefix

def sort_and_compress_vcf(input_vcf):
    """
    Process VCF file:
    1. Sort VCF while preserving headers
    2. Compress with bgzip
    3. Index with tabix
    
    Args:
        input_vcf (str): Path to input VCF file
    """
    prefix = get_prefix(input_vcf)
    sorted_vcf = f"{prefix}_sorted.vcf"
    compressed_vcf = f"{prefix}_sorted.vcf.gz"
    
    # Command to sort VCF while preserving headers
    sort_cmd = f"cat {input_vcf} | awk '$1 ~ /^#/ {{print $0;next}} {{print $0 | \"sort -k1,1V -k2,2n\"}}' > {sorted_vcf}"
    
    # Command to compress with bgzip
    bgzip_cmd = f"bgzip -c {sorted_vcf} > {compressed_vcf}"
    
    # Command to index with tabix
    tabix_cmd = f"tabix -p vcf {compressed_vcf}"
    
    try:
        # Execute sort command
        print(f"Sorting VCF file: {input_vcf}")
        subprocess.run(sort_cmd, shell=True, check=True)
        
        # Execute bgzip command
        print(f"Compressing sorted VCF: {sorted_vcf}")
        subprocess.run(bgzip_cmd, shell=True, check=True)
        
        # Execute tabix command
        print(f"Creating tabix index for: {compressed_vcf}")
        subprocess.run(tabix_cmd, shell=True, check=True)
        
        print("\nProcess completed successfully!")
        print(f"Generated files:")
        print(f"1. Sorted VCF: {sorted_vcf}")
        print(f"2. Compressed VCF: {compressed_vcf}")
        print(f"3. Tabix index: {compressed_vcf}.tbi")
        
    except subprocess.CalledProcessError as e:
        print(f"Error occurred during execution: {e}", file=sys.stderr)
        sys.exit(1)

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description='Sort, compress, and index VCF file')
    parser.add_argument('-i', '--input', required=True, help='Input VCF file')
    
    # Parse arguments
    args = parser.parse_args()
    
    # Check if input file exists
    if not os.path.exists(args.input):
        print(f"Error: Input file {args.input} does not exist", file=sys.stderr)
        sys.exit(1)
    
    # Process the VCF file
    sort_and_compress_vcf(args.input)

if __name__ == '__main__':
    main()
