#!/usr/bin/env python3

import argparse
import re
from collections import OrderedDict

def parse_fasta_headers(fasta_file):
    """Parse FASTA file to get chromosome names and lengths."""
    chrom_info = OrderedDict()
    current_chrom = None
    length = 0
    
    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                if current_chrom:
                    chrom_info[current_chrom] = length
                current_chrom = line.strip().split()[0][1:]  # Remove '>' and get first word
                length = 0
            else:
                length += len(line.strip())
    
        # Don't forget to add the last chromosome
        if current_chrom:
            chrom_info[current_chrom] = length
            
    return chrom_info

def create_chrom_mapping(chrom_info):
    """Create mapping for common chromosome naming variations."""
    mapping = {}
    for chrom in chrom_info.keys():
        # Remove 'chr' prefix if present
        base_name = chrom.replace('chr', '')
        # Map both with and without 'chr' prefix
        mapping[base_name] = chrom
        mapping[f'chr{base_name}'] = chrom
        # Handle special case for chromosome names
        if base_name.upper() in ['X', 'Y', 'M', 'MT']:
            mapping[base_name.upper()] = chrom
            mapping[base_name.lower()] = chrom
    return mapping

def process_vcf_header(header_lines, chrom_info):
    """Process VCF header lines, adding or updating contig information."""
    # Store existing header lines without contig entries
    new_header_lines = []
    contig_pattern = re.compile(r'##contig=<ID=[^,]+,length=\d+>')
    
    # Keep all non-contig header lines
    for line in header_lines:
        if not contig_pattern.match(line):
            new_header_lines.append(line)
    
    # Add contig lines in the proper location (before #CHROM line)
    final_header_lines = []
    chrom_line_found = False
    
    for line in new_header_lines:
        if line.startswith('#CHROM'):
            # Add all contig lines before the #CHROM line
            for chrom, length in chrom_info.items():
                final_header_lines.append(f'##contig=<ID={chrom},length={length}>')
            chrom_line_found = True
        final_header_lines.append(line)
    
    return final_header_lines

def correct_vcf_file(input_vcf, output_vcf, reference_fasta):
    """Main function to correct VCF file."""
    # Parse reference genome information
    chrom_info = parse_fasta_headers(reference_fasta)
    chrom_mapping = create_chrom_mapping(chrom_info)
    
    header_lines = []
    with open(input_vcf, 'r') as in_f, open(output_vcf, 'w') as out_f:
        # Process header lines
        for line in in_f:
            if line.startswith('#'):
                header_lines.append(line.strip())
                continue
                
            # Process and write header
            if header_lines:
                corrected_header = process_vcf_header(header_lines, chrom_info)
                for header_line in corrected_header:
                    out_f.write(f'{header_line}\n')
                header_lines = None
            
            # Process variant lines
            fields = line.strip().split('\t')
            if len(fields) < 8:  # Minimum VCF fields
                continue
                
            chrom = fields[0]
            if chrom in chrom_mapping and chrom_mapping[chrom] in chrom_info:
                # Update chromosome name to reference format
                fields[0] = chrom_mapping[chrom]
                out_f.write('\t'.join(fields) + '\n')
            elif chrom in chrom_info:
                # Chromosome name already matches reference
                out_f.write('\t'.join(fields) + '\n')

def main():
    parser = argparse.ArgumentParser(description='Correct VCF contig names and add contig headers based on reference genome')
    parser.add_argument('input_vcf', help='Input VCF file')
    parser.add_argument('output_vcf', help='Output VCF file')
    parser.add_argument('reference_fasta', help='Reference FASTA file')
    
    args = parser.parse_args()
    
    correct_vcf_file(args.input_vcf, args.output_vcf, args.reference_fasta)

if __name__ == '__main__':
    main()
