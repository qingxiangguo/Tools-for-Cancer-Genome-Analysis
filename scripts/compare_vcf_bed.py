#!/usr/bin/env python3

import sys
from collections import defaultdict

def parse_vcf(vcf_file):
    """Parse VCF file and count unique SVs by type and position"""
    sv_counts = defaultdict(int)
    unique_positions = set()
    
    with open(vcf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            chrom = fields[0]
            pos = fields[1]
            info = dict(item.split('=', 1) if '=' in item else (item, True) 
                       for item in fields[7].split(';'))
            
            svtype = info.get('SVTYPE')
            if svtype:
                key = f"{chrom}_{pos}_{svtype}"
                unique_positions.add(key)
                sv_counts[svtype] += 1
                
    return sv_counts, len(unique_positions)

def parse_bed(bed_file):
    """Parse BED file and count SVs by type"""
    sv_counts = defaultdict(int)
    unique_positions = set()
    
    with open(bed_file, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            chrom = fields[0]
            start = fields[1]
            sv_type = fields[3]
            
            # Convert BED types to VCF types
            type_mapping = {
                'deletion': 'DEL',
                'insertion': 'INS',
                'inversion': 'INV',
                'tandem duplication': 'DUP',
                'translocation copy-paste': 'TRA'
            }
            
            sv_type = type_mapping.get(sv_type, sv_type)
            key = f"{chrom}_{start}_{sv_type}"
            unique_positions.add(key)
            sv_counts[sv_type] += 1
            
    return sv_counts, len(unique_positions)

def main():
    vcf_file = "nstd106.GRCh38.variant_true_call.vcf"
    bed1_file = "visor_bed_haplotype1.bed"
    bed2_file = "visor_bed_haplotype2.bed"
    
    # Parse files
    vcf_counts, vcf_unique = parse_vcf(vcf_file)
    bed1_counts, bed1_unique = parse_bed(bed1_file)
    bed2_counts, bed2_unique = parse_bed(bed2_file)
    
    # Combine bed1 and bed2 counts for homozygous variants
    total_bed_counts = defaultdict(int)
    for sv_type in set(list(bed1_counts.keys()) + list(bed2_counts.keys())):
        total_bed_counts[sv_type] = bed1_counts[sv_type] + bed2_counts[sv_type]
    
    # Print comparison
    print("\nComparison Results:")
    print("-" * 50)
    print(f"Total unique positions in VCF: {vcf_unique}")
    print(f"Total unique positions in BED1: {bed1_unique}")
    print(f"Total unique positions in BED2: {bed2_unique}")
    print("\nCounts by SV type:")
    print("Type\tVCF\tBED_Total")
    for sv_type in sorted(set(list(vcf_counts.keys()) + list(total_bed_counts.keys()))):
        print(f"{sv_type}\t{vcf_counts[sv_type]}\t{total_bed_counts[sv_type]}")

if __name__ == "__main__":
    main()
