#!/usr/bin/env python3

import sys
import argparse
import re
import random
from collections import defaultdict

def parse_args():
    parser = argparse.ArgumentParser(description="Convert custom VCF to VISOR HACk BED format with zygosity assignment.")
    parser.add_argument("-i", "--input_vcf", required=True, help="Input VCF file.")
    parser.add_argument("-o", "--output_prefix", required=True, help="Output prefix for BED files (haplotype1 and haplotype2).")
    parser.add_argument("--homozygous_ratio", type=float, default=0.3, help="Ratio of homozygous variants (default: 0.3).")
    return parser.parse_args()

def generate_random_sequence(length):
    """Generate random DNA sequence of specified length"""
    return ''.join(random.choices('ATCG', k=length))

def parse_translocation(alt):
    """Parse different types of translocation notation in ALT field"""
    # Pattern 1: C[12:10767[ or A[1:224014549[
    pattern1 = r'[ATCG]?\[([^:]+):(\d+)\['
    # Pattern 2: [Y:10136631[T
    pattern2 = r'\[([^:]+):(\d+)\][ATCG]?'
    # Pattern 3: ]4:71592079]G
    pattern3 = r'\]([^:]+):(\d+)\][ATCG]?'
    
    for pattern in [pattern1, pattern2, pattern3]:
        match = re.search(pattern, alt)
        if match:
            target_chrom, target_pos = match.groups()
            return target_chrom, int(target_pos)
    return None, None

def process_insertion(pos, info_dict, alt):
    """Process insertion and return appropriate sequence"""
    # 如果有SEQ字段，直接使用
    if 'SEQ' in info_dict:
        return info_dict['SEQ']
    
    # 如果有SVLEN，使用SVLEN生成随机序列
    if 'SVLEN' in info_dict and info_dict['SVLEN'] != '.' and info_dict['SVLEN'] != '0':
        try:
            length = abs(int(info_dict['SVLEN']))
            return generate_random_sequence(length)
        except ValueError:
            pass
    
    # 如果有END，使用END-POS+1作为长度生成随机序列
    if 'END' in info_dict and info_dict['END'] != '.':
        try:
            end = int(info_dict['END'])
            length = end - int(pos) + 1
            return generate_random_sequence(length)
        except ValueError:
            pass
            
    # 默认使用100bp长度
    return generate_random_sequence(100)

def main():
    args = parse_args()
    stats = defaultdict(int)
    
    hap1_bed = args.output_prefix + '_haplotype1.bed'
    hap2_bed = args.output_prefix + '_haplotype2.bed'
    log_file = args.output_prefix + '_conversion.log'
    
    with open(hap1_bed, 'w') as hap1_out, \
         open(hap2_bed, 'w') as hap2_out, \
         open(log_file, 'w') as log:
        
        with open(args.input_vcf, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                
                stats['total'] += 1
                fields = line.strip().split('\t')
                chrom = fields[0]
                pos = int(fields[1])
                ref = fields[3]
                alt = fields[4]
                info = dict(item.split('=', 1) if '=' in item else (item, True) 
                           for item in fields[7].split(';'))
                
                svtype = info.get('SVTYPE')
                if not svtype:
                    log.write(f"Warning: No SVTYPE found in line: {line}")
                    continue
                
                # 分配单倍型
                if random.random() < args.homozygous_ratio:
                    haplotype = 'both'
                    stats['homozygous'] += 1
                else:
                    haplotype = random.choice(['hap1', 'hap2'])
                    stats['heterozygous'] += 1
                
                # 处理不同类型的SV
                bed_entry = None
                
                if svtype == 'DEL':
                    end = info.get('END', str(pos + abs(int(info.get('SVLEN', '100')))))
                    bed_entry = [chrom, str(pos), end, 'deletion', 'None', '0']
                
                elif svtype == 'INS':
                    insertion_seq = process_insertion(pos, info, alt)
                    bed_entry = [chrom, str(pos), str(pos+1), 'insertion', insertion_seq, '0']
                
                elif svtype == 'INV':
                    end = info.get('END', str(pos + abs(int(info.get('SVLEN', '100')))))
                    bed_entry = [chrom, str(pos), end, 'inversion', 'None', '0']
                
                elif svtype == 'DUP':
                    end = info.get('END', str(pos + abs(int(info.get('SVLEN', '100')))))
                    bed_entry = [chrom, str(pos), end, 'tandem duplication', '2', '0']
                
                elif svtype == 'TRA':
                    target_chrom, target_pos = parse_translocation(alt)
                    if target_chrom and target_pos:
                        hap = 'h1' if haplotype == 'hap1' else 'h2'
                        bed_entry = [chrom, str(pos), str(pos), 
                                   'translocation copy-paste',
                                   f"{hap}:{target_chrom}:{target_pos}:forward", 
                                   '0']
                
                if bed_entry:
                    bed_line = '\t'.join(bed_entry) + '\n'
                    stats['processed'] += 1
                    stats[f'processed_{svtype}'] += 1
                    
                    if haplotype == 'both':
                        hap1_out.write(bed_line)
                        hap2_out.write(bed_line)
                    elif haplotype == 'hap1':
                        hap1_out.write(bed_line)
                    else:
                        hap2_out.write(bed_line)
                else:
                    log.write(f"Warning: Could not process line: {line}")
    
    # 写入统计信息到日志文件
    with open(log_file, 'a') as log:
        log.write(f"\nProcessing Statistics:\n")
        log.write(f"Total variants: {stats['total']}\n")
        log.write(f"Processed variants: {stats['processed']}\n")
        log.write(f"Homozygous variants: {stats['homozygous']}\n")
        log.write(f"Heterozygous variants: {stats['heterozygous']}\n")
        log.write("\nProcessed by type:\n")
        for key in stats:
            if key.startswith('processed_'):
                log.write(f"{key.replace('processed_', '')}: {stats[key]}\n")

if __name__ == "__main__":
    main()
