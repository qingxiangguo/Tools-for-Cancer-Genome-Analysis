#!/usr/bin/env python 
# -*- coding: utf-8 -*-
import argparse
import os
import sys

def parse_args():
    parser = argparse.ArgumentParser(description='Filter VCF file by SV length.')
    parser.add_argument('-i', '--input', required=True, help='Input VCF file.')
    parser.add_argument('-o', '--output', required=True, help='Output VCF file.')
    parser.add_argument('-l', '--length', type=int, required=True, help='Minimum SV length.')
    return parser.parse_args()

def filter_vcf(input_file, output_file, min_length):
    with open(input_file, 'r') as fin, open(output_file, 'w') as fout:
        for line in fin:
            if line.startswith('#'):
                fout.write(line)
                continue
            info = {x.split('=')[0]:x.split('=')[1] for x in line.split('\t')[7].split(';') if '=' in x}
            svlen = info.get('SVLEN')
            svtype = info.get('SVTYPE')
            if svlen is not None:
                svlen = abs(int(svlen))
            if svtype in ['BND', 'TRA'] or svlen is None or svlen >= min_length:
                fout.write(line)

if __name__ == '__main__':
    args = parse_args()
    filter_vcf(args.input, args.output, args.length)
