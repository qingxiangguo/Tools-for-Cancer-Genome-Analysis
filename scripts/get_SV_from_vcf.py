#!/usr/bin/env python 
# -*- coding: utf-8 -*-
import sys

def print_help():
    print("Usage: get_SV_from_vcf.py input.vcf INS/DEL/TRA/DUP/INV > output.vcf")
    print("Extract SV type from a VCF file using the specified type.")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Error: wrong number of arguments.")
        print_help()
        sys.exit(1)

# Read the quality threshold from the command line argument
svtype = str(sys.argv[2])

# Open the input VCF file
with open(sys.argv[1], 'r') as f:
    # Iterate through each line of the file
    for line in f:
        # Skip comment lines
        if line.startswith('#'):
            print(line.strip())
        # Split the line into fields
        else:
            fields = line.strip().split('\t')
            # Extract the QUAL field
            if '<' + svtype + '>' == str(fields[4]):
            # Print the line if it pass
                print(line.strip())