#!/home/qgn1237/2_software/mambaforge/envs/mamba666/bin/python
# -*- coding: utf-8 -*-
import sys

def print_help():
    print("Usage: filter_vcf_based_on_quality.py input.vcf 20 > output.vcf")
    print("Filters a VCF file using the specified threshold.")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Error: wrong number of arguments.")
        print_help()
        sys.exit(1)

# Read the quality threshold from the command line argument
threshold = float(sys.argv[2])

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
            qual = float(fields[5])
        # Check if the QUAL is above the threshold
            if qual >= threshold:
            # Print the line if it passes the filter
                print(line.strip())