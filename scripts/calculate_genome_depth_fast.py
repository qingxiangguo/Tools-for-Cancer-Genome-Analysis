#!/home/qgn1237/2_software/mambaforge/envs/mamba666/bin/python
# _*_ coding=utf-8 _*_
import pysam
import sys

# read in the BAM file as the first argument
bam_file = sys.argv[1]

# open the BAM file using pysam
samfile = pysam.AlignmentFile(bam_file, "rb")

# initialize a variable to store the total length of the mapped reads
total_length = 0

for read in samfile.fetch():
    if read.is_mapped:
        # add the query length of the mapped read to the total length
        total_length += read.query_length

# calculate the genome coverage:
genome_cov = total_length/3209286105

with open(f"{bam_file}.txt", "w") as outfile:
    outfile.write("The genome coverage is {}".format(genome_cov))