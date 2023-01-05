#!/home/qgn1237/2_software/mambaforge/envs/mamba_py_39/bin/python
# _*_ coding=utf-8 _*_
from collections import Counter
import pandas as pd
import matplotlib.pyplot as plt
import argparse
import os
import re

# Set arguments
parser = argparse.ArgumentParser(description='Process VCF file and plot the normalized stacked frequency of defined SV in chromosomes')
parser.add_argument('-i', '--input', required=True, metavar='', help='VCF file to process')
args = parser.parse_args()
input_file = args.input

# get the path of input file
path = os.path.abspath(input_file)
print(path)
match = re.search(r"(\w+)/(\w+)/SURVIVOR", path)
if match:
    cell_name = match.group(1)
    depth_name = match.group(2)
else:
    print("Not matched")

# Set list that will be used later
chromo_correct_list = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6',
                       'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12',
                       'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18',
                       'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']

chr_len_list = [237.1, 231.5, 188.9, 180.9, 173.3, 162.5,
                152.5, 137.5, 131.9, 127.4, 128.5, 126.4,
                109.1, 102.2, 96.9, 86.0, 79.1, 77.0,
                56.3, 61.5, 44.7, 48.3, 148.1, 54.7]

del_raw_list = []
dup_raw_list = []
ins_raw_list = []
inv_raw_list = []
tra_raw_list = []


# Define a function to judge
def whether_in_correct(chr_str):
    if chr_str in chromo_correct_list:
        return True
    else:
        return False


with open(input_file, "r") as vcf_file:
    for line in vcf_file:
        if line[0] == "#":
            continue
        fields = line.strip().split('\t')
        info = fields[7]
        if f'SVTYPE=DEL' in info:
            if whether_in_correct(fields[0]):
                del_raw_list.append(fields[0])
            else:
                continue
        elif f'SVTYPE=DUP' in info:
            if whether_in_correct(fields[0]):
                dup_raw_list.append(fields[0])
            else:
                continue
        elif f'SVTYPE=INS' in info:
            if whether_in_correct(fields[0]):
                ins_raw_list.append(fields[0])
            else:
                continue
        elif f'SVTYPE=INV' in info:
            if whether_in_correct(fields[0]):
                inv_raw_list.append(fields[0])
            else:
                continue
        elif f'SVTYPE=TRA' in info:
            if whether_in_correct(fields[0]):
                tra_raw_list.append(fields[0])
            else:
                continue


# Start to transform raw list
def count_frequencies(lst):
    frequencies = [lst.count(element) for element in chromo_correct_list] 
    return frequencies


def get_chr_normalize_list(lst):
        freq_list = count_frequencies(lst)
        new_list = [x / y for x, y in zip(freq_list, chr_len_list)]
        return new_list


del_normal_list = get_chr_normalize_list(del_raw_list)
dup_normal_list = get_chr_normalize_list(dup_raw_list)
ins_normal_list = get_chr_normalize_list(ins_raw_list)
inv_normal_list = get_chr_normalize_list(inv_raw_list)
tra_normal_list = get_chr_normalize_list(tra_raw_list)

# Start to make chromosome stack plot
index_list = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6',
              'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12',
              'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18',
              'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']

# Create pandas dataframe
df = pd.DataFrame({"Deletion": del_normal_list, "Duplication": dup_normal_list,
                   "Insertion": ins_normal_list, "Inversion": inv_normal_list,
                   "Translocation": tra_normal_list}, index=index_list)
# Create stacked barplot, rot=0, means do not rotate the x axis
df.plot.bar(stacked=True, rot=0, width=0.6)

plt.xlabel('Chromosome', fontsize=18)
plt.ylabel('Normalized ratio (No. of SV events/Mbp of sequences)', fontsize=18)

plt.title(f'Relative chromosomal distribution of the SV events for {cell_name} {depth_name}', fontsize=22)

plt.xticks(size=16, rotation=30)
plt.yticks(size=16)

plt.legend(fontsize=18)

# Save the plot as a PNG file with the desired size
# plt.gcf use the plt from matplotlib.pyplot to get the picture self and then resize
plt.gcf().set_size_inches(16, 8)

# Save the plot as a PNG file
plt.savefig(f"{cell_name}_{depth_name}_normalized_chromosome.png", dpi=900)
