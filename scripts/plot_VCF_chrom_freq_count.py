#!/home/qgn1237/2_software/mambaforge/envs/mamba666/bin/python
# _*_ coding=utf-8 _*_
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import argparse
import os
import re

# Set arguments
parser = argparse.ArgumentParser(description='Process VCF file and plot the frequency of defined SV in chromosomes') 
parser.add_argument('-i', '--input', required=True, metavar='', help='VCF file to process')
parser.add_argument('-v', '--variant', required=True, metavar='', help='Variant type to extract, e.g., INS, DEL, INV, DUP, TRA')

args = parser.parse_args()
input_file = args.input
variant_type = args.variant

#get the path of input file
path = os.path.abspath(input_file)
print(path)
match = re.search(r"(\w+)/(\w+)/SURVIVOR", path)
if match:
    cell_name = match.group(1)
    depth_name = match.group(2)
else:
    print("Not matched")

# Initialize one empty lists
chromo_distribution_list = []
chromo_correct_list = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6',
                       'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12',
                       'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18',
                       'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']

# Def a function to judge
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
        if f'SVTYPE={variant_type}' in info:
            if whether_in_correct(fields[0]):
                chromo_distribution_list.append(fields[0])
            else:
                continue

# Prepare pandas dataframe from previous list
data = {'chr_info' : chromo_distribution_list}
df = pd.DataFrame.from_dict(data)

# Start to plot using seaborn
if variant_type == 'DEL':
    sns.countplot(data=df, x='chr_info', color='#1f77b4', order=chromo_correct_list)
    plt.title(f'Chromosome distribution of deletion for {cell_name} {depth_name}', fontsize=22)
elif variant_type == 'INS':
    sns.countplot(data=df, x='chr_info', color='#2ba02b', order=chromo_correct_list) 
    plt.title(f'Chromosome distribution of insertion for {cell_name} {depth_name}', fontsize=22)
elif variant_type == 'INV':
    sns.countplot(data=df, x='chr_info', color='#d62728', order=chromo_correct_list)             
    plt.title(f'Chromosome distribution of inversion for {cell_name} {depth_name}', fontsize=22)
elif variant_type == 'DUP':
    sns.countplot(data=df, x='chr_info', color='#ff7f0f', order=chromo_correct_list)     
    plt.title(f'Chromosome distribution of duplication for {cell_name} {depth_name}', fontsize=22)
elif variant_type == 'TRA':
    sns.countplot(data=df, x='chr_info', color='#9467bd', order=chromo_correct_list)     
    plt.title(f'Chromosome distribution of translocation for {cell_name} {depth_name}', fontsize=22)
else:
    print("The variant type must be one of the INS, DEL, INV, DUP, TRA")
        
plt.xlabel('Chromosome', fontsize=16)

plt.ylabel('Number of variations', fontsize=16)

plt.xticks(size = 14, rotation=30)

plt.yticks(size = 14)

# Save the plot as a PNG file with the desired size
plt.gcf().set_size_inches(14, 9)

# Save the plot as a PNG file
plt.savefig(f"{cell_name}_{depth_name}_{variant_type}_chromosome.png", dpi=900)        

    