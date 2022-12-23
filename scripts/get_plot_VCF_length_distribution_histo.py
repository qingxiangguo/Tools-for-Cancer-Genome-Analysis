#!/home/qgn1237/2_software/mambaforge/envs/mamba666/bin/python
# _*_ coding=utf-8 _*_
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import argparse

# Set arguments
parser = argparse.ArgumentParser(description = 'Process VCF file and plot the length distribution of your defined variant type')
parser.add_argument('-i', '--input', required = True, metavar='', help='VCF file to process')
parser.add_argument('-v', '--variant', required=True, metavar='', help='Variant type to extract, e.g., INS, DEL, INV, TRA, DUP')

args = parser.parse_args()
input_file = args.input
variant_type = args.variant 

# Initialize two empty lists
name_list = []
length_list = []

with open(input_file, 'r') as vcf_file:
    for line in vcf_file:
        if line[0] == '#': # Skip annotation lines
            continue
        fields = line.strip().split('\t')
        info_field = fields[7]
        if f'SVTYPE={variant_type}' in info_field:
            name_list.append(fields[2])
            # Should take care of the type, int and str
            length = abs(int(info_field.split('SVLEN=')[1].split(';')[0])) # Split and split, get the length
            length_list.append(length)
            
# Prepare pandas dataframe from previous list
data = {'SV_name' : name_list, 'SV_length' : length_list}
df = pd.DataFrame.from_dict(data) # Create pandas dataframe from dictionary

# Start to plot using seaborn
if variant_type == 'DEL':
    sns.histplot(data=df, x='SV_length', color='#1f77b4')
elif variant_type == 'INS':
    sns.histplot(data=df, x='SV_length', color='#2ba02b') 
elif variant_type == 'INV':
    sns.histplot(data=df, x='SV_length', color='#d62728')             
elif variant_type == 'TRA':
    sns.histplot(data=df, x='SV_length', color='#9467bd') 
elif variant_type == 'DUP':
    sns.histplot(data=df, x='SV_length', color='#ff7f0f')     
else:
    print("The variant type must be one of the INS, DEL, INV, TRA, DUP")

# Save the plot as a PNG file with the desired size
plt.gcf().set_size_inches(14, 8)

# Save the plot as a PNG file
plt.savefig('plot.png', dpi=900)