#!/home/qgn1237/2_software/mambaforge/envs/mamba666/bin/python
# _*_ coding=utf-8 _*_
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import sys

# Read the data from the Excel file
excel_file = sys.argv[1]

df = pd.read_excel(excel_file)

# Set the color palette for the plot
palette = sns.color_palette("Set3", n_colors=6)

# Plot the data as a bar plot
sns.barplot(data=df, x='depth', y='SV_no', hue='cell_line', palette=palette)

plt.xlabel('Sequencing depth', fontsize=18)

plt.ylabel('Number of variations', fontsize=18)

plt.title('Number of SVs in cell lines (PacBio) across depths', fontsize=22)

plt.xticks(size = 16)

plt.yticks(size = 16)

plt.legend(fontsize=18)

# Save the plot as a PNG file with the desired size
plt.gcf().set_size_inches(14, 8)

# Save the plot as a PNG file
plt.savefig('plot.png', dpi=900)