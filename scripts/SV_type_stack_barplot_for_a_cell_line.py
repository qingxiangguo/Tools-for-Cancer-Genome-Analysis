#!/home/qgn1237/2_software/mambaforge/envs/mamba666/bin/python
# _*_ coding=utf-8 _*_
"""
Surprisingly, seaborn was not a good choice for plotting stacked barplot, so we use pandas
"""
import pandas as pd
import matplotlib.pyplot as plt
import sys

# Create the list for the number for each SV type at each depth
del_list = [10735, 12642, 13638, 14246, 14815, 15411]
dup_list = [1783, 3228, 4056, 4418, 4684, 4874]
ins_list = [16806, 19956, 21467, 22333, 23258, 24263]
inv_list = [634, 1058, 1514, 1885, 2388, 3004]
tra_list = [5316, 10148, 15570, 19796, 25426, 32066]

index_list = ["5X", "10X", "15X", "20X", "25X", "max"]

# Create pandas dataframe
df = pd.DataFrame({"Deletion" : del_list, "Duplication" : dup_list, "Insertion" : ins_list, "Inversion" : inv_list, "Translocation" : tra_list}, index = index_list)

# Create stacked barplot, rot=0, means do not rotate the x axis
df.plot.bar(stacked=True, rot=0, ylim=(0, 84000))

plt.xlabel('Sequencing depth', fontsize=18)

plt.ylabel('Number of variations', fontsize=18)

plt.title('Number of SVs in 22Rv1 (PacBio) across depths', fontsize=22)

plt.xticks(size = 16)

plt.yticks(size = 16)

plt.legend(fontsize=18)
# Save the plot as a PNG file with the desired size
# plt.gcf use the plt from matplotlib.pyplot to get the picture self and then resize  
plt.gcf().set_size_inches(12, 8)

# Save the plot as a PNG file
plt.savefig('22Rv1_SV_type_across_depth.png', dpi=900)