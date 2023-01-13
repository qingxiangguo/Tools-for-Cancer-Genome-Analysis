#!/usr/bin/env python 
# -*- coding: utf-8 -*-
import sys
import random

random_number = random.randint(0, 1000)

outputfile = open(f"cmd_list_{random_number}", "w")

with open(sys.argv[1], 'r') as f:
    for line in f:
        line = line.strip()
        outputfile.write(f"prefetch {line} -O /home/qgn1237/qgn1237"
                         f"/4_single_cell_SV_chimera/3_feature_WGA_chimera"
                         f"_PC3_22Rv1_short_reads/PC3_NGS_MDA_MALBAC_WGS/PC3_MALBAC"
                         
                         f" && fasterq-dump --split-files /home/qgn1237/qgn1237"
                         f"/4_single_cell_SV_chimera/3_feature_WGA_chimera"
                         f"_PC3_22Rv1_short_reads/PC3_NGS_MDA_MALBAC_WGS/PC3_MALBAC/{line} -o"
                         f"/home/qgn1237/qgn1237/4_single_cell_SV_chimera/3_feature_WGA_chimera"
                         f"_PC3_22Rv1_short_reads/PC3_NGS_MDA_MALBAC_WGS/PC3_MALBAC/{line} \n")
        
outputfile.close()        