#!/usr/bin/env python 
# -*- coding: utf-8 -*-
import sys

outputfile = open("cmd_list", "w")

with open(sys.argv[1], 'r') as f:
    for line in f:
        line.strip()
        outputfile.write(f"prefetch {line} -O /home/qgn1237/qgn1237  \
                         /4_single_cell_SV_chimera/3_feature_WGA_chimera \
                         _PC3_22Rv1_short_reads/PC3_NGS_MDA_MALBAC_WGS/PC3_bulk \
                         && fasterq-dump --split-files /home/qgn1237/qgn1237  \
                         /4_single_cell_SV_chimera/3_feature_WGA_chimera \
                         _PC3_22Rv1_short_reads/PC3_NGS_MDA_MALBAC_WGS/PC3_bulk/{line}")
        
outputfile.close()        