#!/usr/bin/env python
# -*- coding: utf-8 -*-
import pandas as pd
import sys

def process_fai(fai_files):
    dfs = []
    for fai_file in fai_files:
        df = pd.read_csv(fai_file, sep="\t", header=None, usecols=[0,1], names=["chromosome", "length"])
        dfs.append(df)
        
    merged_df = pd.concat(dfs) # Merge dataframe
    max_df = merged_df.groupby("chromosome")["length"].max().reset_index() # Group by chromosome and get the max length in each group
    
    return max_df

fai_files = sys.argv[1:]  # Get the parameter
max_df = process_fai(fai_files)

max_df.to_csv("maxdims.tsv", sep="\t", header=False, index=False)
