# _*_ coding=utf-8 _*_
import glob
import os

# Get the route of your target file
for f in glob.glob("/home/qgn1237/qgn1237/2_raw_data/smooth_seq_95_sc_K562_SMRT/*/*.sra"):
        os.remove(f)  # Remove the file according the route
