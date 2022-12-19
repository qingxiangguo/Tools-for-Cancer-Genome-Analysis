#!/home/qgn1237/2_software/mambaforge/envs/mamba666/bin/python
# _*_ coding=utf-8 _*_
import os
import sys

# Get the suffix from the command line argument
suffix = sys.argv[1]

# Walk through all the directories and subdirectories
for root, dirs, files in os.walk('.'):
    # Check each file in the current directory
    for file in files:
        # If the file has the specified suffix, print the full path
        if file.endswith(suffix):
            print(os.path.join(root, file))