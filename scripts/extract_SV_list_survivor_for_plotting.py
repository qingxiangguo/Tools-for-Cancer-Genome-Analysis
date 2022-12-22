#!/home/qgn1237/2_software/mambaforge/envs/mamba666/bin/python
# _*_ coding=utf-8 _*_
"""
The input is the ouput of survivor like:
Processing: 35274
Parsing done:
Tot     DEL     DUP     INS     INV     TRA
35274   10735   1783    16806   634     5316
Processing: 47032
Parsing done:
Tot     DEL     DUP     INS     INV     TRA
47032   12642   3228    19956   1058    10148
Processing: 56245
Parsing done:
Tot     DEL     DUP     INS     INV     TRA
56245   13638   4056    21467   1514    15570
Processing: 62678
Parsing done:
Tot     DEL     DUP     INS     INV     TRA
62678   14246   4418    22333   1885    19796
"""
import sys

del_list = []
dup_list = []
ins_list = []
inv_list = []
tra_list = []

# Read the input text file
with open(sys.argv[1], 'r') as f:
    for line in f:
        line.strip()
        if line.startswith('Tot'):
            fields = next(f).strip().split('\t') # The next line should be strip()
            del_list.append(fields[1])
            dup_list.append(fields[2])
            ins_list.append(fields[3])
            inv_list.append(fields[4])
            tra_list.append(fields[5])

print("del_list = ", "[", ', '.join(str(x) for x in del_list), "]")
print("dup_list = ", "[", ', '.join(str(x) for x in dup_list), "]")
print("ins_list = ", "[", ', '.join(str(x) for x in ins_list), "]")
print("inv_list = ", "[", ', '.join(str(x) for x in inv_list), "]")
print("tra_list = ", "[", ', '.join(str(x) for x in tra_list), "]")

