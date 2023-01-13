#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import os

suffix = sys.argv[1]
delete_list = []

for root, dirs, files in os.walk('.'):
    for file in files:
        if file.endswith(suffix):
            path = os.path.join(root, file)
            delete_list.append(path)
            
for _ in delete_list:
    answer = input(f"Delete the file in {_}? [y/n] ").lower()
    if answer == 'y' or answer == 'yes':
        os.remove(_)
        print(f"{path} deleted successfully!")
    else:
        sys.exit()