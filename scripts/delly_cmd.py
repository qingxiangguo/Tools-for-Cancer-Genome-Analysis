#!/home/qgn1237/2_software/mambaforge/envs/mamba666/bin/python
# -*- coding: utf-8 -*-

import sys
import os

# Check if the correct number of arguments are provided
if len(sys.argv) != 2:
    print("Error: Incorrect number of arguments provided.")
    print("Usage: generate_cmd.py <input_file>")
    sys.exit()

# Get the absolute path and basename of the input file
input_file = sys.argv[1]
input_file_abs_path = os.path.abspath(input_file)
input_file_dir_name = os.path.dirname(input_file_abs_path)
input_file_base_name = os.path.basename(input_file_abs_path).split(".")[0]
print(input_file_dir_name , "dirname")
# Generate the command using the template and the arguments
cmd1 = f"delly call -g /projects/b1171/twp7981/database/genome/hg38.fa {input_file_abs_path} > {input_file_dir_name}/delly/{input_file_base_name}_delly.vcf"
cmd2 = f"bcftools filter -i 'FILTER=\"PASS\" && INFO/PE>=1' {input_file_dir_name}/delly/*_delly.vcf > {input_file_dir_name}/delly/{input_file_base_name}_filtered.vcf"

# Write the command to the cmd_list file in append mode
with open("cmd_list", "w") as f:
    f.write("#!/bin/bash\n"
            "#SBATCH --account=b1042\n"
            "#SBATCH --partition=genomics ## Required: (buyin, short, normal, long, gengpu, genhimem, etc)\n"
            "#SBATCH --time=01:00:00 ## Required: How long will the job need to run (remember different partitions have restrictions on this parameter)\n"
            "#SBATCH --nodes=1 ## how many computers/nodes do you need (no default)\n"
            "#SBATCH --ntasks-per-node=1 ## how many cpus or processors do you need on per computer/node (default value 1)\n"
            "#SBATCH --mem=30G ## how much RAM do you need per computer/node (this affects your FairShare score so be careful to not ask for more than you need))\n"
            "#SBATCH --job-name=allen_sam1 ## When you run squeue -u  this is how you can identify the job\n"
            "#SBATCH --output=output.log ## standard out and standard error goes to this file\n\n"
            "# A regular comment in Bash\n"
            "/home/qgn1237/2_software/mambaforge/bin/mamba init\n"
            "source ~/.bashrc\n"
            "mamba activate mamba666\n"
            )
    f.write(cmd1 + " &&" "\n\n")
    f.write(cmd2 + "\n")

print("Command generated and written to cmd_list file.")