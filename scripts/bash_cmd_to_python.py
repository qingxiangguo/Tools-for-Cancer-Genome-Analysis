#!/home/qgn1237/2_software/mambaforge/envs/mamba666/bin/python
# _*_ coding=utf-8 _*_
import glob
import subprocess

for path in glob.glob("/home/qgn1237/working/*/*/SURVIVOR"):
    cmd_1 = "cd " + path + "; " \
        + "get_plot_VCF_length_distribution_histo.py -i merged_filtered.vcf -v DEL;" \
        + "get_plot_VCF_length_distribution_histo.py -i merged_filtered.vcf -v INS;" \
        + "get_plot_VCF_length_distribution_histo.py -i merged_filtered.vcf -v INV;" \
        + "get_plot_VCF_length_distribution_histo.py -i merged_filtered.vcf -v DUP"
    subprocess.call(cmd_1, shell=True)




