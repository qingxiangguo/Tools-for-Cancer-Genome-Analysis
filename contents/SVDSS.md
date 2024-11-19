# The installation and usage of SVDSS

## 1. About
A SV caller for long-read single-molecular sequencing data.

## 2. Usage

# You can generate the slurm.sh by SVDSS_step1_generator.py

./SVDSS_step1_generator.py --bam ../minimap2/visor_ONT.bam \
    --reference /path/to/reference.fa \
    --index /path/to/index.fmd \
    --time 03:40:00
```


