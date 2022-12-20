# The installation and usage of SURVIVOR

## 1. About

 Simulating/evaluating SVs, merging and comparing SVs within and among samples, and includes various methods to reformat or summarize SVs

## 2. Installation and Usage

### mamba

```bash
mamba install -c bioconda survivor=1.0.7
```

### 2.1 Merge VCF file

```bash
# Get the absolute path of VCF files to merge
readline -f SV*.vcf > list_vcf

SURVIVOR merge list_vcf 500 1 1 0 0 50 merged_filtered.vcf 
```

Here the filtering parameters of SURVIVOR are the following:
Max distance between breakpoints: 500
Minimum number of supporting caller: 1
Take the type into account (1==yes, else no): 1
Take the strands of SVs into account (1==yes, else no): 0
Estimate distance based on the size of SV (1==yes, else no): 0

### 2.3 Output files

merged_filtered.vcf
