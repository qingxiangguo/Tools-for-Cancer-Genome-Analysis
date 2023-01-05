# The installation and usage of BCFtools

## 1. About

 Utilities for variant calling and manipulating VCFs and BCFs

## 2. Installation and Usage

### mamba

```bash
mamba install -c bioconda bcftools=1.15
```

### 2.1 filter INFO field, like, filter the results of SVs detected by delly 

```bash
bcftools filter -i 'FILTER="PASS" && INFO/PE>=1' 22Rv1.5X.vcf > test.vcf
# This will include the variants that meet the conditions at the same time
```