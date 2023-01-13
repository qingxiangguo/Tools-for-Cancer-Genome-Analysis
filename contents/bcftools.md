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

### 2.2 Intersect two vcf

```bash
bcftools sort 10X_long.vcf > 10X_long.sort.vcf

bgzip -c 10X_long.vcf > 10X_long.vcf.gz

bcftools sort 10X_long.vcf.gz > 10X_long.sor.vcf.gz 

bcftools index 10X_long.sort.vcf.gz > 10X_long.sort.vcf.gz.csi

# You should to build indexes for them:

# Try to run isec:

bcftools isec isec.a.vcf.gz isec.b.vcf.gz -p dir
```

I don't like bcftools isec, it doesn't work well. It needs the two vcf files to be exactly the same.
I recommend surpyvor, a SURVIVOR wrapper.