# The installation and usage of Delly

## 1. About

https://github.com/dellytools/delly

## 2. Installation and Usage

### mamba

```bash
mamba install -c bioconda delly
```

### 2.1 Usage 1

```bash
# The bam file has to be indexed first
delly call -g /projects/b1171/twp7981/database/genome/hg38.fa /projects/b1171/qgn1237/4_single_cell_SV_chimera/2_test_5_prostate_depth_SV_impact/NGS_data/22Rv1/max_depth/22Rv1.max.bam > /projects/b1171/qgn1237/4_single_cell_SV_chimera/2_test_5_prostate_depth_SV_impact/NGS_data/22Rv1/max_depth/delly/22Rv1_delly_INS.vcf &&
s
bcftools filter -i 'FILTER="PASS"' /projects/b1171/qgn1237/4_single_cell_SV_chimera/2_test_5_prostate_depth_SV_impact/NGS_data/22Rv1/max_depth/delly/*_delly_INS.vcf > /projects/b1171/qgn1237/4_single_cell_SV_chimera/2_test_5_prostate_depth_SV_impact/NGS_data/22Rv1/max_depth/delly/22Rv1_INS_filtered.vcf
```

Output files: 22Rv1.5X.vcf

### 2.2 Simple usage

```bash
delly call -g /projects/b1171/qgn1237/1_my_database/GRCh38_p13/GRCh38.p13.genome.fa input.bam > delly.vcf
```