# The installation and usage of manta

## 1. About

 Manta calls structural variants (SVs) and indels from mapped paired-end sequencing reads.

## 2. Installation and Usage

### mamba

```bash
mamba install -c bioconda delly
```

### 2.1 run the manta config script

```bash
# The bam file has to be indexed first
delly call -t INS -g /projects/b1171/twp7981/database/genome/hg38.fa /projects/b1171/qgn1237/4_single_cell_SV_chimera/2_test_5_prostate_depth_SV_impact/NGS_data/22Rv1/max_depth/22Rv1.max.bam > /projects/b1171/qgn1237/4_single_cell_SV_chimera/2_test_5_prostate_depth_SV_impact/NGS_data/22Rv1/max_depth/delly/22Rv1_delly_INS.vcf &&

bcftools filter -i 'FILTER="PASS"' /projects/b1171/qgn1237/4_single_cell_SV_chimera/2_test_5_prostate_depth_SV_impact/NGS_data/22Rv1/max_depth/delly/*_delly_INS.vcf > /projects/b1171/qgn1237/4_single_cell_SV_chimera/2_test_5_prostate_depth_SV_impact/NGS_data/22Rv1/max_depth/delly/22Rv1_INS_filtered.vcf
```

Output files: 22Rv1.5X.vcf