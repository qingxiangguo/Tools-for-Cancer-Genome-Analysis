# The installation and usage of Qualimap

## 1. About

Bam file QC.

## 2. Installation

```
mamba install -c bioconda qualimap
```

## 3. Usage

```bash
qualimap bamqc -bam /home/qgn1237/qgn1237/4_single_cell_SV_chimera/20230624_1st_PC3_bulk_ONT_DNA_MK1c/data_QC_assessment_comparison/PC3_Pacbio_hifi_DNA/PC3_hifi_DNA_3X_downsampled/PC3_3X.bam -c -outdir /home/qgn1237/qgn1237/4_single_cell_SV_chimera/20230624_1st_PC3_bulk_ONT_DNA_MK1c/data_QC_assessment_comparison/PC3_Pacbio_hifi_DNA/PC3_hifi_DNA_3X_downsampled/qualimap --java-mem-size=70G
```