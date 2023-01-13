# The installation and usage of SVIM

## 1. About

SVIM is able to detect and classify the following six classes of structural variation: deletions, insertions, inversions, tandem duplications, interspersed duplications and translocations

## 2. Installation and Usage

### mamba

```bash
mamba install -c bioconda svim
```

### 2.1 Assume that you've already index and sort the BAM file first

```bash
svim alignment /home/qgn1237/qgn1237/4_single_cell_SV_chimera/1_smooth_seq_95_sc_K562_SMRT/SRR11951439/svim /projects/b1171/qgn1237/4_single_cell_SV_chimera/1_smooth_seq_95_sc_K562_SMRT/SRR11951439/SRR11951439_sort.bam ~/qgn1237/1_my_database/GRCh38_p13/GRCh38.p13.genome.fa
```

### 2.2 Output files

SRR11951439_sort.var.vcf

### 2.3 Filtering

Since the output SVIM is unfiltered, we have to filter them manually
This is very very important since SVIM ouput almost everything.

```bash
for dir in *depth/; do cd "$dir"; cd svim; filter_vcf_based_on_quality.py variants.vcf 2 > filtered_variant.vcf; cd ../..; done
```

Or you can do it with BCFtools

```bash
bcftools view -i 'QUAL >= 10' variants.vcf'.
# or you can do
filter_vcf_based_on_quality.py variants.vcf 2 > filtered_variant.vcf
```

For high-coverage datasets (>40x), we would recommend a threshold of 10-15. For low-coverage datasets, the threshold should be lower (>3-5). For 30 I choose 8.

Or you can do not do this

```bash
for dir in *depth/; do cd "$dir"; cd svim; filter_vcf_based_on_quality.py variants.vcf 2 > filtered_variant.vcf; cd ../..; done
# For 30x, the value is 7
# For 25x, the value is 6
# for 20x, the value is 5
# for 15x, the value is 4
# for 10x, the value is 3
# for 5x, the value is 2
```