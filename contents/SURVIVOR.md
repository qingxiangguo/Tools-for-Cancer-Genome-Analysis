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
rl ../../pbsv/pbsv.var.vcf
 ../../cutesv/cutesv.vcf ../../svim/variants.vcf ../../debreak/debreak.vcf ../../svdss/svs_poa.vcf ../../sniffles/sniffles.vcf > list_vcf

SURVIVOR merge list_vcf 10 1 0 0 0 50 merged_filtered.vcf 
```

Max distance between breakpoints: 10
Minimum number of supporting caller: 1
Take the type into account (1==yes, else no): 0
Take the strands of SVs into account (1==yes, else no): 0
Estimate distance based on the size of SV (1==yes, else no): 0
Minimum size of SVs to be taken into account.

Output files: merged_filtered.vcf

### 2.2 Get SV statistics

```bash
SURVIVOR stats merged_filtered.vcf -1 -1 -1 vcf_summary > merged_filtered.vcf.stat
```

This will output: merged_filtered.vcf.stat  vcf_summary_CHR vcf_summary vcf_summarysupport

### 2.3 Get overlap between software or samples

```bash
perl -ne 'print "$1\n" if /SUPP_VEC=([^,;]+)/' merged_filtered.vcf | sed -e 's/\(.\)/\1 /g' > merged_filtered_overlapp.txt
```

```R
#here is the R code
library(VennDiagram)

t = read.table("merged_filtered_overlapp.txt", header = FALSE)

dst = data.matrix(t(t))

venn.diagram(list(PBSV = which(t[,1]==1), SVIM = which(t[,2]==1)) , fill = c("orange" ,"blue"), alpha = c(0.5, 0.5), cex = 2, lty = 2, filename = "my_software_overlapp.png");

#exit R 
quit()
```

### 2.3 Benchmarking against ground truth

```bash
SURVIVOR merge list_vcf 500 1 0 0 0 0 merged_filtered.vcf

grep -v '^#' merged_filtered.vcf | grep 'SUPP_VEC=11' | wc -l > TP
grep -v '^#' merged_filtered.vcf | grep 'SUPP_VEC=01' | wc -l > FP
grep -v '^#' merged_filtered.vcf | grep 'SUPP_VEC=10' | wc -l > FN

```