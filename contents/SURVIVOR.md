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