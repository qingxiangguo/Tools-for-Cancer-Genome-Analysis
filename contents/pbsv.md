# The installation and usage of PBSV
# 1. About
pbsv calls insertions, deletions, inversions, duplications, and translocations. Both single-sample calling and joint (multi-sample) calling are provided. pbsv is most effective for:

insertions 20 bp to 10 kb  
deletions 20 bp to 100 kb  
inversions 200 bp to 10 kb  
duplications 20 bp to 10 kb  
translocations between different chromosomes or further than 100kb apart on a single chromosome  

# 2. Installation and Usage
## mamba
```
mamba install -c bioconda PBSV=2.8.0
```

## 2.1 Discover signatures
It is highly recommended to provide one tandem repeat annotation .bed file of your reference to pbsv discover via --tandem-repeats. This increases sensitivity and recall. 

```
pbsv discover --tandem-repeats ~/qgn1237/1_my_database/GRCh38_p13/tandem_repeats/human_GRCh38_no_alt_analysis_set.trf.bed ./SRR11951439_sort.bam ./SRR11951439_sort.svsig.gz
```

## 2.2 Call structural variants and assign genotypes
```

```
-a means outputing in SAM format; -x means a preset depends on your task aim.

## 2.3 Output files



