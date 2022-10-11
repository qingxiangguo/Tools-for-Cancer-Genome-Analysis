# Tools-and-tricks-for-Cancer-Genome-Analysis
Installation and usage for various tools for cancer genomics

# Contributors
Qingxiang (Allen) Guo  
Postdoctoral Fellow  
Northwestern University, Feinberg School of Medicine
qingxiang.guo@northwestern.edu

# Introduction
In this section, I provide the installation and usage for a wide range of bioinformatics tools, especially for cancer genomics. This repo will be kept updating. Feedback or experience is warmly welcomed.

# Tools
## Data transfer
### [fasterq-dump](/contents/fasterq.md)

## Splice unware aligner
### [BWA-MEM](/contents/bwa.md)
### [BWA-MEM2](/contents/bwa2.md)

## RNA-seq aligner (splice aware)
### [STAR](/contents/STAR.md)

##  Manipulating alignments
### [Samtools](/contents/samtools.md)
### [Picard](/contents/picard.md)

## Indel calling
### [transindel](/contents/transindel.md)

## Gene fusion analysis
### [Arriba](/contents/arriba.md)

# Tricks
## Find and load R in Northwestern quest  
You can see which versions of R are available on Quest, and which version is the default, with the command  
```
module spider R
```

You can make a particular version of R available to use by typing the full module name with the version included as listed in the output  
```
module load R/4.2.0
```
