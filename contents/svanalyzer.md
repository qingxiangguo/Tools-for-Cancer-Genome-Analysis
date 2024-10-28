# The installation and usage

## 1. About

SVbenchmark compares a set of “test” structural variants in VCF format to a known truth set (also in VCF format) and outputs estimates of sensitivity and specificity.

## 2. Installation and Usage

### mamba

```bash
conda create -n svanalyzer
conda activate svanalyzer
conda install -c bioconda svanalyzer
```

### 2.1 Benchmarking VCF file

```bash
svanalyzer benchmark --ref <reference FASTA file> --test <VCF-formatted file of variants to test> --truth <VCF-formatted file of true variants>
```