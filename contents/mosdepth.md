# The installation and usage of mosdepth
# 1. About
fast BAM/CRAM depth calculation for WGS, exome, or targeted sequencing.

# 2. Installation and Usage
## mamba
```
mamba install -c bioconda mosdepth
```

## 2.1 Depth for WGS

```
mosdepth -n --fast-mode SRR11951439 SRR11951439_sort.bam
# -n means dont output per-base depth. skipping this output will speed execution
```

## 2.2 Output
The useful output is SRR11951439.mosdepth.summary.txt






