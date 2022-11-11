# The installation and usage of Minimap2
# 1. About
the use-case, dataset and the running machine.
# 2. Installation and Usage
## mamba
```
mamba install -c bioconda minimap2
```

## 2.1 Build genome index, so you don't need to build this everytime
The mmi file is the index
```
minimap2 -d /home/qgn1237/qgn1237/1_my_database/GRCh38_p13/minimap2_index/GRCh38.p13.genome.mmi /projects/b1171/qgn1237/1_my_database/GRCh38_p13/GRCh38.p13.genome.fa -t 6
```

## 2.2 Running mapping job with BWA-MEM

## 2.3 Output files

# 3. Citation
