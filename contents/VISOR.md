# The installation and usage of VISOR

## 1. About

VISOR (opens new window)is an efficient and versatile command-line application capable to simulate structural variants and small/single-nucleotide variants in a haplotype-resolved manner. 

## 2. Installation and Usage

bseq.c:1:10: fatal error: zlib.h: No such file or directory

```bash

git clone https://github.com/davidebolo1993/VISOR.git
cd VISOR
pip install -r requirements.txt
python setup.py install
VISOR --help

# You need to switch to python 3.9 and let the script install the right version of minimap2 and samtools for you
```

### 2.1 Generate simulated genome with Hack mode (random SVs)

The input SVs should be in BED FORMAT

chr1	10000000	11000000	deletion	None	5
chr1	20000000	21000000	inversion	None	0
chr1	30000000	31000000	tandem duplication	2	10
chr1	40000000	41000000	inverted tandem duplication	2	10
chr1	50999999	51000000	insertion	TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT	2
chr1	80000000	81000000	translocation cut-paste	h1:chr1:180000000:forward	5  # Means ABC to BAC or ACB in a single chrom
chr1	90000000	91000000	translocation copy-paste	h1:chr1:190000000:reverse	0  # Means interspersed duplication
chr1	100000000	101000000	reciprocal translocation	h2:chr2:200000000:forward:forward	3

VISOR can convert real-world VCF into BED, or you can randomly generate by your self.

I want to simulate SVs considering both haplotypes.

```bash
module load R/4.3.0

Add a line in randomregion.r, it will install in your local library

.libPaths("~/R/library")

#this script requires a R version >= 3.5.0. It will try to install the needed packages if not already in your package list

cut -f1,2 ~/qgn1237/1_my_database/GRCh38_p13/GRCh38.p13.genome.fa.fai > chrom.dim.tsv # From index

# Need bedtools to be installed and in path

Rscript randomregion.r -d /projects/b1171/qgn1237/6_SV_VCF_merger/VISOR_simulation/chrom.dim.tsv -n 18000 -l 200000 -v 'deletion,insertion,tandem duplication,inverted tandem duplication,translocation copy-paste,translocation cut-paste,reciprocal translocation,inversion' -r '35:35:5:5:5:5:5:5' | sortBed > HACk.random.bed

# Generate 18000 non-overlapping random variants, with mean length 200 Kb, with a fixed ratio.

```


