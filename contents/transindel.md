# The installation and usage of transindel
# 1. About
transIndel is used to detect indels (insertions and deletions) from DNA-seq or RNA-seq data by parsing chimiric alignments from BWA-MEM.
# 2. Installation and Usage
## Prerequisites
Samtools/1.0 or newer (http://www.htslib.org/) Python 3.6 or newer (https://www.python.org/)  
Python packages:  
Pysam/0.13.0 or newer (https://pypi.org/project/pysam)  
HTSeq/0.6.1 or newer (https://pypi.python.org/pypi/HTSeq)  

## Getting Soure Code
```
git clone git://github.com/cauyrd/transIndel.git
cd transIndel
```

## 2.1 Build new BAM file with redefined CIGAR string

For DNA-seq:

```
python transIndel_build_DNA.py -i input_bam_file -o output_bam_file [options]
```

For RNA-seq

```
python transIndel_build_RNA.py -i input_bam_file -r reference_genome_fasta -g gtf_file -o output_bam_file [options]
```

```
-h, --help            show this help message and exit
-i INPUT, --input INPUT
                    Input BAM file
-o OUTPUT, --output OUTPUT
                    Output BAM file
-r REF, --ref REF     reference genome used for analyzing RNA-seq data
-g GTF, --gtf GTF     gene annotatino file used for analyzing RNA-seq data
-s SPLICE_BIN, --splice_bin SPLICE_BIN
                    splice site half bin size (default: 20)
-m MAPQ, --mapq MAPQ  minimal MAPQ of read from BAM file for supporting
                    Indel (default: 15)
-l LENGTH, --length LENGTH
                     Maximum deletion length to be detected (default:1000000)
-v, --version         show program's version number and exit
```

<b>Output files</b>  
your_output_bam_file			:BAM file for CIGAR string redefinement.

transIndel generates the following optional fields in output BAMs  

OA  original representative alignment; format: (pos,CIGAR)  
JM  splicing junction reads; infered from GTF or splicing motif (used in RNA-seq BAM)  

## 2.2 Call indel
```
python transIndel_call.py -i input_bam_from_transIndel_build -o output_vcf_filename_prefix [options]	
```

```
-h, --help            show this help message and exit
-i INPUT, --input INPUT
                        Input BAM file
-o OUTPUT, --output OUTPUT
                    output VCF file prefix
-c AO, --ao AO        minimal observation count for Indel (default: 4)
-d DEPTH, --depth DEPTH
                 minimal depth to call Indel (default: 10)
-f VAF, --vaf VAF     minimal variant allele frequency (default: 0.1)
-l LENGTH, --length LENGTH
                 minimal Indel length (>=1) to report (default: 10)
-m MAPQ, --mapq MAPQ  minimal MAPQ of read from BAM file to call Indel
                 (default: 15)
-t REGION             Limit analysis to targets listed in the BED-format
                 FILE or a samtools region string
-v, --version         show program's version number and exit
```

<b>Output files</b>  
output_vcf_file   			:Reported Indels with VCF format

