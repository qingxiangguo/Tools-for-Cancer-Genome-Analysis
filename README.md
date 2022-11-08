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

# Other tips
## Find and load R in Northwestern quest  
You can see which versions of R are available on Quest, and which version is the default, with the command  
```
module spider R
```

You can make a particular version of R available to use by typing the full module name with the version included as listed in the output  
```
module load R/4.2.0
```

## Batch download SRA files from NCBI
Go to the NCBI SRA site, select the SRA file you need and go to "Run selector" to got batch list in this way.

Use perl one-liner to process the header

```
perl -p -i -e 's/(\S+)\n/"$1", /g' list
```

Revise the script and batch download with the python script:

```
# _*_ coding=utf-8 _*_
import subprocess

# Insert your SRR accession numbers here
sra_accession_number = ["SRR11951439", "SRR11951443", "SRR11951444", "SRR11951445", "SRR11951446", "SRR11951447", "SRR11951448"]

# This will download the .sra files to your custom file (will create directory if not present)
for sra_id in sra_accession_number:
    print("Current downloading" + sra_id)
    prefetch_cmd = "prefetch " + sra_id + " -O /home/qgn1237/qgn1237/2_raw_data/smooth_seq_95_sc_K562_SMRT" # Note the space between commands
    print("The running command is " + prefetch_cmd)
    subprocess.call(prefetch_cmd, shell=True)

# This will create the fastq file from the downloaded sra file
for sra_id in sra_accession_number:
    print("Creating fastq file for" + sra_id)
    fasterq_dump_cmd = "fasterq-dump /home/qgn1237/qgn1237/2_raw_data/smooth_seq_95_sc_K562_SMRT/" + {sra_id} + " -O /home/qgn1237/qgn1237/2_raw_data/smooth_seq_95_sc_K562_SMRT" # Use the fstring here
    print("The running command is " + fasterq_dump_cmd)
    subprocess.call(fasterq_dump_cmd, shell=True)
```
