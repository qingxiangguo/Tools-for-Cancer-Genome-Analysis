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

## Management
### [Conda and Mamba](/contents/conda.md)

## DNA and RNA-seq aligner (splice aware)
### [Minimap2](/contents/minimap2.md)

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

## Genome structural variation analysis
### [PBSV](/contents/pbsv.md)

## Gene fusion analysis - RNA-seq level
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
sra_accession_number = ["SRR11951439", "SRR11951443", "SRR11951444", "SRR11951445", "SRR11951446", "SRR11951447", "SRR11951448",
                        "SRR11951449", "SRR11951450", "SRR11951451", "SRR11951452", "SRR11951453", "SRR11951454", "SRR11951455",
                        "SRR11951456", "SRR11951457", "SRR11951458", "SRR11951459", "SRR11951460", "SRR11951461", "SRR11951462",
                        "SRR11951463", "SRR11951464", "SRR11951465", "SRR11951466", "SRR11951467", "SRR11951468", "SRR11951469",
                        "SRR11951470", "SRR11951471", "SRR11951472", "SRR11951473", "SRR11951474", "SRR11951475", "SRR11951476",
                        "SRR11951477", "SRR11951478", "SRR11951479", "SRR11951480", "SRR11951481", "SRR11951482", "SRR11951483",
                        "SRR11951484", "SRR11951485", "SRR11951486", "SRR11951487", "SRR11951488", "SRR11951489", "SRR11951490",
                        "SRR11951491", "SRR11951492", "SRR11951493", "SRR11951494", "SRR11951495", "SRR11951496", "SRR11951497",
                        "SRR11951498", "SRR11951499", "SRR11951500", "SRR11951501", "SRR11951502", "SRR11951503", "SRR11951504",
                        "SRR11951505", "SRR11951506", "SRR11951507", "SRR11951508", "SRR11951509", "SRR11951510", "SRR11951511",
                        "SRR11951512", "SRR11951513", "SRR11951514", "SRR11951515", "SRR11951516", "SRR11951517", "SRR11951519",
                        "SRR11951520", "SRR11951521", "SRR11951522", "SRR11951523", "SRR11951524", "SRR11951525", "SRR11951526",
                        "SRR11951527", "SRR11951528", "SRR11951530", "SRR11951531", "SRR11951532", "SRR11951533", "SRR11951534",
                        "SRR11951535", "SRR11951536", "SRR11951537", "SRR11951542"]

# This will download the .sra files to ~/ncbi/public/sra/ (will create directory if not present)
for sra_id in sra_accession_number:
    print("Current downloading" + sra_id)
    prefetch_cmd = "prefetch " + sra_id + " -O /home/qgn1237/qgn1237/2_raw_data/smooth_seq_95_sc_K562_SMRT" # Note the space between commands
    print("The running command is " + prefetch_cmd)
    subprocess.call(prefetch_cmd, shell=True)

# This will create the fastq file from the downloaded sra file
for sra_id in sra_accession_number:
    print("Creating fastq file for" + sra_id)
    fasterq_dump_cmd = "fasterq-dump /home/qgn1237/qgn1237/2_raw_data/smooth_seq_95_sc_K562_SMRT/" + sra_id + " -O /home/qgn1237/qgn1237/2_raw_data/smooth_seq_95_sc_K562_SMRT/" + sra_id # Use the fstring here
    print("The running command is " + fasterq_dump_cmd)
    subprocess.call(fasterq_dump_cmd, shell=True)
```

## Batch remove a certain type of files from a directory
```
python batch_delete_all_sra_files.py
```

## Use Mamba or Conda environment in Slurm Batch sbumit in NU Quest
You have to resource the mamba by yourself
```
#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics ## Required: (buyin, short, normal, long, gengpu, genhimem, etc)
#SBATCH --time=03:00:00 ## Required: How long will the job need to run (remember different partitions have restrictions on this parameter)
#SBATCH --nodes=1 ## how many computers/nodes do you need (no default)
#SBATCH --ntasks-per-node=6 ## how many cpus or processors do you need on per computer/node (default value 1)
#SBATCH --mem=35G ## how much RAM do you need per computer/node (this affects your FairShare score so be careful to not ask for more than you need))
#SBATCH --job-name=allen_download ## When you run squeue -u  this is how you can identify the job
#SBATCH --output=output.log ## standard out and standard error goes to this file

# A regular comment in Bash
/home/qgn1237/2_software/mambaforge/bin/mamba init

source ~/.bashrc

mamba activate mamba666

minimap2 -d /home/qgn1237/qgn1237/1_my_database/GRCh38_p13/minimap2_index/GRCh38.p13.genome.mmi /projects/b1171/qgn1237/1_my_database/GRCh38_p13/GRCh38.p13.genome.fa -t 6
```

## Batch create directories from list
```
python ./batch_create_directory_from_list.py
```

## Cancel a job in NU Quest
```scancel -u NETID ```
