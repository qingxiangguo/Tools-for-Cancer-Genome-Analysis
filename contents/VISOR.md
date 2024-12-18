# The installation and usage of VISOR

## 1. About

VISOR (opens new window)is an efficient and versatile command-line application capable to simulate structural variants and small/single-nucleotide variants in a haplotype-resolved manner. 

## 2. Installation and Usage

bseq.c:1:10: fatal error: zlib.h: No such file or directory

```bash
mamba create -y -n visorenv python=3.8
mamba activate visorenv
mamba install -c bioconda samtools=1.11 minimap2 htslib=1.11 bedtools
pip install --upgrade cython
git clone https://github.com/davidebolo1993/VISOR.git
cd VISOR
pip install -r requirements.txt
python setup.py install

git clone https://github.com/davidebolo1993/VISOR.git
cd VISOR
pip install -r requirements.txt
python setup.py install
VISOR --help

# You need to switch to python 3.9 and let the script install the right version of minimap2 and samtools for you
```

### 2.1 Generate simulated genome with Hack mode (random SVs) in only one haplotype (so all SVs will be 1/0 or 0/1)

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

Add a line in randomregion.r, it will install in your local library, to solve the sudo problem

.libPaths("~/R/library")

# Install the packages
# You have to comment the line in .bashrc to change from fish to bash shell to use the module function.

Rscript randomregion.r -h1

#this script requires a R version >= 3.5.0. It will try to install the needed packages if not already in your package list

cut -f1,2 ~/qgn1237/1_my_database/GRCh38_p13/GRCh38.p13.genome.fa.fai > chrom.dim.tsv # From index

# Need bedtools to be installed and in path

Rscript randomregion.r -d /projects/b1171/qgn1237/6_SV_VCF_merger/VISOR_simulation/chrom.dim.tsv -n 18000 -l 10000 -v 'deletion,insertion,tandem duplication,inverted tandem duplication,translocation copy-paste,translocation cut-paste,reciprocal translocation,inversion' -r '35:35:5:5:5:5:5:5' | sortBed > HACk.random.bed

# Generate 18000 non-overlapping random variants, with mean length 200 Kb, with a fixed ratio.

# Sometimes you need to convert VCF to  BED file, the tips is
vcf to BED, insertion point-1, start-1, end point unchanged
For example

chr1    16725233        pbsv.BND.chr1:16725233-chr1:234776445   C       [chr1:234776445[C       .       PASS    SVTYPE=BND;CIPOS=-4,14;MATEID=pbsv.BND.chr1:234776445-chr1:16725233;MATEDIST=218051212  GT:AD:DP        0/1:2,4:6

BED:

chr1	234776444	248956422	translocation cut-paste	h1:chr1:16725232:reverse	5

chr1    234776445       pbsv.BND.chr1:234776445-chr1:16725233   G       [chr1:16725233[G        .       PASS    SVTYPE=BND;CIPOS=-4,8;MATEID=pbsv.BND.chr1:16725233-chr1:234776445;MATEDIST=218051212   GT:AD:DP        1/1:1,4:5

BED:

chr1	16725232   	234776445	inversion	None	5

```

However, VISOR is not 0-based, the BED file is not a true BED file! VISOR BED is VCF style!
So to convert VCF to VISOR BED-like file, the tips is
insertion point unchanged, start unchanged, end point unchanged, nothing needs to be changed!
This is really convenient, but you have to aware that VISOR BED is a fake BED.

### 2.2 Generate simulated genome with Hack mode (random SVs) in two haplotype (so all SVs will be 1/0 or 0/1 or 1/1)

The weird part that is #reciprocal translocations are placed by default on a haplotype different then the one specified with the -i parameter (default to 1 -that is, h1 in the final BED). Other translocations types are placed on the same haplotype.

So if you want to simulate SVs on both haplotypes, you need to do run randomregion.r two times:

```bash
Rscript randomregion.r -d /projects/b1171/qgn1237/6_SV_VCF_merger/VISOR_simulation/chrom.dim.tsv -i 1 -n 18000 -l 10000 -v 'deletion,insertion,tandem duplication,inverted tandem duplication,translocation copy-paste,translocation cut-paste,reciprocal translocation,inversion' -r '35:35:5:5:5:5:5:5' | sortBed > HACk.random.h1.bed

Rscript randomregion.r -d /projects/b1171/qgn1237/6_SV_VCF_merger/VISOR_simulation/chrom.dim.tsv -i 2 -n 18000 -l 10000 -v 'deletion,insertion,tandem duplication,inverted tandem duplication,translocation copy-paste,translocation cut-paste,reciprocal translocation,inversion' -r '35:35:5:5:5:5:5:5' | sortBed > HACk.random.h2.bed

# However, in the h1 file, all the Svs will be on h1, only reciprocal translocations is on h2
# But in h2 file, all the SVs will be on h2, only reciprocal translocations is on hc3806b6_0
# h3 is a mistake ,there is no h3
# GL000250.2      1927905 1937904 reciprocal translocation        h3:chr8:112551381:reverse:forward       9
# We need to correct the h3 manually to h1

sed 's/h3/h1/g' HACk.random.h2.bed > HACk.random.h2.corrected.bed

# Then I want to add some decoy SVs only to haplotype 1

cat special_decoy.bed >> HACk.random.h1.bed
# Add some decoy

sortBed -g ~/qgn1237/1_my_database/GRCh38_p13/GRCh38.p13.genome.fa.fai -i HACk.random.h1.bed  > HACk.random.h1.sorted.bed

```

### 2.3 HACk generates a FASTA haplotype

```bash
# You need to switch to python 3.9, the environment you properly install VISOR
VISOR HAck -b HACk.random.h1.sorted.bed HACk.random.h2.corrected.bed -g ~/qgn1237/1_my_database/GRCh38_p13/GRCh38.p13.genome.fa -o simulated_genome

# I got some error:
bedtools: /lib64/libstdc++.so.6: version `GLIBCXX_3.4.20' not found (required by bedtools)
bedtools: /lib64/libstdc++.so.6: version `CXXABI_1.3.8' not found (required by bedtools)
bedtools: /lib64/libstdc++.so.6: version `GLIBCXX_3.4.21' not found (required by bedtools)
bedtools: /lib64/libstdc++.so.6: version `GLIBCXX_3.4.20' not found (required by bedtools)
bedtools: /lib64/libstdc++.so.6: version `CXXABI_1.3.8' not found (required by bedtools)

# The gcc is not compatible with bedtools installed in main environment path
# Install in current environment by mamba again

mamaba install -c bioconda bedtools

# Run again, I got another error
[17/05/2023 11:02:08][Error] Line 943: column 2 (chromosome start) contains an invalid coordinate (lower than chromosome start or greater than chromosome end)

# You should check the decoy coordinate typo, and aware that the VISOR BED is not a true BED format
# Finally work.
```

Now I have a folder:

h1.fa  h1.fa.fai  h2.fa  h2.fa.fai

Next you need to create a BED file containing the max length of each chromo between h1, h2, and ref genome.

If you want to simulate sequencing data generated from the entire reference genome (not just a particular chromosome), you need to first create a BED file describing the region you want to simulate. Since structural variants, such as insertions and deletions, may cause changes in the length of the chromosome, we need to consider all possible lengths. The maximum value is taken to ensure that the sequencing data we simulate will cover all possible regions and that a region of variation will not be missed due to improper length selection.

```bash
# To convert a VCF to Visor BED file
python vcf2visor_bed_input.py -i nstd106.GRCh38.variant_true_call.vcf -o visor_bed --homozygous_ratio 0.3
# Then compare the output results with original VCF file
python compare_vcf_bed.py
# you will get  visor_bed_haplotype1.bed and  visor_bed_haplotype2.bed
# If the chromosome format of input VCF is different from Genome fasta, you will have problem
# Then you need to change them manually

cut -f1 visor_bed_haplotype1.bed | sort | uniq > bed_chroms_all.txt
cut -f1 visor_bed_haplotype2.bed | sort | uniq >> bed_chroms_all.txt
grep "translocation" visor_bed_haplotype1.bed | cut -f5 | sed 's/h[12]:chr//; s/:.*//' >> bed_chroms_all.txt
grep "translocation" visor_bed_haplotype2.bed | cut -f5 | sed 's/h[12]:chr//; s/:.*//' >> bed_chroms_all.txt


sort genome_chroms.txt > genome_chroms_sorted.txt
sort bed_chroms_all.txt > bed_chroms_sorted.txt

comm -23 bed_chroms_sorted.txt genome_chroms_sorted.txt > diff_chroms.txt

# Delete those bad names
awk 'NR==FNR{a[$1]=1;next} {found=0; for(chr in a) if(index($0,chr)>0){found=1;break}; if(!found)print}' diff_chroms.txt visor_bed_haplotype2.bed > visor_bed_haplotype2.clean.bed
```

```bash
# To make this BED file, run:
VISOR_LASeR_BED_generator.py simulated_genome/h1.fa.fai simulated_genome/h2.fa.fai ~/qgn1237/1_my_database/GRCh38_p13/GRCh38.p13.genome.fa.fai

# You will get maxdims.tsv, this tells which part of genome the virtual sequencing generated.
# create a BED to simulate reads from whole genome, without coverage fluctuations (that is, capture bias value in 4th column is 100.0) and without normal contamination (that is, purity value in 5th column is 100.0)
awk 'OFS=FS="\t"''{print $1, "1", $2, "100.0", "100.0"}' maxdims.tsv > shorts.laser.simple.bed
```

### Run simulation and mapping

```bash
# Let's simulated fastq reads, so you can control the mapping steps
# You have to make sure that minimap2, samtools is installed in this environment

Then you may get error of 

import pysam
  File "/home/qgn1237/2_software/mambaforge/envs/mamba_py_39/lib/python3.9/site-packages/pysam/__init__.py", line 4, in <module>
    from pysam.libchtslib import *
ImportError: libcrypto.so.3: cannot open shared object file: No such file or directory

mamba install -c bioconda samtools openssl=3.0
mamba config --set channel_priority flexible

#!/bin/bash -l
#SBATCH --account=b1171
#SBATCH --partition=b1171
# set the number of nodes
#SBATCH --nodes=1
#SBATCH --ntasks=24
##SBATCH --exclusive
#SBATCH --mem=100G
# set max wallclock time
#SBATCH --time=334:00:00
# mail alert at start, end and abortion of execution
#SBATCH --mail-type=ALL

# run the application
cd $SLURM_SUBMIT_DIR
/home/qgn1237/2_software/mambaforge/bin/mamba init
source ~/.bashrc
mamba activate visorenv

VISOR LASeR -g ~/qgn1237/1_my_database/GRCh38_p13/GRCh38.p13.genome.fa -s ../../SV_bed_for_all/simulated_genome -b ../../SV_bed_for_all/100VAF.laser.simple.bed -o pacbio_fastq --threads 24 --coverage 5 --fastq --read_type pacbio

# The result is: r.fq  sim.srt.bam  sim.srt.bam.bai
```

# Generating Specific VAF SVs with VISOR

## Introduction

This tutorial guides you through using VISOR to simulate Structural Variants (SVs) with a specific Variant Allele Frequency (VAF) in a human diploid genome.

## Steps

1. **Create Variant BED File for a Single Haplotype**:
   Generate a BED file with SVs for one haplotype using the `randomregion.r` script. Ensure the SVs are the same for both haplotypes to simulate a 1/1 genotype.

   ```bash
   Rscript randomregion.r -d /projects/b1171/qgn1237/6_SV_VCF_merger/VISOR_simulation/chrom.dim.tsv -n 18000 -l 10000 -v 'deletion,insertion,tandem duplication,reciprocal translocation,inversion' -r '40:40:10:5:5' | sortBed > HACk.random.h1.bed
   cp HACk.random.h1.bed HACk.random.h2.bed
   ```

2. **Generate FASTA Haplotypes**:
   Run VISOR HAck to create two FASTA haplotypes with all 1/1 genotype SVs:

   ```bash
   VISOR HAck -b ../../SV_bed_for_all/HACk.random.h1.bed ../../SV_bed_for_all/HACk.random.h2.bed -g ~/qgn1237/1_my_database/GRCh38_p13/GRCh38.p13.genome.fa -o simulated_genome
   ```

3. **Prepare BED for Sequencing Simulation**:
   Create a BED file for SHORtS/LASeR with the desired purity (e.g., 80% for 80% VAF):

   ```bash
   VISOR_LASeR_BED_generator.py simulated_genome/h1.fa.fai simulated_genome/h2.fa.fai ~/qgn1237/1_my_database/GRCh38_p13/GRCh38.p13.genome.fa.fai
   awk 'OFS=FS="\t"''{print $1, "1", $2, "100.0", "100.0"}' maxdims.tsv > 100VAF.laser.simple.bed
   ```

4. **Run VISOR LASeR for Sequencing Simulation**:
   Simulate sequencing reads using LASeR:

   ```bash
   VISOR LASeR -g genome.fa -s simulated_genome -b 100VAF.laser.simple.bed -o ont_fastq --threads 12 --coverage 10 --fastq --read_type nanopore
   ```

## Conclusion

By following these steps, you can simulate SVs with a specific VAF in a human diploid genome using VISOR, providing a valuable tool for genetic research and analysis.


# Resuming an Interrupted VISOR LASeR Run

If your VISOR LASeR simulation was interrupted, you can manually resume it by identifying unfinished regions and continuing the simulation for those specific regions. Below are the steps and commands to help you achieve this.

## Steps to Resume Simulation

### 1. Identify Unfinished Regions

#### 1.1 Calculate Region Lengths

Extract information from your original BED file (`100VAF.laser.simple.bed`) and calculate the length of each region:

```bash
awk '{len=$3-$2; print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"len}' 100VAF.laser.simple.bed > regions_with_length.txt
```

This will create `regions_with_length.txt` with the following columns:

- Chromosome (`$1`)
- Start position (`$2`)
- End position (`$3`)
- Capture bias (`$4`)
- Purity (`$5`)
- Region length (`$6`)

#### 1.2 Compute Total Aligned Bases per Region

For each region, calculate the total number of bases aligned in your existing BAM file (`sim.srt.bam`):

```bash
> coverage_per_region.txt  # Create or clear the file
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' regions_with_length.txt | while read chr start end cap_bias purity length; do
    echo "Processing $chr:$start-$end"
    total_bases=$(samtools view sim.srt.bam "$chr:$start-$end" | awk '{sum += length($10)} END {print sum}')
    total_bases=${total_bases:-0}  # Default to 0 if variable is empty
    echo -e "$chr\t$start\t$end\t$cap_bias\t$purity\t$length\t$total_bases" >> coverage_per_region.txt
done
```

This will generate `coverage_per_region.txt` with:

- Columns from `regions_with_length.txt`
- Total bases aligned (`$7`)

#### 1.3 Calculate Average Coverage per Region

Compute the average coverage for each region:

```bash
awk '{coverage=$7/$6; print $0"\t"coverage}' coverage_per_region.txt > coverage_with_coverage.txt
```

Now, `coverage_with_coverage.txt` includes:

- All previous columns
- Calculated coverage (`$8`)

#### 1.4 Generate BED File of Unfinished Regions

Identify regions where the coverage is less than the desired threshold (e.g., 30x):

```bash
desired_coverage=30
awk -v dc=$desired_coverage '{if ($8 < dc*0.95) print $1"\t"$2"\t"$3"\t"$4"\t"$5}' coverage_with_coverage.txt > remaining_regions.bed
```

The `remaining_regions.bed` file will contain regions needing further simulation, formatted correctly for VISOR LASeR.

### 2. Resume VISOR LASeR Simulation

Run VISOR LASeR using the new BED file and specify a new output directory:

```bash
VISOR LASeR \
  -g ~/qgn1237/1_my_database/GRCh38_p13/GRCh38.p13.genome.fa \
  -s ../../Simulation_data/downloaded_SV_data_and_simulated_genome/visor_simulated_genome/ \
  -b remaining_regions.bed \
  -o pacbio_fastq_remaining \
  --threads 24 \
  --coverage 30 \
  --fastq \
  --read_type pacbio
```

### 3. Merge Simulation Results

#### 3.1 Merge FASTQ Files

Combine the original and newly generated FASTQ files:

```bash
cat pacbio_fastq/*.fq pacbio_fastq_remaining/*.fq > combined_reads.fq
```

If your FASTQ files are compressed (`.fq.gz`), use:

```bash
zcat pacbio_fastq/*.fq.gz pacbio_fastq_remaining/*.fq.gz | gzip > combined_reads.fq.gz
```

#### 3.2 Merge BAM Files

Ensure both BAM files are sorted and indexed before merging:

```bash
samtools merge -@ 24 combined_sim.srt.bam sim.srt.bam pacbio_fastq_remaining/sim.srt.bam
samtools index combined_sim.srt.bam
```

### 4. Verify Final Coverage

#### 4.1 Recalculate Coverage on Merged BAM

Repeat the coverage calculation steps on the merged BAM file:

```bash
> combined_coverage_per_region.txt  # Create or clear the file
awk '{len=$3-$2; print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"len}' 100VAF.laser.simple.bed > combined_regions_with_length.txt

awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' combined_regions_with_length.txt | while read chr start end cap_bias purity length; do
    echo "Processing $chr:$start-$end"
    total_bases=$(samtools view combined_sim.srt.bam "$chr:$start-$end" | awk '{sum += length($10)} END {print sum}')
    total_bases=${total_bases:-0}
    echo -e "$chr\t$start\t$end\t$cap_bias\t$purity\t$length\t$total_bases" >> combined_coverage_per_region.txt
done

awk '{coverage=$7/$6; print $0"\t"coverage}' combined_coverage_per_region.txt > combined_coverage_with_coverage.txt
```

#### 4.2 Identify Any Remaining Low-Coverage Regions

Check for regions still below the desired coverage:

```bash
desired_coverage=30
awk -v dc=$desired_coverage '{if ($8 < dc*0.95) print $1"\t"$2"\t"$3"\t"coverage}' combined_coverage_with_coverage.txt
```

If no regions are listed, all regions have achieved the desired coverage.
