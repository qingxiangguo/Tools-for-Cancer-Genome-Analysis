# Samtools

## Filtering by the MAPQ (60) and output BAM file

```bash
samtools view -@ 8 -bh -q 60 test.bam -o test_q60.bam
```

\# -h       include header in SAM output
\# -b       output BAM
\# -q INT   only include reads with mapping quality >= INT [0]

## Summarize the BAM file

```bash
samtools flagstat SRR9736820_STAR_q30_sorted.bam
```

Output:

105606654 + 0 in total (QC-passed reads + QC-failed reads) # Total reads
0 + 0 secondary
221416 + 0 supplementary
0 + 0 duplicates
105606654 + 0 mapped (100.00% : N/A)
105385238 + 0 paired in sequencing
52692619 + 0 read1
52692619 + 0 read2
105015090 + 0 properly paired (99.65% : N/A)
105385238 + 0 with itself and mate mapped
0 + 0 singletons (0.00% : N/A)
276594 + 0 with mate mapped to a different chr
276594 + 0 with mate mapped to a different chr (mapQ>=5)

## Index the BAM file

```bash
samtools index SRR9736820_STAR_q30_sorted_noduplicate.bam SRR9736820_STAR_q30_sorted_noduplicate.bam.bai
```

## Sort the SAM file into a sorted BAM file

```bash
 -@ is the number of precessors
samtools sort SRR11951439.sam -o SRR11951439_sort.bam -@ 8
```

## Get the overall genome coverage from a sorted bam

$3 means the sixth column, NR means row number > 1
We add all the covered base and total bases respectively

```bash
samtools coverage SRR11951439_sort.bam | awk 'NR>1 {suma+=$5; sumb+=$3 } END { print "GenomeCoverageAverage = ",suma/sumb}'
```

## Estimate the overall genome depth (X) from a bam

```bash

samtools depth -a  *bamfile*  |  awk '{sum+=$3} END { print "Average = ",sum/{$genome_total_base_size}}'
```

The total size can be calculated like this:

```bash
samtools view -H *bamfile* | grep -P '^@SQ' | cut -f 3 -d ':' | awk '{sum+=$1} END {print sum}'
```

This command get the @SQ line from sam/bam, get the third line and sum them, then you get all the total genome size

## Subsample BAM file to a desired fraction

```bash
samtools view -s 0.156 -b ../22Rv1.bam -@ 8 > 22Rv1_5X.bam
```

## Samtools index

```bash
samtools index in.bam out.bam.bai
```

## Samtools merge

```bash
samtools merge -o PC3_10cells_exp_1_all.bam 1.bam 2.bam 3.bam -@ 6
```
