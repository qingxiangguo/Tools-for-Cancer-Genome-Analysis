# Samtools
## Filtering by the MAPQ (60) and output BAM file

```
samtools view -@ 8 -bh -q 60 test.bam -o test_q60.bam
```

\# -h       include header in SAM output
\# -b       output BAM
\# -q INT   only include reads with mapping quality >= INT [0]

## Summarize the BAM file
```
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
```
samtools index SRR9736820_STAR_q30_sorted_noduplicate.bam SRR9736820_STAR_q30_sorted_noduplicate.bam.bai
```

## Sort the SAM file into a sorted BAM file
```
# -@ is the number of precessors
samtools sort SRR11951439.sam -o SRR11951439_sort.bam -@ 8
```

## Get the overall genome coverage from a sorted bam
```
# $6 means the sixth column, NR means row number > 1, n++ means the index add 1 everytime
samtools coverage SRR11951439_sort.bam | awk 'NR>1 {sum+=$6; n++} END { print "GenomeAverageCoverage= ",sum/n}'
```




