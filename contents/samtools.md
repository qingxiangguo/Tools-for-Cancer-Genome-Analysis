# Samtools
## Filtering by the MAPQ (60) and output BAM file

```
samtools view -@ 8 -bh -q 60 test.bam -o test_q60.bam
```

\# -h       include header in SAM output
\# -b       output BAM
\# -q INT   only include reads with mapping quality >= INT [0]
