# Picard
## Help
```
java -jar ./picard.jar -h
```

## Mark and remove duplicates
```
java -jar /home/qgn1237/2_software/picard/picard.jar MarkDuplicates\
            I=/projects/b1171/qgn1237/3_scanneo2_pipeline/1_STAR_align/SRR9736820_STAR_q30.bam \
            O=/projects/b1171/qgn1237/3_scanneo2_pipeline/1_STAR_align/SRR9736820_STAR_q30_noduplicate.bam \
            M=marked_dup_metrics.txt REMOVE_DUPLICATES=true
```
