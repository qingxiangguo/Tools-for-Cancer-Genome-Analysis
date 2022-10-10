# Picard
## Help
```
java -jar ./picard.jar -h
```
## Sort files by coordinate
BAM is compressed. Sorting helps to give a better compression ratio because similar sequences are grouped together.
```
java -jar /home/qgn1237/2_software/picard/picard.jar SortSam \
            I=/projects/b1171/qgn1237/3_scanneo2_pipeline/1_STAR_align/SRR9736820_STAR_q30.bam \
            O=/projects/b1171/qgn1237/3_scanneo2_pipeline/1_STAR_align/SRR9736820_STAR_q30_sorted.bam \
            SORT_ORDER=coordinate
```


## Mark and remove duplicates  
It need the BAM file to be sorted (whether by coordinate or query name).

```
java -jar /home/qgn1237/2_software/picard/picard.jar MarkDuplicates\
            I=/projects/b1171/qgn1237/3_scanneo2_pipeline/1_STAR_align/SRR9736820_STAR_q30_sorted.bam \
            O=/projects/b1171/qgn1237/3_scanneo2_pipeline/1_STAR_align/SRR9736820_STAR_q30_sorted_noduplicate.bam \
            M=marked_dup_metrics.txt REMOVE_DUPLICATES=true
```
