# Picard
## Help
```
java -jar ./picard.jar -h
```

## Mark and remove duplicates
```
java -jar /home/qgn1237/2_software/picard/picard.jar MarkDuplicates\
            I=input.bam \
            O=filtered_duplicates.bam \
            M=marked_dup_metrics.txt REMOVE_DUPLICATES=true
```
