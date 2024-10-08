
# Using featureCounts for Transcriptome Quality Assessment

## 1. About

featureCounts is a high-performance read summarization program that counts mapped reads for genomic features such as genes, exons, promoters, and genomic bins.

## 2. Installation and Usage

### 2.1 Installation

Install featureCounts using conda:

```bash
conda install -c bioconda subread
```

Or download from the [Subread website](http://subread.sourceforge.net/).

### 2.2 Basic Usage

Run featureCounts:

```bash
featureCounts -a reference.gtf -o transcriptome_counts.txt sample.bam
```

### 2.3 Analyzing Output for Transcriptome Quality

1. Examine summary file:
   ```bash
   cat transcriptome_counts.txt.summary
   ```

2. Count expressed genes:
   ```bash
   awk 'NR>2 {if ($7 > 0) count++} END {print "Number of detected expressed genes:", count}' transcriptome_counts.txt
   ```

3. Calculate average reads per gene:
   ```bash
   awk 'NR>2 {sum+=$7; count++} END {print "Average number of reads per gene:", sum/count}' transcriptome_counts.txt
   ```

## 3. Interpreting Results

1. Summary file:
   - "Assigned" should ideally be 60-80% or higher.
   - High "Unassigned_NoFeatures" may indicate annotation mismatch.
   - High "Unassigned_MultiMapping" suggests duplicate sequences or gene family issues.

2. Expressed genes:
   - 15,000-20,000 for human samples is generally good.

3. Average reads per gene:
   - 30-50 reads per gene is typically sufficient.
   - Higher numbers indicate better sequencing depth.

4. Median vs Average expression:
   - Compares distribution skew.


