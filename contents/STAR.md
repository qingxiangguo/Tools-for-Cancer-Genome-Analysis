# The installation and usage of STAR
# 1. About
RNA-seq aligner.  STAR can discover non-canonical splices and chimeric (fusion) transcripts, and is also capable of mapping full-length RNA sequences. Extremely fast (also does splice alignment, requires at least 30 Gb memory.
# 2. Installation and Usage
## Download the latest release from and uncompress it
```
# Get latest STAR source from releases
wget https://github.com/alexdobin/STAR/archive/2.7.10a.tar.gz
tar -xzf 2.7.10a.tar.gz
cd STAR-2.7.10a
# Compile
cd STAR/source
make STAR
```

## 2.1 Build genome index
The genome indexes are saved to disk and need only be generated once for each genome/annotation combination.  It is strongly recommended that users generate their own genome indexes with most up-to-date assemblies and annotations.Very importantly, chromosome names in the annotations GTF file have to match chromosome names in the FASTA genome sequence files.

```
STAR --runThreadN 6 --runMode genomeGenerate --genomeDir ~/reference/index/STAR/mm10/ --genomeFastaFiles ~/reference/genome/mm10/GRCm38.p5.genome.fa --sjdbGTFfile ~/annotation/mm10/gencode.vM13.annotation.gtf --sjdbOverhang 100
```

--runThreadN option defines the number of threads to be used for genome generation, it has
to be set to the number of available cores on the server node.

--runMode genomeGenerate option directs STAR to run genome indices generation job.

--genomeDir specifies path to the directory

--genomeFastaFiles specifies one or more FASTA files with the genome reference sequences.

--sjdbGTFfile specifies the path to the file with annotated transcripts in the standard GTF format. STAR will extract splice junctions from this file and use them to greatly improve accuracy of the mapping. While this is optional, and STAR can be run without annotations, using annotations is highly recommended whenever they are available.

## 2.2 Running mapping jobs

```
STAR --runThreadN 20 --genomeDir ~/reference/index/STAR/mm10/ --readFilesIn SRR3589959_1.fastq SRR3589959_2.fastq --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ./SRR3589959
```

--runThreadN NumberOfThreads

--genomeDir /path/to/genomeDir

--readFilesIn /path/to/read1 [/path/to/read2 ]

--outFileNamePrefix /path/to/output/dir/prefix You can change the file prefixes

--outSAMtype BAM Unsorted output unsorted Aligned.out.bam file

--outSAMtype BAM SortedByCoordinate 

--outSAMtype BAM SortedByCoordinate output sorted by coordinate Aligned.sortedByCoord.out.bam file, similar to samtools sort command. If this option causes problems, it is recommended to reduce --outBAMsortingThreadN from the default 6 to lower values (as low as 1)

## 2.3 Output files

SRR3589959Aligned.sortedByCoord.out.bam   output sorted by coordinate Aligned.sortedByCoord.out.bam file, similar to samtools sort command.

SRR3589959Log.final.out   summary mapping statistics after mapping job is complete, very useful for quality control. The statistics are calculated for each read (single- or paired-end) and then summed or averaged over all reads.

SRR3589959Log.out   main log file with a lot of detailed information about the run. This file is most useful for troubleshooting and debugging.

SRR3589959Log.progress.out    reports job progress statistics, such as the number of processed reads, % of mapped reads etc. It is updated in 1 minute intervals.

SRR3589959SJ.out.tab   SJ.out.tab contains high confidence collapsed splice junctions in tab-delimited format. 

<b> Splice junctions </b> Note that STAR defines the junction start/end as intronic bases, while many other software define them as exonic bases. The columns have the following meaning:

column 1: chromosome

column 2: first base of the intron (1-based)  

column 3: last base of the intron (1-based)

column 4: strand (0: undefined, 1: +, 2: -)  

column 5: intron motif: 0: non-canonical; 1: GT/AG, 2: CT/AC, 3: GC/AG, 4: CT/GC, 5: AT/AC, 6: GT/AT  

column 6: 0: unannotated, 1: annotated in the splice junctions database. Note that in 2-pass mode, junctions detected in the 1st pass are reported as annotated, in addition to annotated junctions from GTF.  

column 7: number of uniquely mapping reads crossing the junction  

column 8: number of multi-mapping reads crossing the junction  

column 9: maximum spliced alignment overhang  

# 3. Citation
Spliced Transcripts Alignment to a Reference © Alexander Dobin, 2009-2022 https://www.ncbi.nlm.nih.gov/pubmed/23104886
