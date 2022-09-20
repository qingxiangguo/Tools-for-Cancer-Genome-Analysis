# The installation and usage of Arriba
# 1. About
Arriba is a command-line tool for the detection of gene fusions from RNA-Seq data. It is based on the ultrafast STAR aligner. Arriba does not require to reduce the STAR parameter --alignIntronMax to detect fusions arising from focal deletions. Reducing this parameter impairs mapping of reads to genes with long introns and may affect expression quantification, hence.  
Apart from gene fusions, Arriba can detect other structural rearrangements with potential clinical relevance, including viral integration sites, internal tandem duplications, whole exon duplications, intragenic inversions, enhancer hijacking events involving immunoglobulin/T-cell receptor loci, translocations affecting genes with many paralogs such as DUX4, and truncations of genes (i.e., breakpoints in introns or intergenic regions).  

# 2. Installation and quick start
Arriba has only a single prerequisite: STAR (version >=2.7.10a recommended). Download and install the tool according to the developers' instructions and make it available in your $PATH.  
```
wget https://github.com/suhrig/arriba/releases/download/v2.3.0/arriba_v2.3.0.tar.gz
tar -xzf arriba_v2.3.0.tar.gz
cd arriba_v2.3.0 && make # or use precompiled binaries
```
Arriba requires an assembly in FastA format, gene annotation in GTF format, and a STAR index built from the two. You can use your preferred assembly and annotation, as long as their coordinates are compatible with hg19/hs37d5/GRCh37 or hg38/GRCh38 or mm10/GRCm38 or mm39/GRCm39. If you use another assembly, then the coordinates in the blacklist will not match and the predictions will contain many false positives.  

ENCODE annotation is recommended over RefSeq due to more comprehensive annotation of immunoglobulin/T-cell receptor loci and splice sites.

If you do not already have the files and a STAR index, you can use the script download_references.sh. It downloads the files to the current working directory and builds a STAR index. Run the script without arguments to see a list of available files. Note that this step requires ~45 GB of RAM and 8 cores (can be adjusted by setting the environment variable THREADS).
```
./download_references.sh hs37d5viral+GENCODE19
```

The download file contains a script run_arriba.sh, which demonstrates the usage of Arriba (see also section Workflow). We recommend that you use this as a guide to integrate Arriba into your existing STAR-based RNA-Seq pipeline. 

Run the demo script with 8 threads.
```
./run_arriba.sh STAR_index_hs37d5viral_GENCODE19/ GENCODE19.gtf hs37d5viral.fa database/blacklist_hg19_hs37d5_GRCh37_v2.3.0.tsv.gz database/known_fusions_hg19_hs37d5_GRCh37_v2.3.0.tsv.gz database/protein_domains_hg19_hs37d5_GRCh37_v2.3.0.gff3 8 test/read1.fastq.gz test/read2.fastq.gz
```

<b>Output files</b>
* fusions.tsv: This is the main output file of Arriba and contains a list of fusion candidates sorted from highest to lowest confidence. Low-confidence fusions should be considered speculative unless there is additional evidence, e.g., when the fusion is a recurrent event in the given entity. For cohort-analyses, you might want to ignore low-confidence fusions.
* fusions.discarded.tsv: This file lists all events which were classified as artifacts or which are observed in healthy tissue (e.g., read-through transcripts and non-canonically spliced transcripts).

# 3. Usage
## 3.1 The logic of Arriba
* Detection of chimeric reads must be enabled in STAR by specifying the parameter --chimSegmentMin. In addition, the parameter --chimOutType WithinBAM must be specified to cause STAR to report chimeric reads as supplementary alignments in the main output file Aligned.out.sam. Old versions of STAR (or when STAR is run with --chimOutType SeparateSAMold) wrote supplementary alignments to a separate file named Chimeric.out.sam. Arriba is compatible with this mode of use (see parameter -c), but it is deprecated, because STAR might not support it anymore in the future.  
* Arriba extracts the supplementary alignments from the given input file(s). The supplementary alignments represent evidence about translocations, inversions, duplications, and deletions larger than the usual intron size (as defined by the parameter --alignIntronMax). In order to find fusions arising from deletions smaller than the maximum intron size, Arriba also extracts alignments which cross the boundaries of annotated genes. Once all alignments have been extracted, it applies a set of filters to remove artifacts and transcripts observed in healthy tissue. The final output is a list of fusion predictions which pass all of Arriba's filters.  
  
## 3.2 Demo script  
The following parameters related to chimeric alignment are recommended for improved sensitivity (requires STAR version 2.7.6a or higher):

```
--outFilterMultimapNmax 50 \
--peOverlapNbasesMin 10 \
--alignSplicedMateMapLminOverLmate 0.5 \
--alignSJstitchMismatchNmax 5 -1 5 5 \
--chimSegmentMin 10 \
--chimOutType WithinBAM HardClip \
--chimJunctionOverhangMin 10 \
--chimScoreDropMax 30 \
--chimScoreJunctionNonGTAG 0 \
--chimScoreSeparation 1 \
--chimSegmentReadGapMax 3 \
--chimMultimapNmax 50
```  

Arriba does not care if the BAM files are sorted/indexed or not. For maximum speed, STAR's output should be piped directly to Arriba, such that chimeric reads are extracted while STAR is still running. 

A complete call of STAR in conjunction with Arriba may look like this:

```
STAR \
    --runThreadN 8 \
    --genomeDir /path/to/STAR_index --genomeLoad NoSharedMemory \
    --readFilesIn read1.fastq.gz read2.fastq.gz --readFilesCommand zcat \
    --outStd BAM_Unsorted --outSAMtype BAM Unsorted --outSAMunmapped Within --outBAMcompression 0 \
    --outFilterMultimapNmax 50 --peOverlapNbasesMin 10 --alignSplicedMateMapLminOverLmate 0.5 --alignSJstitchMismatchNmax 5 -1 5 5 \
    --chimSegmentMin 10 --chimOutType WithinBAM HardClip --chimJunctionOverhangMin 10 --chimScoreDropMax 30 \
    --chimScoreJunctionNonGTAG 0 --chimScoreSeparation 1 --chimSegmentReadGapMax 3 --chimMultimapNmax 50 |
arriba \
    -x Aligned.out.bam \
    -o fusions.tsv -O fusions.discarded.tsv \
    -a /path/to/assembly.fa -g /path/to/annotation.gtf \
    -b /path/to/blacklist.tsv.gz -k /path/to/known_fusions.tsv.gz -t /path/to/known_fusions.tsv.gz -p /path/to/protein_domains.gff3
```

Note: In this example, the same file is passed to the parameters -k and -t, because it is used for two purposes: applying sensitive filtering parameters to known fusions (-k) and tagging known fusions in the tags column (-t). However, it is possible to use different files for these two parameters if a user wants to separate the two tasks.

## 3.3 Inputfiles
Arriba takes the main output file of STAR (Aligned.out.bam) as input (parameter -x). If STAR was run with the parameter --chimOutType WithinBAM, then this file contains all the information needed by Arriba to find fusions. When STAR was run with the parameter --chimOutType SeparateSAMold, the main output file lacks chimeric alignments. Instead, STAR writes them to a separate output file named Chimeric.out.sam. In this case, the file needs to be passed to Arriba via the parameter -c in addition to the main output file Aligned.out.bam.





# 4. Citation
Sebastian Uhrig, Julia Ellermann, Tatjana Walther, Pauline Burkhardt, Martina Fröhlich, Barbara Hutter, Umut H. Toprak, Olaf Neumann, Albrecht Stenzinger, Claudia Scholl, Stefan Fröhling and Benedikt Brors: Accurate and efficient detection of gene fusions from RNA sequencing data. Genome Research. March 2021 31: 448-460; Published in Advance January 13, 2021. doi: 10.1101/gr.257246.119
