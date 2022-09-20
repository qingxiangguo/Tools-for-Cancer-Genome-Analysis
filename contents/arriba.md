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
### 3.3.1 Alignments  
Arriba takes the main output file of STAR (Aligned.out.bam) as input (parameter -x). If STAR was run with the parameter --chimOutType WithinBAM, then this file contains all the information needed by Arriba to find fusions. When STAR was run with the parameter --chimOutType SeparateSAMold, the main output file lacks chimeric alignments. Instead, STAR writes them to a separate output file named Chimeric.out.sam. In this case, the file needs to be passed to Arriba via the parameter -c in addition to the main output file Aligned.out.bam.

<b>Deletion or Splicing?</b> Arriba extracts alignments which cross the boundaries of annotated genes, because these alignments might arise from focal deletions. In RNA-Seq data deletions of up to several hundred kb are hard to distinguish from splicing. They are represented identically as gapped alignments, because the sizes of many introns are in fact of this order of magnitude.

STAR applies a rather arbitrary measure to decide whether a gapped alignment arises from splicing or from a genomic deletion: The parameter --alignIntronMax <b>(the size of a normal intron, if larger, it might be a deletion)</b> determines what gap size is still assumed to be a splicing event and introns are used to represent these gaps. Only gaps larger than this limit are classified as potential evidence for genomic deletions and are stored as chimeric alignments. Most STAR-based fusion detection tools only consider chimeric alignments as evidence for gene fusions and are blind to focal deletions, hence. 

As a workaround, these tools recommend reducing the value of the parameter --alignIntronMax <b>(to make it more sensible to a deletion)</b>. But this impairs the quality of alignment, because it reduces the scope that STAR searches to find a spliced alignment <b>(but too small intron size will sacrifice the splicing detection)</b>. To avoid compromising the quality of alignment for the sake of fusion detection, the only solution would be to run STAR twice - once with settings optimized for regular alignment and once for fusion detection. This would double the runtime. In contrast, Arriba does not require to reduce the maximum intron size. It employs a more sensible criterion to distinguish splicing from deletions: Arriba considers all those reads as potential evidence for deletions that span the boundary of annotated genes.

### 3.3.2 Assembly  
Arriba takes the assembly as input (parameter -a).  Any assembly can be used as long as its coordinates are compatible with one of the supported assemblies (hg19/hs37d5/GRCh37, hg38/GRCh38, mm10/GRCm38, mm39/GRCm39).

### 3.3.3 Annotation  
GENCODE annotation is recommended over RefSeq annotation, because the former has a more comprehensive annotation of transcripts and splice-sites, which boosts the sensitivity. The file must be provided in GTF format.

### 3.3.4 Blacklist  
It is strongly advised to run Arriba with a blacklist (parameter -b). The blacklist removes recurrent alignment artifacts and transcripts which are present in healthy tissue. 

### 3.3.5 Known fusions  
Arriba can be instructed to be particularly sensitive towards events between certain gene pairs by supplying a list of gene pairs (parameter -k). 

### 3.3.6 Tags
Arriba can be supplied with a list of user-defined tags using the parameter -t. This feature is useful to annotate known oncogenic fusions, for example. It is constructed in a way that it can be passed as arguments to the parameters -k and -t alike.

### 3.3.7 Protein domains
Protein domain annotation can be passed to Arriba via the parameter -p. The column retained_protein_domains of Arriba's output file is then populated accordingly.The file must be in GFF3 format.

### 3.3.8 Structural variant calls from WGS
If whole-genome sequencing (WGS) data is available, the sensitivity and specificity of Arriba can be improved by passing a list of structural variants detected from WGS to Arriba via the parameter -d. If the structural variant calls were obtained from whole-exome sequencing (WES) data rather than WGS data, the filter no_genomic_support should be disabled, since WES has poor coverage in most regions of the genome, such that many structural variants are missed.Two file formats are accepted: a simple four-column format and the standard Variant Call Format (VCF). The format is detected automatically.

In case of the simple format, the file must contain four columns separated by tabs. The first two columns contain the breakpoints of the structural variants in the format CONTIG:POSITION. The last two columns contain the orientation of the breakpoints. The accepted values are:

* downstream or +: the fusion partner is fused downstream of the breakpoint, i.e., at a coordinate higher than the breakpoint
(This means the fusion gene is from the breakpoint to the right side)

* upstream or -: the fusion partner is fused at a coordinate lower than the breakpoint
(This means the fusion gene is from the left side to the breakpoint)

Example:
```
1:54420491  6:9248349   +   -
20:46703288 20:46734546 -   +
17:61499820 20:45133874 +   +
3:190967119 7:77868317  -   -
```

## 3.4 Outputfiles
### 3.4.1 fusions.tsv
It should be highly enriched for true predictions. Details are below:

* gene1 and gene2 : Gene1 contains the gene which makes up the 5' end of the transcript and gene2 the gene which makes up the 3' end. If a breakpoint is in an intergenic region, Arriba lists the closest genes upstream and downstream from the breakpoint, separated by a comma. The numbers in parentheses after the closest genes state the distance to the genes. If no genes are annotated for a contig (e.g., for viral genomes), the column contains a dot (.).

* strand1(gene/fusion) and strand2(gene/fusion) :  Each of these columns contains two values seperated by a slash. The strand before the slash reflects the strand of the gene according to the gene annotation supplied to Arriba via the parameter -g. 

If the breakpoint is in an intergenic region, the value is .. (because it contains no gene)

The value after the slash reflects the strand that is transcribed. This does not necessarily match the strand of the gene, namely when the sense strand of a gene serves as the template for transcription （<b>sometines the non-template (sense) strand can also be transcribed, so it is not guranteed that the strand of a transcript must be the same with the gene </b>）. Occassionally, the strand that is transcribed cannot be predicted reliably. In this case, Arriba indicates the lack of information as a dot (.). 

* breakpoint1 and breakpoint2 :  The columns contain the coordinates of the breakpoints in gene1 and gene2, respectively. 

* site1 and site2 :  These columns add information about the location of the breakpoints. Possible values are: 5' UTR, 3' UTR, UTR (overlapping with a 5' UTR as well as a 3' UTR), CDS (coding sequence), exon, intron, and intergenic. 








# 4. Citation
Sebastian Uhrig, Julia Ellermann, Tatjana Walther, Pauline Burkhardt, Martina Fröhlich, Barbara Hutter, Umut H. Toprak, Olaf Neumann, Albrecht Stenzinger, Claudia Scholl, Stefan Fröhling and Benedikt Brors: Accurate and efficient detection of gene fusions from RNA sequencing data. Genome Research. March 2021 31: 448-460; Published in Advance January 13, 2021. doi: 10.1101/gr.257246.119
