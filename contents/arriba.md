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


# 3. Usage

# 4. Citation
Sebastian Uhrig, Julia Ellermann, Tatjana Walther, Pauline Burkhardt, Martina Fröhlich, Barbara Hutter, Umut H. Toprak, Olaf Neumann, Albrecht Stenzinger, Claudia Scholl, Stefan Fröhling and Benedikt Brors: Accurate and efficient detection of gene fusions from RNA sequencing data. Genome Research. March 2021 31: 448-460; Published in Advance January 13, 2021. doi: 10.1101/gr.257246.119
