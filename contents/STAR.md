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



# 3. 
# 4. Citation
Spliced Transcripts Alignment to a Reference © Alexander Dobin, 2009-2022 https://www.ncbi.nlm.nih.gov/pubmed/23104886
