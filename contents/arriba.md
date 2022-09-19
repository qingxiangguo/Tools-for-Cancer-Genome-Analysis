# The installation and usage of Arriba
# 1. About
Arriba is a command-line tool for the detection of gene fusions from RNA-Seq data. It is based on the ultrafast STAR aligner. Arriba does not require to reduce the STAR parameter --alignIntronMax to detect fusions arising from focal deletions. Reducing this parameter impairs mapping of reads to genes with long introns and may affect expression quantification, hence.  
Apart from gene fusions, Arriba can detect other structural rearrangements with potential clinical relevance, including viral integration sites, internal tandem duplications, whole exon duplications, intragenic inversions, enhancer hijacking events involving immunoglobulin/T-cell receptor loci, translocations affecting genes with many paralogs such as DUX4, and truncations of genes (i.e., breakpoints in introns or intergenic regions).  

# 2. Installation and quick start
Arriba has only a single prerequisite: STAR (version >=2.7.10a recommended). Download and install the tool according to the developers' instructions and make it available in your $PATH.
'''
wget https://github.com/suhrig/arriba/releases/download/v2.3.0/arriba_v2.3.0.tar.gz
tar -xzf arriba_v2.3.0.tar.gz
cd arriba_v2.3.0 && make # or use precompiled binaries
'''
# 3. Usage

# 4. Citation
Sebastian Uhrig, Julia Ellermann, Tatjana Walther, Pauline Burkhardt, Martina Fröhlich, Barbara Hutter, Umut H. Toprak, Olaf Neumann, Albrecht Stenzinger, Claudia Scholl, Stefan Fröhling and Benedikt Brors: Accurate and efficient detection of gene fusions from RNA sequencing data. Genome Research. March 2021 31: 448-460; Published in Advance January 13, 2021. doi: 10.1101/gr.257246.119
