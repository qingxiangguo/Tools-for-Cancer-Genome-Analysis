# StringTie Usage Guide

Assemble transcripts and estimate their abundances:

```bash
stringtie input.bam -o output.gtf -G reference.gtf -p 8

stringtie ../minimap2/VCaP_PCLC_input.bam \
    -o VCaP_PCLC_input.gtf \
    -G ~/qgn1237/1_my_database/GRCh38_p13/gencode.v41.chr_patch_hapl_scaff.annotation.gtf \
    -L \
    -p 8
```