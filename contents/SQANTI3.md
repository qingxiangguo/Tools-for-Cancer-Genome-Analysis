# SQANTI Usage Guide

Evaluate transcript quality and structure:
```
# You need to fix the GFT file produced by Stringtie first
awk '$7 != "."' VCaP_PCLC_input.gtf > VCaP_PCLC_input_fixed.gtf


sqanti3_qc.py input.gtf reference_annotation.gtf reference_genome.fa --coverage coverage.bam --report pdf --output_dir output_directory

sqanti3_qc.py ./VCaP_direct_cDNA_fixed.gtf gencode.v38.annotation.gtf hg38.fa --report pdf -t 4 --force_id_ignore --isoform_hits
```

## Output Explanation
- **FSM (Full-Splice Match)**: Transcripts that match a reference annotation across the entire length including all splice junctions.
- **ISM (Incomplete-Splice Match)**: Transcripts that match some but not all splice junctions of any reference transcript, indicating potential degradation or incomplete assembly.
