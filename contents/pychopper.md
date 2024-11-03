mamba install pychopper

```bash
# For cDNA-PCR, whole QC process
pychopper -k PCS114 -m edlib -r pychopper_report.pdf -u unclassified.fastq -w rescued.fastq -S statistics.tsv -A alignment_scores.bed -t 4 ../raw.fastq ./full_length.fastq

pychopper -k PCS114 -x DCS114 -m edlib -r pychopper_report_2.pdf -u unclassified_2.fastq -w rescued_2.fastq -S statistics_2.tsv -A alignment_scores_2.bed -t 4 ./unclassified.fastq ./VCaP_PCLC_output_full_length_2.fastq

cat VCaP_PCLC_output_full_length.fastq VCaP_PCLC_output_full_length_2.fastq rescued.fastq rescued_2.fastq > combined_full_length.fastq

cutadapt -g TTTTTTTTCCTGTACTTCGTTCAGTTACGTATTGCT -a AGCAATACGTAACTGAACGAAGTACAGGAAAAAAAA -g TTTCTGTTGGTGCTGATATTGCTGCCATTACGGCCGGG -a CCCGGCCGTAATGGCAGCAATATCAGCACCAACAGAAA -a "A{20}" -g "T{20}" -a "T{20}" -g "A{20}"  --error-rate=0.1 --times=2 --poly-a -o final_trimmed_full_length.fastq ./combined_full_length.fastq

minimap2 -ax splice -k14 --MD -t 8 -Y -R '@RG\tID:PC310cells\tPL:ont\tLB:library\tSM:PC310cells' /home/qgn1237/qgn1237/1_my_database/GRCh38_p13/minimap2_index/GRCh38.p13.genome.mmi ../pychopper/final_trimmed_full_length.fastq | samtools sort -@ 8 -m 2G -O BAM -o VCaP_PCLC_input_trimmed_full_length.bam && samtools index VCaP_PCLC_input_trimmed_full_length.bam  VCaP_PCLC_input_trimmed_full_length.bam.bai
```