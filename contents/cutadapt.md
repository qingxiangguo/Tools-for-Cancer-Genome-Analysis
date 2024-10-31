# The usage

## 1. About

Cutadapt finds and removes adapter sequences, primers, poly-A tails and other types of unwanted sequence from your high-throughput sequencing reads.

## 2. Usage

```bash
# for nanopore cDNA-PCR:
cutadapt -g TTTTTTTTCCTGTACTTCGTTCAGTTACGTATTGCT -a AGCAATACGTAACTGAACGAAGTACAGGAAAAAAAA -g TTTCTGTTGGTGCTGATATTGCTGCCATTACGGCCGGG -a CCCGGCCGTAATGGCAGCAATATCAGCACCAACAGAAA -a "A{20}" -g "T{20}" -a "T{20}" -g "A{20}"  --error-rate=0.1 --times=2 --poly-a -o cleaned_unmapped.fa ./raw_unmapped.fa
```