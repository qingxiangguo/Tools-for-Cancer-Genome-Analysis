# SeqKit: A Toolkit for FASTA/Q Files Manipulation

SeqKit is a versatile command-line toolkit for manipulating FASTA/Q files.

## Installation

To install SeqKit, run:
```
conda install -c bioconda seqkit
```

## Common Commands

**Get basic statistics of FASTA/Q files:**
```bash
seqkit stats example.fasta
```

**Extract sequences by name or pattern from FASTA/Q files:**
```bash
seqkit grep -n -r -p "pattern" example.fasta
```

**Split FASTA/Q files into multiple parts:**
```bash
seqkit split2 --by-part 10 -O output_directory yourfile.fastq
```

**Convert FASTQ to FASTA:**
```bash
seqkit fq2fa example.fastq > example.fasta
```

**Sort sequences in FASTA/Q files by name or sequence length:**
```bash
seqkit sort --by-length example.fasta
```

**Extract sub-sequences by region from sequences:**
```bash
seqkit subseq --chr chr1 --start 100 --end 500 example.fasta
```

**Find and remove duplicate sequences:**
```bash
seqkit rmdup -s example.fasta
```

**Convert DNA sequences to RNA:**
```bash
seqkit seq -t dna2rna example.fasta
```

**Concatenate multiple FASTA/Q files:**
```bash
seqkit concat example1.fasta example2.fasta > merged.fasta
```

**Randomly sample sequences from FASTA/Q files:**
```bash
seqkit sample -n 100 example.fasta
```
