# The installation and usage of Duplex tools

# 1. About

for duplex calling of guppy, needs the output of first round of guppy simplex calling.

# 2. Installation and Usage

## Installation

```bash
pip install duplex-tools
```

## Usage

```bash
duplex_tools pairs_from_summary /path/to/your/output/basecalled_fastq/sequencing_summary.txt /path/to/your/output/
```

Above step will generate pair_ids.txt.

```bash
duplex_tools filter_pairs /path/to/your/output/pair_ids.txt /path/to/your/simplex_fastq_output/basecalled_fastq
```

Get pair_ids_filtered.txt, it is for downstream duplex basecalling.
