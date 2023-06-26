# The installation and usage of Porechop

## 1. About

A tool for finding and removing adapters from Oxford Nanopore reads

## 2. Installation

```
git clone https://github.com/rrwick/Porechop.git
cd Porechop
python3 setup.py install
porechop -h
```

## 3. Usage

```bash
porechop -i input_reads.fastq -o output_reads.fastq -t 24
```