# BBMap repair.sh Guide

## Purpose
- Fixes unpaired FASTQ files by matching read pairs and ensuring proper pairing

## Basic Usage
```bash
repair.sh in1=input_R1.fastq.gz in2=input_R2.fastq.gz out1=fixed_R1.fastq.gz out2=fixed_R2.fastq.gz
```

## Key Parameters
- `in1`: Forward reads input file
- `in2`: Reverse reads input file  
- `out1`: Fixed forward reads output
- `out2`: Fixed reverse reads output

## Example
```bash
repair.sh in1=merged_R1.fastq.gz in2=merged_R2.fastq.gz out1=HG001_fixed_R1.fastq.gz out2=HG001_fixed_R2.fastq.gz
```

## Features
- Automatically identifies and matches read pairs
- Handles gzipped or uncompressed files
- Preserves read order
- Removes orphaned reads
