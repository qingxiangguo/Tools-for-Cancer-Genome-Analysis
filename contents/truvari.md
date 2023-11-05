# Truvari Tutorial

## Introduction

Truvari is a command-line toolkit for structural variant (SV) comparison. It's used for benchmarking SVs against a known reference, merging SV calls, and annotating differences.

## Preparing Data

Ensure your VCF files are bgzipped, sorted, and indexed. Here's how you can prepare your data:

```bash
# Compress VCF file
bgzip -c input.vcf > input.vcf.gz

# Sort the compressed VCF
bcftools sort input.vcf.gz -Oz -o sorted_input.vcf.gz

# Index the sorted VCF
tabix -p vcf sorted_input.vcf.gz
```

## Benchmarking SVs

To benchmark SVs, use the `truvari bench` command:

```bash
# Running Truvari benchmarking
truvari bench -b base.vcf.gz -c comp.vcf.gz -o output_directory
```

Replace `base.vcf.gz` with your reference VCF file and `comp.vcf.gz` with the VCF file you want to compare.

## Results

After running the benchmarking, Truvari will generate an `output_directory` containing:

- `summary.txt`: A summary of the stats between the base and comparison.
- `tp-base.vcf.gz` and `tp-call.vcf.gz`: True positives from the base and call sets.
- `fn.vcf.gz`: False negatives.
- `fp.vcf.gz`: False positives.

To check the results, you can look into `summary.json`:

```bash
cat output_directory/summary.json
```

## Advanced Options

Truvari has several advanced options for fine-tuning the comparison, such as adjusting size filters and match criteria. Consult the [official Truvari documentation](https://github.com/spiralgenetics/truvari) for detailed usage and examples.
