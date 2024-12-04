# SVanalyzer Quick Guide: Benchmarking and Merging Structural Variants

This guide provides essential commands for using SVanalyzer's `benchmark` and `merge` functions to analyze structural variants (SVs).

---

## Benchmarking Structural Variants

Compare test SVs against a truth set to estimate sensitivity and specificity.

### Usage

```bash
svanalyzer benchmark --ref <reference.fasta> --test <test.vcf> --truth <truth.vcf>
```

### Key Options

- `--ref`: Reference FASTA file (required).
- `--test`: VCF file with test variants (required).
- `--truth`: VCF file with true variants (required).
- `--maxdist`: Maximum positional distance between variants to consider a match (default: 100,000).
- `--normshift`: Maximum normalized shift allowed (default: 0.2).
- `--normsizediff`: Maximum normalized size difference (default: 0.2).
- `--normdist`: Maximum normalized edit distance (default: 0.2).
- `--minsize`: Minimum variant size to include (default: 0).
- `--prefix`: Prefix for output files (default: `benchmark`).

### Example

```bash
svanalyzer benchmark \
  --ref reference.fasta \
  --test sample_test.vcf \
  --truth sample_truth.vcf \
  --prefix benchmark_results
```

---

## Merging Structural Variants

Group SVs by similarity to create a unified set of variants.

### Usage

**Merge from a single VCF file:**

```bash
svanalyzer merge --ref <reference.fasta> --variants <variants.vcf> --prefix <output_prefix>
```

**Merge from multiple VCF files:**

```bash
svanalyzer merge --ref <reference.fasta> --fof <vcf_list.txt> --prefix <output_prefix>
```

### Key Options

- `--ref`: Reference FASTA file.
- `--variants`: VCF file containing variants to merge.
- `--fof`: File listing paths to VCF files (one per line).
- `--prefix`: Prefix for output files (default: `merged`).
- `--maxdist`: Maximum distance between variants for merging (default: 2000).
- `--reldist`: Maximum normalized edit distance (default: 0.2).
- `--relsizediff`: Maximum normalized size difference (default: 0.2).
- `--relshift`: Maximum normalized shift (default: 0.2).
- `--seqspecific`: Only include variants with explicit REF and ALT sequences.

### Example

**Merging from multiple VCF files:**

1. **Create a file listing your VCFs (e.g., `vcf_list.txt`):**

   ```
   caller1.vcf
   caller2.vcf
   caller3.vcf
   ```

2. **Run the merge command:**

   ```bash
   svanalyzer merge \
     --ref reference.fasta \
     --fof vcf_list.txt \
     --prefix merged_results
   ```

---

## Notes

- **Reference Consistency:** Ensure the reference FASTA file matches the one used for generating your VCFs.
- **Parameter Tuning:** Adjust `--maxdist`, `--reldist`, etc., based on your data characteristics.
- **Sequence-Specific Variants:** Use `--seqspecific` to focus on variants with defined REF and ALT sequences.
- **Output Files:** The `--prefix` option sets the prefix for all output files generated.
