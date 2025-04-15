# 🧬 Modkit Methylation Analysis Cheatsheet

**Tool:** [`modkit`](https://github.com/nanoporetech/modkit)  
**Platform:** Oxford Nanopore WGS with DNA modifications (5mC/5hmC/6mA)  
**Data:** BAM with MM/ML tags, aligned to GRCh38  

---

## 🔧 1. Setup (If using standalone binary)

```bash
# Navigate to modkit binary directory
cd /path/to/modkit/

# Make binary available
export PATH=$PWD:$PATH
```

---

## 📊 2. Summarize Methylation Calls

**Command:**

```bash
modkit summary aligned.bam > methylation_summary.txt
```

**Output preview (tabular):**

```
# total_reads_used     10042
base  code  pass_count  pass_frac  all_count  all_frac
C     m     1404887     0.045      2116171    0.061
C     h     180168      0.0058     512969     0.015
```

- `m` = 5mC  
- `h` = 5hmC  
- `a` = 6mA  
- `-` = canonical base

---

## 📌 3. Generate CpG Methylation Profile (bedMethyl format)

```bash
modkit pileup aligned.bam cpg_methylation.bed \
  --ref reference.fa \
  --preset traditional \
  --threads 16
```

- `--preset traditional`: focuses on CpG 5mC, combines strands, filters 5hmC
- Output columns include position, mod counts, fractions, and coverage.

---

## 📈 4. Generate bedGraph for Visualization (5mC)

```bash
modkit pileup aligned.bam output_dir \
  --ref reference.fa \
  --preset traditional \
  --threads 16 \
  --bedgraph --prefix sample_m5C
```

- Outputs:
  - `sample_m5C.5mC_pos.bg`
  - `sample_m5C.5mC_neg.bg`
- Can be converted to `.bw` via `bedGraphToBigWig` for IGV/jBrowse

---

## 🧬 5. Region-Specific Methylation (Optional)

### CpG Islands

```bash
modkit pileup --region-bed cpg_islands.bed \
  aligned.bam cpg_islands_methylation.bed \
  --ref reference.fa \
  --mod m5C \
  --threads 16
```

### Promoters

```bash
modkit pileup --region-bed promoters.bed \
  aligned.bam promoter_methylation.bed \
  --ref reference.fa \
  --mod m5C \
  --threads 16
```

---

## 📂 Output Files Overview

| File                          | Description                                |
|-------------------------------|--------------------------------------------|
| `methylation_summary.txt`     | Global counts & fractions for each mod     |
| `cpg_methylation.bed`         | Genome-wide methylation stats (bedMethyl)  |
| `sample_m5C.*.bg`             | Strand-specific bedGraph (5mC)             |
| `*_region_methylation.bed`    | Targeted methylation for region of interest|
