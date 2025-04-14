# 🧬 CTAT-LR-Fusion Cheatsheet (Long-read Fusion Detection)

## 🧰 Installation via Singularity (HPC-friendly)

```bash
# Load Singularity (required on most HPC systems)
module load singularity

# Pull the CTAT-LR-Fusion Singularity image
singularity pull docker://trinityctat/ctat_lr_fusion

# Result: ctat_lr_fusion_latest.sif
```

---

## 📥 Download and Unpack CTAT Genome Library

```bash
# Download CTAT genome lib (GRCh38 + gencode v44)
wget -c -t 100 https://.../GRCh38_gencode_v44_CTAT_lib_*.tar.gz

# Unpack the archive
tar -xzvf GRCh38_gencode_v44_CTAT_lib_*.tar.gz

# Result: genome lib directory ready (contains fasta, gtf, etc.)
```

---

## 🚀 Run CTAT-LR-Fusion (Kickstart mode from aligned BAM)

```bash
singularity exec -e -B $(pwd) -B /path/to/ctat_genome_lib \
  ctat_lr_fusion_latest.sif \
  ctat-LR-fusion \
  --LR_bam your_aligned_longread.bam \
  --genome_lib_dir /path/to/ctat_genome_lib \
  --CPU 8 \
  --vis
```

> ✅ Use kickstart mode if you already have minimap2-aligned BAM.

---

## 🧬 Run CTAT-LR-Fusion (Standard mode from FASTQ)

```bash
singularity exec -e -B $(pwd) -B /path/to/ctat_genome_lib \
  ctat_lr_fusion_latest.sif \
  ctat-LR-fusion \
  -T your_reads.fastq.gz \
  --genome_lib_dir /path/to/ctat_genome_lib \
  --CPU 8 \
  --prep_reference \
  --vis
```

> ⚠️ `--prep_reference` is only needed once per genome lib.

---

## 📄 Output Files Overview

| File | Description |
|------|-------------|
| `*.fusion_predictions.tsv` | Final filtered fusion calls |
| `*.fusion_predictions.abridged.tsv` | Simplified fusion call list |
| `*.fusion_inspector_web.html` | Interactive fusion visualization |
| `*.fusion_predictions.preliminary.tsv` | Raw, unfiltered fusions |

