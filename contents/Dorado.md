# Dorado Basecalling Cheat Sheet (v0.9.1)

> If you've ever felt like the command line was gaslighting you — you're not alone.  
> This guide skips the surprises and gives you working examples for WGS, cDNA, and direct-RNA.

This guide shows how to use [Dorado](https://github.com/nanoporetech/dorado) to basecall Oxford Nanopore data with optional **modification calling**, using real tested commands.

---

## 🧬 1. Whole-Genome Sequencing (WGS) with DNA Modifications

Call **5mC/5hmC in all sequence contexts** and **6mA** using pre-downloaded models.

```bash
/home/qgn1237/2_software/dorado-0.9.1-linux-x64/bin/dorado basecaller \
  /home/qgn1237/2_software/dorado-0.9.1-linux-x64/bin/dna_r10.4.1_e8.2_400bps_hac@v5.0.0 \
  --modified-bases-models \
  /home/qgn1237/2_software/dorado-0.9.1-linux-x64/bin/dna_r10.4.1_e8.2_400bps_hac@v5.0.0_5mC_5hmC@v3,\
/home/qgn1237/2_software/dorado-0.9.1-linux-x64/bin/dna_r10.4.1_e8.2_400bps_hac@v5.0.0_6mA@v3 \
  --emit-moves \
  --device cuda:all \
  ../merged.pod5 \
  > WGS_mod_calls.bam
  ```

✅ Detects **cytosine (C)** methylation across genome (5mC/5hmC)  
✅ Detects **adenine (A)** methylation (6mA)  
✅ Compatible with downstream DMR, modkit, or IGV visualization.

## 🔬 2. direct-cDNA Basecalling with poly(A) Tail Estimation

Use `--estimate-poly-a` to measure poly(A) tail lengths from cDNA reads.

```bash
/home/qgn1237/2_software/dorado-0.9.1-linux-x64/bin/dorado basecaller \
  /home/qgn1237/2_software/dorado-0.9.1-linux-x64/bin/dna_r10.4.1_e8.2_400bps_hac@v5.0.0 \
  ../merged.pod5 \
  --estimate-poly-a \
  --device cuda:all \
  > cDNA_calls.bam
```

✅ Adds `pt:i` tag to each read for estimated poly(A) tail length  
✅ Recommended for PCS/PCB kits in full-length transcript profiling

## 🧫 3. direct-cDNA Basecalling with RNA Modifications

```bash
/gpfs/home/qgn1237/2_software/dorado-0.9.1-linux-x64/bin/dorado basecaller \
  /home/qgn1237/2_software/dorado-0.9.1-linux-x64/bin/dna_r10.4.1_e8.2_400bps_hac@v5.0.0 ../pod5/ \
  --estimate-poly-a \
  --device cuda:all \
  > out.bam 
```

## 🧫 4. direct-RNA Basecalling with RNA Modifications

Example: Call **m6A (in DRACH motifs)** and **m5C** from direct-RNA data.

```bash
/home/qgn1237/2_software/dorado-0.9.1-linux-x64/bin/dorado basecaller \
  /home/qgn1237/2_software/dorado-0.9.1-linux-x64/bin/rna004_130bps_sup@v5.1.0/ \
  ../merged.pod5 \
  --modified-bases m5C m6A_DRACH \
  --emit-moves \
  --device cuda:all \
  > RNA_mod_calls.bam
```

✅ Supports SQK-RNA004 direct-RNA kits  
✅ Output contains `MM` / `ML` tags for modification probability  
✅ Compatible with modkit and nanocompore downstream tools

## 5. Dorado aligner, better than minimap2 aloneo

```bash
mamba activate mamba666

/gpfs/home/qgn1237/2_software/dorado-0.9.1-linux-x64/bin/dorado aligner \
  ~/qgn1237/1_my_database/GRCh38_p13/GRCh38.p13.genome.fa \
  ./PCa_50904_WGS_calls_5mCG.bam \
  | samtools sort -@ 16 -o ./50904_WGS_meth_aligned.bam

samtools index -@ 16 ./50904_WGS_meth_aligned.bam
```

---

## 🛠️ Useful Tips

**Convert BAM to FASTQ**:

```bash
samtools fastq calls.bam > output.fastq
```

**View modification tags**:

```bash
samtools view calls.bam | head -n 5
```

**Extract reads by modification type (e.g., 6mA only)**:

Use [`modkit`](https://github.com/nanoporetech/modkit) or convert to BED with `modbam2bed`.

---

## ✅ Best Practices Summary

| Task | Recommended Dorado options |
|------|-----------------------------|
| WGS + DNA modifications | `--modified-bases-models` + downloaded models |
| cDNA + polyA tail length | `--estimate-poly-a` |
| direct-RNA + base mods | `--modified-bases m6A_DRACH m5C` |
| All basecalling | Always add `--emit-moves` for downstream support |
| GPU usage | `--device cuda:all` for multi-GPU, or `cuda:0` to specify |
