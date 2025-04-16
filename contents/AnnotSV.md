# 🧬 AnnotSV Cheat Sheet (Minimal)

## 🔧 Installation
```bash
conda install -c bioconda annotsv
```

## 📥 Download Human Annotations (GRCh38)
```bash
cd /your/path/
bash INSTALL_annotations.sh 
```

> This creates `/your/path/AnnotSV_annotations/`

## ✅ Basic Usage
```bash
AnnotSV \
  -SVinputFile merged_SV.vcf \
  -outputFile merged_SV_annotated.tsv \
  -genomeBuild GRCh38 \
  -annotationsDir /your/path/AnnotSV_annotations \
  -annotationMode full \
  -overwrite 1
```

## 🔑 Key Parameters

| Parameter           | Description                            |
|---------------------|----------------------------------------|
| `-SVinputFile`      | Input SV VCF file                      |
| `-outputFile`       | Output TSV file                        |
| `-genomeBuild`      | Use `GRCh38` or `GRCh37`               |
| `-annotationsDir`   | Path to downloaded annotation files    |
| `-annotationMode`   | Use `full` for all annotation layers   |
| `-overwrite`        | Set `1` to overwrite existing files    |
