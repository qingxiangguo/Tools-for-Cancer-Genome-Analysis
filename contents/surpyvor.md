
This software can merge and intersect with the logic of SUrVIVOR, considering
the breakpoints.

```bash
mamba create -n surpyvor python==3.7

mamba install -c bioconda surpyvor

# Intersect SV vcf

surpyvor highconf --variants 15X_long.vcf 15X_NGS.vcf -d 1000 -o 15X_isec.vcf
```