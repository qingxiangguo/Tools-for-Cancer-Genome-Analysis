# The installation and usage of manta

## 1. About

 Manta calls structural variants (SVs) and indels from mapped paired-end sequencing reads.

## 2. Installation and Usage

### mamba

```bash
# Note that manta need python2.7
mamba create --name mamba_py2 python=2.7
mamba activate mamba_py2
mamba install -c bioconda manta
# The exe file of manta is configManta.py

# Or 

wget https://github.com/Illumina/manta/releases/download/v${MANTA_VERSION}/manta-${MANTA_VERSION}.release_src.tar.bz2
tar -xjf manta-${MANTA_VERSION}.release_src.tar.bz2
mkdir build && cd build
# Ensure that CC and CXX are updated to target compiler if needed, e.g.:
#     export CC=/path/to/cc
#     export CXX=/path/to/c++
../manta-1.6.0.release_src/configure --jobs=4 --prefix=/home/qgn1237/2_software/manta/
make -j4 install
```

### 2.1 run the manta config script

```bash
# The bam file has to be indexed first
# Step 1, run configManta.py first

configManta.py --tumorBam=22Rv1.5X.bam --referenceFasta=/projects/b1171/twp7981/database/genome/hg38.fa --runDir=/home/qgn1237/working/NGS_data/22Rv1/5X_depth/manta

# Successfully created workflow run script
# Step 2, run runWorkflow.py
/projects/b1171/qgn1237/4_single_cell_SV_chimera/2_test_5_prostate_depth_SV_impact/NGS_data/22Rv1/5X_depth/manta/runWorkflow.py
```

Output files: tumorSV.vcf


```bash
python MANTA_steps_generator.py \
    --bam sample.bam \
    --reference reference.fa \
    --work-dir /path/to/manta/output
```

```bash
It's really important to make sure the input bam file is non-redundant.

samtools rmdup -S HG002_NGS.bam  removed.bam
samtools index removed.bam
```