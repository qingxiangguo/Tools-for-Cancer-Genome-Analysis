# The installation and usage of SVABA

## 1. About

SvABA is a method for detecting structural variants in sequencing data using genome-wide local assembly

## SVaba Installation Guide

### Prerequisites

Before installing SVaba, ensure you have:
- Git
- A C++ compiler (GCC recommended)
- Make
- CMake (version 3.10 or higher)
- HTSlib (version 1.16 recommended)

### Step 1: Install HTSlib

HTSlib is a required dependency for SVaba. Here's how to install it in your user directory:

```bash
# Create a software directory
mkdir -p ~/software
cd ~/software

# Download HTSlib
wget https://github.com/samtools/htslib/releases/download/1.16/htslib-1.16.tar.bz2

# Extract the archive
tar -xjf htslib-1.16.tar.bz2
cd htslib-1.16

# Configure with user directory prefix
./configure --prefix=$HOME/.local

# Compile and install
make
make install
```

### Step 2: Install CMake (if version < 3.10)

SVaba requires CMake version 3.10 or higher. If your system version is older:

```bash
# Create installation directory
cd ~/software
mkdir cmake_install
cd cmake_install

# Download CMake
wget https://github.com/Kitware/CMake/releases/download/v3.25.1/cmake-3.25.1.tar.gz

# Extract and enter directory
tar -xzf cmake-3.25.1.tar.gz
cd cmake-3.25.1

# Configure and install
./bootstrap --prefix=$HOME/.local
make
make install

# Add to PATH
echo 'export PATH=$HOME/.local/bin:$PATH' >> ~/.bashrc
source ~/.bashrc   # For bash users
# Or for fish shell users:
# fish_add_path $HOME/.local/bin
```

### Step 3: Install SVaba

Now we can proceed with SVaba installation:

```bash
# Clone SVaba repository
git clone --recursive https://github.com/walaj/svaba
cd svaba

# Create and enter build directory
mkdir build
cd build

# Configure with CMake
cmake .. -DHTSLIB_DIR=$HOME/.local

# Compile
make
```

### Step 4: Set Up Environment

To make SVaba accessible from anywhere, add it to your PATH:

For Bash users:
```bash
echo 'export PATH=/path/to/svaba/build:$PATH' >> ~/.bashrc
source ~/.bashrc
```

For Fish shell users:
```fish
fish_add_path /path/to/svaba/build
```

### Testing the Installation

To verify your installation:

```bash
svaba --help
```

### Basic Usage Example

```bash
# Run tumor/normal analysis on Chr22 with 4 cores
svaba -t tumor.bam -n normal.bam -k 22 -G ref.fa -a test_id -p 4
```

## 2. Usage

### 2.1 Assume that you've already index and sort the BAM file first, you also have a BWM indexed genome

```bash
bwa index -a bwtsw /projects/b1171/qgn1237/1_my_database/GRCh38_p13/GRCh38.p13.genome.fa

svaba run -t in.bam -p 8 -a your_tag -G /projects/b1171/qgn1237/1_my_database/GRCh38_p13/GRCh38.p13.genome.fa
```

### 2.2 Output files

SRR11951439_sort.var.vcf

```bash
python SVABA_steps_generator.py --bam input.bam --reference reference.fa --tag my_analysis

python SVABA_steps_generator.py --bam input.bam --reference reference.fa --tag my_analysis --threads 16 --mem 100G --time 24:00:00
```