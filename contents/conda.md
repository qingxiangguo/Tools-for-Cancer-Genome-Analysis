# The installation and usage of Conda (Anaconda) and Mamba
# 1. About
Package management system and environment management system
# 2. Installation and Usage
## 2.1 Download file and install
```
wget https://repo.anaconda.com/archive/Anaconda3-2022.10-Linux-x86_64.sh

bash Anaconda3-2022.10-Linux-x86_64.sh

# Press Enter

# Type yes, press Enter

# Specify the place you want to install:
/home/qgn1237/2_software/anaconda3

# Do you wish the installer to initialize Anaconda3 by running conda init?

# Type No, press Enter, in this way, you will not automatically go to conda base environment
```

## 2.2 Create an Environment, for example, you want to add a series of bio-softwares in an independent environment

It's easier to create stand-alone environments for project management, do not install any packages in the "base" environment, keep it clean unless you know how the package will affect the native environment

check how many environment you have
```
conda env list
```

Create a environment named "biosoft" with specified Python version
```
conda create --name biosoft python=3.10
```

Then activate this environment
```
conda activate biosoft
```

You can see now you are in the environment: (biosoft) [XXX@XXXX 2_software]$

You can search whether and which version of the software can be installed in the bioconda channel
```
conda search -c bioconda minimap2
```

```
# Name                       Version           Build  Channel             
minimap2                    2.0.r191               0  bioconda            
minimap2                    2.1.r311               0  bioconda            
minimap2                       2.1.1               0  bioconda            
minimap2                         2.3               0  bioconda
```

Install the specified version:
```
conda install -c bioconda minimap2=2.24
```

However, this may report some errors due to python version, then change the python version to lower:
```
conda install python=3.7
```

Exit the current environment
```
conda deactivate
```

Install mamba in the base environment
```
wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh
```

Install mamba the same way as conda

Remove the former environment
```
conda remove -n biosoft --all
```

Go to the bin directory of mambaforge， because you don't use auto init, you have to do this every time when you want to use mamba
```
./mamba init

source ~/.bashrc
```

Create a new environment with mamba
``
mamba create -n mamba666 python==3.7
``

Activate this environment
```
mamba activate mamba666
```

Install like conda
```
mamba install -c bioconda minimap2
```

```
minimap2 --version  # You succeed!
```

Check the packages and softwares you install
```
conda list
```

Uninstall software
```
mamba uninstall samtools
```
