# Tools-and-tricks-for-Cancer-Genome-Analysis

Installation and usage for various tools for cancer genomics

## Contributors

Qingxiang (Allen) Guo  
Postdoctoral Fellow  
Northwestern University, Feinberg School of Medicine
<qingxiang.guo@northwestern.edu>

## Introduction

A collection of installation and usage guides, focusing on cancer genomics tools. Written by someone who survived one too many 'comprehensive' manuals and lived to tell a simpler tale.

## Tools

## Data transfer

### [fasterq-dump](/contents/fasterq.md)

## Configuration and management

### [Conda and Mamba](/contents/conda.md)

### [Git](/contents/git.md)

### [Homebrew](/contents/homebrew.md)

### [Poetry](/contents/poetry.md)

### [Pre-commit](/contents/pre-commit.md)

### [nox](/contents/nox.md)

### [Singularity](/contents/singularity.md)

## DNA and RNA-seq aligner (splice aware)

### [Minimap2](/contents/minimap2.md)

## Splice unware aligner

### [BWA-MEM](/contents/bwa.md)

### [BWA-MEM2](/contents/bwa2.md)

## RNA-seq aligner (splice aware)

### [STAR](/contents/STAR.md)

## Manipulating and analyzing alignments

### [Samtools](/contents/samtools.md)

### [Picard](/contents/picard.md)

### [mosdepth](/contents/mosdepth.md)

### [BCFtools](/contents/bcftools.md)

### [Qualimap](/contents/qualimap.md)

### [BEDtools](/contents/bedtools.md)

### [RSeQC](/contents/RSeQC.md)

### [featureCounts](/contents/featureCounts.md)

### [Stringtie](/contents/stringtie.md)

### [SQANTI3](/contents/SQANTI3.md)

## Indel calling

### [transindel](/contents/transindel.md)

## Genome structural variation detection (NGS)

### [lumpy](/contents/lumpy.md)

### [Delly](/contents/delly.md)

### [Manta](/contents/manta.md)

### [SVABA](/contents/SVABA.md)

## Genome structural variation detection (long read)

### [SVIM](/contents/SVIM.md)

### [PBSV](/contents/pbsv.md)

### [Sniffles2](contents/sniffles2.md)

### [cuteSV](/contents/cuteSV.md)

### [SVDSS](/contents/SVDSS.md)

### [DeBreak](/contents/DeBreak.md)

## Genome structural variation VCF manipulation and downstream analysis

### [SURVIVOR](/contents/SURVIVOR.md)

### [surpyvor](/contents/surpyvor.md)

### [Truvari](/contents/truvari.md)

### [svanalyzer](/contents/svanalyzer.md)

### [AnnotSV](/contents/AnnotSV.md)

## Genome structural variation simulation

### [VISOR](/contents/VISOR.md)

## Gene fusion analysis - RNA-seq level

### [Arriba](/contents/arriba.md)

### [ctat_lr_fusion](/contents/CTAT_LR_Fusion.md)

## Nanopore ONT analysis

### [guppy](/contents/guppy.md)

### [Dorado](/contents/Dorado.md)

### [Chopper](/contents/Chopper.md)

### [Porechop](/contents/porechop.md)

### [Porechop_ABI](/contents/porechop_ABI.md)

### [pod5](/contents/pod5.md)

### [NanoPlot](/contents/nanoplot.md)

### [PycoQC](/contents/pycoQC.md)

### [Duplex tools](/contents/duplex_tools.md)

### [Cutadapt](/contents/cutadapt.md) 

## Dealing with FastQ/FASTA

### [Seqkit](/contents/Seqkit.md)

### [FastQC](/contents/FastQC.md)

### [BBMap](/contents/bbmap.md)

## modification analysis

### [Modkit](/contents/modkit.md)

# Other small tricks and tips

# Find and load R in Northwestern quest  

You can see which versions of R are available on Quest, and which version is the default, with the command  

```bash
module spider R
```

You can make a particular version of R available to use by typing the full module name with the version included as listed in the output  

```bash
module load R/4.2.0
```

# Batch download SRA files from NCBI

Go to the NCBI SRA site, select the SRA file you need and go to "Run selector" to got batch list in this way.

Use perl one-liner to process the header

```bash
perl -p -i -e 's/(\S+)\n/"$1", /g' list
```

Revise the script and batch download with the python script:

```bash
# _*_ coding=utf-8 _*_
import subprocess

# Insert your SRR accession numbers here
sra_accession_number = ["SRR11951439", "SRR11951443", "SRR11951444", "SRR11951445", "SRR11951446", "SRR11951447", "SRR11951448",
                        "SRR11951449", "SRR11951450", "SRR11951451", "SRR11951452", "SRR11951453", "SRR11951454", "SRR11951455",
                        "SRR11951456", "SRR11951457", "SRR11951458", "SRR11951459", "SRR11951460", "SRR11951461", "SRR11951462",
                        "SRR11951463", "SRR11951464", "SRR11951465", "SRR11951466", "SRR11951467", "SRR11951468", "SRR11951469",
                        "SRR11951470", "SRR11951471", "SRR11951472", "SRR11951473", "SRR11951474", "SRR11951475", "SRR11951476",
                        "SRR11951477", "SRR11951478", "SRR11951479", "SRR11951480", "SRR11951481", "SRR11951482", "SRR11951483",
                        "SRR11951484", "SRR11951485", "SRR11951486", "SRR11951487", "SRR11951488", "SRR11951489", "SRR11951490",
                        "SRR11951491", "SRR11951492", "SRR11951493", "SRR11951494", "SRR11951495", "SRR11951496", "SRR11951497",
                        "SRR11951498", "SRR11951499", "SRR11951500", "SRR11951501", "SRR11951502", "SRR11951503", "SRR11951504",
                        "SRR11951505", "SRR11951506", "SRR11951507", "SRR11951508", "SRR11951509", "SRR11951510", "SRR11951511",
                        "SRR11951512", "SRR11951513", "SRR11951514", "SRR11951515", "SRR11951516", "SRR11951517", "SRR11951519",
                        "SRR11951520", "SRR11951521", "SRR11951522", "SRR11951523", "SRR11951524", "SRR11951525", "SRR11951526",
                        "SRR11951527", "SRR11951528", "SRR11951530", "SRR11951531", "SRR11951532", "SRR11951533", "SRR11951534",
                        "SRR11951535", "SRR11951536", "SRR11951537", "SRR11951542"]

# This will download the .sra files to ~/ncbi/public/sra/ (will create directory if not present)
for sra_id in sra_accession_number:
    print("Current downloading" + sra_id)
    prefetch_cmd = "prefetch " + sra_id + " -O /home/qgn1237/qgn1237/2_raw_data/smooth_seq_95_sc_K562_SMRT" # Note the space between commands
    print("The running command is " + prefetch_cmd)
    subprocess.call(prefetch_cmd, shell=True)

# This will create the fastq file from the downloaded sra file
for sra_id in sra_accession_number:
    print("Creating fastq file for" + sra_id)
    fasterq_dump_cmd = "fasterq-dump /home/qgn1237/qgn1237/2_raw_data/smooth_seq_95_sc_K562_SMRT/" + sra_id + " -O /home/qgn1237/qgn1237/2_raw_data/smooth_seq_95_sc_K562_SMRT/" + sra_id # Use the fstring here
    print("The running command is " + fasterq_dump_cmd)
    subprocess.call(fasterq_dump_cmd, shell=True)
```

# Batch remove a certain type of files from a directory

```bash
python batch_delete_all_sra_files.py
```

# Use Mamba or Conda environment in Slurm Batch sbumit in NU Quest

You have to resource the mamba by yourself

```bash
#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics ## Required: (buyin, short, normal, long, gengpu, genhimem, etc)
#SBATCH --time=03:00:00 ## Required: How long will the job need to run (remember different partitions have restrictions on this parameter)
#SBATCH --nodes=1 ## how many computers/nodes do you need (no default)
#SBATCH --ntasks-per-node=6 ## how many cpus or processors do you need on per computer/node (default value 1)
#SBATCH --mem=35G ## how much RAM do you need per computer/node (this affects your FairShare score so be careful to not ask for more than you need))
#SBATCH --job-name=allen_download ## When you run squeue -u  this is how you can identify the job
#SBATCH --output=output.log ## standard out and standard error goes to this file

# A regular comment in Bash
/home/qgn1237/2_software/mambaforge/bin/mamba init

source ~/.bashrc

mamba activate mamba666

minimap2 -d /home/qgn1237/qgn1237/1_my_database/GRCh38_p13/minimap2_index/GRCh38.p13.genome.mmi /projects/b1171/qgn1237/1_my_database/GRCh38_p13/GRCh38.p13.genome.fa -t 6
```

# Batch create directories from list

```bash
python ./batch_create_directory_from_list.py
```

# Cancel a job in NU Quest

```bash
scancel -u NETID
```

# Kill a job

```bash
scancel 35087 
```

# Solve the error when connecting quest using SSH client: client_global_hostkeys_private_confirm: server gave bad signature for RSA key 0

This might because the cliend delete your key file, and the server can't recognize you.
Run the cmd to regenerate the keys:

```bash
/usr/bin/ssh-keygen -A
```

# Estimate genome coverage from sorted BAM file to exclude files with too low genome coverage

```estimate_bam_genome_coverage.py```

It will produce a coverage_list, exclude single cells that less than 10% coverage, and get a bad list

```second_coloumn_smaller_than.py```

# Calculate the genome depth (total mapped length/total genome base)

```
calculate_genome_depth_fast.py <your BAM file>
```

# How to remain connect with SSH on MAC

```bash
sudo vim /etc/ssh/sshd_config
```

Then edit the content: #ClientAliveInterval 900 #ClientAliveCountMax 10

```bash
sudo launchctl start com.openssh.sshd
```

# Make a alias for your command

```bash
alias sq='squeue | grep "netid"'

alias cdd='cd ../../'

alias ma='mamba activate mamba666'

alias rl='readlink -f'

alias sb='sbatch'

alias lsn='less -SN'

alias vim='lvim'

alias ls='lsd'
```

Add these lines to ~/.bash.rc then source

# Upload seq data directly from Nanopore to NU Quest server

You can treat Nanopore MinION as a server itself. First, login into MinION with wifi, and get the IP address of MinION. Then, at your on PC terminal, type： ssh -x minit@<IP>. The it will ask you for the password. The password is: minit. Then you will login into MinION system.

Navigate the system, cd .. to the higher level of directory and find your data.

Then use the following command to transfer

```bash
rsync -avz /path/to/minion/files qgn1237@quest.northwestern.edu:/home/qgn1237

You can add --append-verify if you want to resume.

# Then it will ask you for the password of northwestern Quest
```

# Use bash to loop make new directory in current directory

```bash
for dir in */; do cd "$dir"; mkdir 10X_depth ; cd ..; done
```

# Generate sbatch cmd list by specifying the input file and the suffix you want, this can sava time and avoid mistake

```bash
# Edit the cmd = f"samtools view -s 0.384 -b {input_file_abs_path} -@ 8 > {input_file_base_name}.{suffix}" in generate_cmd.py to your desired command
./generate_cmd.py ../VCaP.bam 10X.bam 

# Then excute in the directory where you want to store the ouput file
sbatch cmd_list

# You can customize the command as you like

./pbsv2_cmd.py 22Rv1_10X_sort.svsig.gz var.vcf
```

# Add the PATH environment variable

```bash
echo 'PATH=$PATH:/home/qgn1237/1_script' >> ~/.bashrc
source ~/.bashrc
```

# Find the files end with certain suffix in current direcotry

```bash
find_suffix_in_current_directory.py <suffix>
```

# Filter vcf based on quality

```bash
filter_vcf_based_on_quality.py input.vcf 20 > output.vcf
```

# Install and import Dracula theme (my favorite theme) into Iterm2, Chrome Vimium, Chrome

For Iterm2 go to any folder and run command:

```bash
git clone https://github.com/dracula/iterm.git

cd iterm/  # Then go to the profiles panel of Iterm2, select your profiles, --> . colors, and import the  Dracula.itermcolors file, 

# Then select Dracula in color presets done!
```

For Vimium,go into the Vimium addon's preferences, right-click the Vimium icon next to the address bar and select “Manage Extension.” On the Extensions page, paste <https://raw.githubusercontent.com/dracula/vimium/master/vimium-dracula.css> to CSS for Vimium UI

For Chrome, install in <https://chrome.google.com/webstore/detail/dracula-chrome-theme/gfapcejdoghpoidkfodoiiffaaibpaem>

# Copy all the file including hidden files

cp -r ~/OneDrive\ -\ Northwestern\ University/deep_learning_math_theory/. ./

# Use Python script to fulfill the function of Bash "for $dir in ./; do ..."

```bash
bash_cmd_to_python.py
```

# How to configure the local SSH and key to enable rapid login

first edit the file: ~/.ssh/config
Add content:

```bash
Host quest
    HostName quest.northwestern.edu
    User XXX
    Port 22
```

create a SSH key and copy to server

```bash
ssh-keygen
# you will be prompted to enter a passphrase for the key, enter for no password. This is a password for password. Don't use it.

ssh-copy-id quest
# copy your public key to an existing server, it will prompt you for the password of the remote user’s account
# Type in the password (your typing will not be displayed for security purposes) and press ENTER. 
# The utility will connect to the account on the remote host using the password you provided. It will then copy the contents of your 
# ~/.ssh/id_rsa.pub key into a file in the remote account’s home ~/.ssh directory called authorized_keys.

# You should be able to log into the remote host without the remote account’s password
ssh quest
```

# Download from quest to local and upload from local to quest

```bash
rsync -azvhP quest:~/R/ ./local

rsync -azvhP ./Tools-for-Cancer-Genome-Analysis/scripts/plot_VCF_chrom_freq_count.py quest:~/1_script/
```

# Install NeoVim from source code

You need to install new gcc, cmake, and make

```bash
git clone https://github.com/neovim/neovim

cd neovim

mamba install -c conda-forge cmake # Install cmake

make CMAKE_INSTALL_PREFIX=/home/qgn1237/2_software/mambaforge/envs/mamba666/ # Avoid root

make install

# add alias
alias vim='nvim'
```

# Install NeoVim by brew

```bash
brew install neovim
```

# Find a file by name

```bash
find ./ -name "my_software_overlapp*"
```

# Double loop in bash, this can make you more efficient

```bash
for i in 22Rv1 DU145 LNCaP PC3 VCaP; do cd $i; for j in ./*depth; do rl $j/delly/*_filtered.vcf > $j/SURVIVOR/list_vcf; rl $j/manta/results/variants/tumorSV.vcf >> $j/SURVIVOR/list_vcf; done; cd ..; done
```

# Install Cargo (Rust)

```bash
curl https://sh.rustup.rs -sSf | sh
# Type 1 and enter
fish_add_path /Users/qgn1237/.cargo/bin
```

You can use cargo now.

# Usage of Rust Cargo to install rm-improved, which will give you a trash can for deletion

```bash
cargo install rm-improved
# Then you will have 'rip' command
mkdir trash_can
# Remove file
rip t.p 
# Undo the last deletion
rip -u
# Change the graveyard 
export GRAVEYARD=/home/qgn1237/qgn1237/trash_can
# You can also udo specific files
rip -u 
```

# Error: AttributeError: module 'matplotlib' has no attribute 'use'"

```bash
# Uninstall and install matplotlib again
pip uninstall matplotlib
pip install matplotlib
```

# Useful Rust command collection, lsd command, dust command in Rust

The next gen ls command, a more intuitive version of du in rust

```bash
cargo install lsd
cargo install du-dust
# You can then alias the lsd to ls
```

# Useful Rust command collection, fd-find

```bash
cargo install fd-find
```

You can use it like,

```bash
fd vcf  # part of the file name
fd -e vcf # Specific extension
```

# Useful Rust command collection, exa, to check files by tree

```bash
cargo install fd-find
```

You can use it like,

```bash
exa -T
```

# Soft link your file to Onedrive in Windows

```dos
mklink /D "OneDrive path/file_name" "/local/file_name"
```

Then the file will be automatically updated.

# Install Starship prompt in your linux system or mac local system

Starship is the minimal, blazing-fast, and infinitely customizable prompt for any shell! Let's install in quest first.

```bash
curl -O https://starship.rs/install.sh
chmod 755 install.sh
# Install to your own directory
./install.sh -b ~/.local/bin

# If you use bash shell, add the following to the end of ~/.bashrc:
# Before conda init block!
eval "$(starship init bash)"
source ~/.bashrc

# If you use mac Fish, add the following in the center  ~/.config/fish/config.fish:
starship init fish | source 

# Exit and you're done!
```

# Make use of ~/.local/bin

Try to make a bin file if there is not. The /home/qgn1237/.local/bin will by default in your $PATH. So try to link or move the binary file of your different mamba/conda ENV in this directory, so you don't need to shift between ENV now.

For example

```bash
[qgn1237@quser23 ~]$ which vim
alias vim='nvim'
 ~/.local/bin/nvim
```

# Install Fish on MAC and switch from Bash to Fish as default shell

First make sure you have installed Homebrew.

```bash
brew install fish
# Make sure that you've successfully installed Fish
fish --version
# Make Fish as your default shell
# Add /usr/local/bin/fish to your list of shells

# It should be noted that, for the M2 chip mac, homebrew will not install stuff
# in /usr/local/bin/fish , but /opt/homebrew/bin/fish

# Besides, the following command will not work 
sudo echo /usr/local/bin/fish >> /etc/shells
Changing shell for qgn1237.
Password for qgn1237:
chsh: /usr/local/bin/fish: non-standard shell

# This is because, the latter doesn't work because the 
# output redirection >> is (tried to be) applied by the shell 
# before the sudo … is executed, and of course the user shell 
# has no permission to do that.
```

Instead, you should do:

```bash
echo /opt/homebrew/bin/fish | sudo tee -a /etc/shells
chsh -s /opt/homebrew/bin/fish
```

Exit and restart your shell, you'are done.

```bash
# Fish can do many little tricks
# Calculation
math 2/5
0.4
```

# Install Fish shell on a remote cluster that you can not sudo

I've rarely explore the other possibility besides bash, and I don't think comfort zone is a good sign for me. So I decide to switch to Fish, which seems to be more powerful.

```bash
# You can refer to https://yangyangli.top/posts/012-make-a-powerful-ternimal/
# Know your current shell
echo $SHELL
# Go to mamba environment
ma
# Build from source code
wget https://github.com/fish-shell/fish-shell/archive/refs/tags/3.6.0.tar.gz
tar zxf 3.6.0.tar.gz
cd fish-shell-3.6.0/
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=$HOME/.local ..
# This would install fish and its associated files under $HOME/.local. The fish executable would be located in $HOME/.ocal/bin (which you may want to add to your $PATH).
make
make install
```

When you type: fish, it shows:
fish: error while loading shared libraries: libpcre2-32.so.0: cannot open shared object file: No such file or director

This because fish can not find the lib file, and you also don't have the permission
to /usr/local/lib. The solution is to find the path of libpcre2-32.so.0

```bash
❯ find_suffix_in_current_directory.py libpcre2-32.so.0
./2_software/mambaforge/pkgs/pcre2-10.40-hc3806b6_0/lib/libpcre2-32.so.0
./2_software/mambaforge/envs/mamba_py_39/lib/libpcre2-32.so.0
./2_software/mambaforge/envs/mamba666/lib/libpcre2-32.so.0

echo 'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:~/2_software/mambaforge/pkgs/pcre2-10.40-hc3806b6_0/lib/' >> ~/.bashrc
source ~/.bashrc

# Type fish and you're good!
# Now make it as your default shell
vim ~/.config/fish/config.fish
# add 
starship init fish | source
# When you go to fish, the starship will prompt automatically.
```

Next step, let's make you to change to Fish automatically while connect to server.
chsh won't do it without the shell being listed in /etc/shells, and other ways of editing that configuration also require root permissions.

```bash
# Add this to .bashrc
if [[ $- = *i* ]]; then
   exec fish
fi
```

Now you're all good, everytime you login into bash, bashrc will lead you to fish
Then fish will start starship automatically.

# The fish version can not be shown, the git version may be too old

```bash
# Upgrade your git command first
mamba install -c anaconda git
# Uninstall fish
rip /home/qgn1237/.local/etc/fish/
rip /home/qgn1237/.local/share/fish
rip /home/qgn1237/.config/fish
cd /home/qgn1237/.local/bin
rip fish fish_indent fish_key_reader

# Reinstall fish
git clone https://github.com/fish-shell/fish-shell.git
cd fish-shell/
cmake -DCMAKE_INSTALL_PREFIX=$HOME/.local .
make
make install

fish --version
# fish, version 3.6.0-24-g16369a3ab

# You're good now. The problem is because that the git is too old.
# Fish generates the version at build time with the build_tools/git_version_gen.sh script
```

# Permantly add alias (functions) to fish

```fish
alias -s sq='squeue | grep "qgn1237"'

alias -s cdd='cd ../../'

alias -s ma='mamba activate mamba666'

alias -s rl='readlink -f'

alias -s sb='sbatch'

alias -s lsn='less -SN'

alias -s vim='lvim'

alias -s ls='lsd'

alias -s zcp='zellij action close-pane'

# That creates ~/.config/fish/functions/ls.fish which will then be available in any fish session.
```

Or you can edit the fish rc file

```bash
vim ~/.config/fish/config.fish
# Add all the alis to thr fishrc file
```

# Add mamba to Fish shell, so you can use mamba activate mamba666 instead of mamba init first 

You will get error like:

Run 'mamba init' to be able to run mamba activate/deactivate
and start a new shell session. Or use conda to activate/deactivate.

Fish will automatically load your .bashrc file.
But you have to add mamba to fish manually.

```fish
# add the following to ~/.config/fish/config.fish
# Adapted from `function conda` printed by `conda shell.fish hook`, following
# `mamba.sh`.
function mamba --inherit-variable CONDA_EXE
  if [ (count $argv) -lt 1 ]
    $CONDA_EXE
  else
    set -l mamba_exe (dirname $CONDA_EXE)/mamba
    set -l cmd $argv[1]
    set -e argv[1]
    switch $cmd
      case activate deactivate
        eval ($CONDA_EXE shell.fish $cmd $argv)
      case install update upgrade remove uninstall
        $mamba_exe $cmd $argv
        and eval ($CONDA_EXE shell.fish reactivate)
      case '*'
        $mamba_exe $cmd $argv
    end
  end
end 
```

```fish
mamba init fish
# This will add lines to your fish rc file
```

Now you can do:

```fish
mamba activate mamba666
# Backup your rc file
cp ~/.config/fish/config.fish ~/.config/fish/config.fish.bk
```

The idea is to give mamba the power to do mamba activate like conda activate.

# Install Oh-my-fish shell on a remote cluster that you can not sudo

Oh My Fish provides core infrastructure to allow you to install packages which extend or modify the look of your shell. It's fast, extensible and easy to use.

```fish
curl https://raw.githubusercontent.com/oh-my-fish/oh-my-fish/master/bin/install > install
fish install --path=~/.local/share/omf --config=~/.config/omf
# Install aborted: Git version 1.9.5 or greater required; you have 1.8.3.1
mamba activate mamba666 # We have the new git here.

fish install --path=~/.local/share/omf --config=~/.config/omf
# It will succeed this time.
```

# Install fish LOGO function by oh-my-fish

```fish
omf install fish_logo

# Add this to your fish greeting with this function:
function fish_greeting
    fish_logo
end

# Just write it to ~/.config/fish/functions/fish_greeting.fish
# You will see it every time you start a new session.
```

<div align=center>
<img src="img/fish_logo.png">
</div

# Install Oh-My-Fish theme

```fish
# All your available theme
omf theme
omf install bobthefish
```

# Uninstall Oh-My-Fish  

```fish
ma # Go to mamba666 env with higher git version
omf destroy
```

# Install Starship presets

```fish
starship preset pastel-powerline > ~/.config/starship.toml
```

# Set environmental variable in Fish

```fish
set -U fish_user_paths /opt/homebrew/bin $fish_user_paths
```

# A more efficient way to produce commandline and submit jobs in quest (fish shell)

This is a useful script from Tingyou Wang.

```fish
# n means the number of command you want to divide
divide -i cmd_list_332 -n 4 
sub # Write by myself
```

# Use Fish shell to do loop

The fish shell language is more clean and comfortable.

```fish
for dir in (ls -d SRR11563614 SRR11563615 SRR11563616)
    set source_dir "$dir/$dir"
    set target_dir "$dir"
    mv "$source_dir"/*.sra "$target_dir"
    rip "$source_dir"
end
```

# Install Fisher and install functions

A plugin manager for Fish—the friendly interactive shell.

```fish
curl -sL https://git.io/fisher | source && fisher install jorgebucaran/fisher

# Install functions
fisher install laughedelic/fish_logo
```

# Setting vscode as the default editor for text files in Mac Finder

```fish
brew install duti

# To set all text files and code files to vscode instead of xcode use this set of shell commands.

duti -s com.microsoft.VSCode public.json all
duti -s com.microsoft.VSCode public.plain-text all
duti -s com.microsoft.VSCode public.python-script all
duti -s com.microsoft.VSCode public.shell-script all
duti -s com.microsoft.VSCode public.source-code all
duti -s com.microsoft.VSCode public.text all
duti -s com.microsoft.VSCode public.unix-executable all
# this works for files without a filename extension
duti -s com.microsoft.VSCode public.data all

duti -s com.microsoft.VSCode .c all
duti -s com.microsoft.VSCode .cpp all
duti -s com.microsoft.VSCode .cs all
duti -s com.microsoft.VSCode .css all
duti -s com.microsoft.VSCode .go all
duti -s com.microsoft.VSCode .java all
duti -s com.microsoft.VSCode .js all
duti -s com.microsoft.VSCode .sass all
duti -s com.microsoft.VSCode .scss all
duti -s com.microsoft.VSCode .less all
duti -s com.microsoft.VSCode .vue all
duti -s com.microsoft.VSCode .cfg all
duti -s com.microsoft.VSCode .json all
duti -s com.microsoft.VSCode .jsx all
duti -s com.microsoft.VSCode .log all
duti -s com.microsoft.VSCode .lua all
duti -s com.microsoft.VSCode .md all
duti -s com.microsoft.VSCode .php all
duti -s com.microsoft.VSCode .pl all
duti -s com.microsoft.VSCode .py all
duti -s com.microsoft.VSCode .rb all
duti -s com.microsoft.VSCode .ts all
duti -s com.microsoft.VSCode .tsx all
duti -s com.microsoft.VSCode .txt all
duti -s com.microsoft.VSCode .conf all
duti -s com.microsoft.VSCode .yaml all
duti -s com.microsoft.VSCode .yml all
duti -s com.microsoft.VSCode .toml all
```

# Save the ~/.config file to save your shell configuration

```bash
# You can save this whole directory, make your transfer and installation easier
/home/qgn1237/.config
```

# Useful fish command collection - autopair.fish

```fish
# autopair.fish, auto-complete matching pairs in the Fish command line. Type pair () for you.
fisher install jorgebucaran/autopair.fish
```

# Useful fish command collection - z.fish

```fish
z tracks the directories you visit. With a combination of frequency and recency,
it enables you to jump to the directory in mind.

# Install
fisher install jethrokuan/z
z foo: Goes to directory best matching foo.
zo foo: Opens file manager of directory best matching foo.
z -t foo: Goes to most recent directory matching foo.
z -x: Removes the current directory from $Z_DATA.
```

# Calculate the allele frequency from the AD and SAC info field in a structural variants vcf file

chr1    10837   svim.INS.4      N       <INS>   8       PASS    SVTYPE=INS;END=10837;SVLEN=79;SUPPORT=7;STD_SPAN=5.8;STD_POS=40.06      GT:DP:AD        0/1:30:23,7

 It is possible to calculate the allele frequency from the AD info in this specific structural variants VCF file. The AD field lists the depth of coverage for each allele at a particular variant site, in this case, the AD field has two values separated by a comma (23, 7). To calculate the frequency, you can divide the number of reads supporting the alternate allele (7) by the total number of reads (23 + 7 = 30). So the allele frequency would be 7/30 = 0.23 or 23%.

DP4 is another form of SVC

DP4=0,0,14,12;

DP4 reports the number of reads covering the position with the reference allele mapped to forward and reverse strands, followed by alternate allele mapped to forward and reverse strands.

Allele frequency here would be 0 (0/26) for reference allele, and 1.0 (26/26) for alternate allele.

# xcrun: error: invalid active developer path (/Library/Developer/CommandLineTools), missing xcrun at: /Library/Developer/CommandLineTools/usr/bin/xcrun, you can't use your git after upgrading macOS

```bash
xcode-select --install
```

This will pop a dialogue box, Select "Install", and it will download and install the Command Line Tools package and fix the problem.

# Update Yang lab Library

Edit the library.yaml file first, add the book name and STAR

Then

```
make cover
git add ..
git commit -n "n"
git push origin main
```

# Configure your NeoVim

Neovim (and many other applications) doesn't touch your ~/.config/ directory per default. You just need to create the nvim directory

It's different from Vim, you need to:

```bash
mkdir -p ~/.config/nvim

nvim ~/.config/nvim/init.vim # Edit it

# For example, add: set mouse=v
```

# Can't Copy to Clipboard from neovim, what's wrong?

Becasue nvim sets your vim into visual mode whenever you select something with the mouse. And for some mad reason one is not allowed to copy when in visual mode. You can get around it by holding down shift when selecting text not to go into visual mode allowing you to use the copy menu.

Simply edit the ~/.config/nvim/init.vim and add:

set mouse=v

Done!

# Github push error: fatal: Unable to create '/Users/qgn1237/Library/CloudStorage/OneDrive-NorthwesternUniversity/github_project/Computational-Medicine-and-Bioinformatics-Terminology-Database/.git/index.lock': File exists

The error message suggests that there is a lock file present, which typically means that another Git process is running or a previous Git process did not exit cleanly. To resolve the issue, follow these steps:

Make sure there are no other active Git processes running. Close any relevant terminal windows or Git applications that might be accessing the repository.

Remove the lock file manually. In your terminal, navigate to the repository's root directory and run the following command:

```bash
rm -f .git/index.lock
```

# Extract all the text from a powerpoint file

```bash
pip install python-pptx

./extract_text_ppt.py 【20230413】qingxiang_guo_basic_seminar.pptx

# You can find it in the script directory
```

# Filter the structural variants VCF results

1. Filter the SVIM results

```bash
# A SVIM result looks like this:
GL000008.2  1524  svim.DEL.7  AGAATGGGATGGAATGGAATTGAATGATGTGGAGTGGAGTCGGGTGGAGTGGA A 2 PASS  SVTYPE=DEL;END=1576;SVLEN=-52;SUPPORT=2;STD_SPAN=4.24;STD_POS=14.85 GT:DP:AD  ./.:.:.,.

# Filtering the variant calls by QUAL, supporting reads, and PASS, we cannot make a general statement about suitable score cutoffs. # For high-coverage datasets (>40x), we would recommend a threshold of 10-15. For low-coverage datasets, the threshold should be lower.

# A good approach could also be to select a cutoff that returns the expected or desired number of calls. A recent study by the Human Genome Structural Variation Consortium (HGSVC), for instance, estimated the average number of deletions and insertions per individual to be 9,488 and 15,337, respectively (Chaisson et al., 2019)

bcftools filter -i 'QUAL >= 5 && FILTER == "PASS" && INFO/SUPPORT >= 2' input.vcf > output.vcf

This is a strict criteria that can filter 400,000 SVs to 2,000

So you should make trial & error
```

2. Filter the PBSV results

```bash
chr1    669597  pbsv.INS.DUP.0  A       <DUP>   .       PASS    SVTYPE=DUP;END=669679;SVLEN=82  GT:AD:DP:SAC    1/1:1,4:5:1,0,3,1
# PBSV will not have the SUPPORT reads option

bcftools filter -i 'FILTER == "PASS" && FORMAT/AD[0:1] >= 2' input.vcf > output.vcf
# When filtering SVs, the variants with high support can be filtered based on the threshold value of AD [1].
# The [0:1] means the 1 subfield of the 0 sample
```

3. Filter the Sniffles2 results
  
```bash
# Sniffles has already an automatic filter in the output.
# If you want to filter further, the easiest is to ignore the IMPRECISE marked calls 
# Another thing would be to filter for the number of read supports

chr1    710590  Sniffles2.INS.14S0      N       AAGAACTGCCTGCCGGGCGCGTGTGGCTCACGCCTTGTAATTCCCAGCACTTTGGGAGGCCGCAGGCGGGCCGGATCACGAGCGTCAGGAGATCGAGACCATCCCGGCTAAAACGGCTTTGAAAAAACCCCGTCTCTACTAAAAATTACAAAAAATTAGCCCCGTAGTGGCGGGCGCGTTTAGTCCCAGCTACTTTGCGGAGGCTGAGGCAGGAAGAATGGGCGTGAACCCGCGGAGGTGGAGCTTTGCAGTGAGCCTAAGATCCCACTCACTCCAGCCTGCGGCGACAGAGCCAGACTCGTCTCAAAAAAAAAAAAAAAAAAAAAAAAAAAA   31      PASS    PRECISE;SVTYPE=INS;SVLEN=333;END=710590;SUPPORT=3;COVERAGE=2,3,3,3,3;STRAND=+-;AF=1.000;STDEV_LEN=18.502;STDEV_POS=6.083;SUPPORT_LONG=0 GT:GQ:DR:DV     1/1:8:0:3

# DV：Depth of variant-supporting bases

bcftools view -i 'INFO/SUPPORT >= 2 && INFO/IMPRECISE != 1' intput.vcf > output.vcf
```

4. Filter the CuteSV results

```bash
chrY 56833655 cuteSV.DEL.48873 TCCTATTCCATTCCTC T . PASS PRECISE;SVTYPE=DEL;SVLEN=-15;END=56833670;CIPOS=-1,1;CILEN=-0,0;RE=3;RNAMES=NULL;STRAND=+- GT:DR:DV:PL:GQ ./.:.:3:.,.,.:.

bcftools view -i 'FILTER == "PASS" && INFO/PRECISE == 1 && INFO/RE >= 2' input.vcf > output.vcf
```

5. Filter the DeBreak results

```bash
chrY 2464326 DB4000 TTCCCTAGCAATCCGGCCAAGGGCCGCTGATGTGCACACACTGAAGACGTTCCCTAAGTGTGTGGCTAAGGGACTGCTACCATATACACACTGAAGATGTTCCCTAAGAATGTGGGTAAGGGACCGCCGCCATGTTCGCACTGAAGACG N . PASS CHR2=chrY;SVLEN=149;SUPPREAD=4;MAPQ=0;SVMETHOD=DeBreak;PRECISE;SVTYPE=DEL GT 0/1

bcftools view -i 'FILTER == "PASS" && INFO/PRECISE == 1 && INFO/SUPPREAD >= 2' input.vcf > output.vcf
```

6. Filter the SVDSS results

```bash
chrX 147512934 DEL_chrX:147512934-147512961_28 TTGCAGTACAATACACATTGTATTACAC T . PASS VARTYPE=SV;SVTYPE=DEL;SVLEN=-28;END=147512961;WEIGHT=2;COV=4;AS=358;NV=1;CIGAR=33=28D46= GT 0/1

# You can filter the reported SVs by passing the --min-sv-length and --min-cluster-weight options. These options control the minimum length and minimum number of supporting superstrings for the reported SVs. Higher values for --min-cluster-weight will increase precision at the cost of reducing recall. For a diploid 30x coverage sample, --min-cluster-weight 2 produced the best results in our experiments. For a haploid 30x sample, instead, --min-cluster-weight 4 produced the best results.

bcftools view -i 'FILTER == "PASS" && INFO/NV >= 2' input.vcf > output.vcf
```

# Convert BED to VCF coordinate

If BED start and end = 0  8, VCF start =1, end = 8
BED is 0-based coordinate system, VCF is 1-based

# Filter a VCF file based on the SVLEN field

```bash
filter_vcf_based_on_length.py -i input.vcf -o length.vcf -l 50
```

# How to Restore the MK1C Device from an Update Error

This guide aims to assist users facing a specific issue with the MinION Mk1C device. Users may experience problems after attempting to update the device's software version, causing the device to blackout, display a blue screen with various icons, or show an "Unknown" installed version.

<b>Symptoms:</b>

1. After clicking on "Install Update" for a software version upgrade, the device blacks out for an extended period (e.g., over 16 hours).
2. The device occasionally displays a blue screen with various icons, including menu icons, a mail icon, and the test flow cell icon. Despite prompts to swipe up to unlock, the device remains on the blue screen.
3. After rebooting, the interface changes and the "MinKNOW Installed version" now shows "Unknown".
4. In the Host Settings, the Network Settings option is missing, preventing the user from confirming the network connection status or obtaining the device's IP address.
5. When attempting to run the “ssh minit@mc-XXXX” command, the system reports an error: “ssh: Could not resolve hostname mc-XXXX: nodename nor servname provided, or not known”.
6. Attempting to restore using a previous SD card (with an earlier software version) does not resolve the issue.

<b>Resolution:</b>

<b>Step 1: Preparing a USB drive</b>

1. Format a USB key as FAT or exFAT.
2. Create a directory on the USB called `minit_config`.
3. In the `minit_config` directory, create two plain files (without file extensions): `ssh` and `access-point`. These are command files that will enable SSH and put the device into hotspot mode.

<b>Step 2: Enabling SSH and hotspot mode on the device</b>

1. Plug the USB drive into the powered-on Mk1C device. The white LEDs on the screen should stop scrolling and all flash red twice, indicating the script's activation.
2. The device will execute the commands in the `minit_config` folder and capture some system configuration information. Once complete, the lights will flash green, indicating the process is complete and the drive can be safely removed.

<b>Step 3: Connecting to the device</b>

1. Reboot the device.
2. Connect to the device's hotspot using another device (e.g., a computer). The hotspot should have the same name as the Mk1C device (e.g., MC-1XXXX). The password is: WarmButterflyWings98

<b>Step 4: Accessing the device via SSH</b>

1. SSH into the device using the following command: `ssh minit@10.42.0.1`
Note: If you have connected with other Mk1c deivce before, clean the key by :

`ssh-keygen -R 10.42.0.1` and ssh again.

2. Run the following command to get the device's IP address while connected to Ethernet: `ifconfig`. The IP address is listed for `eth0`.

<b>Step 5: Updating the device</b>

If the device has an internet connection, run the following commands to update it:

```
sudo apt clean
sudo apt --fix-broken install
sudo shutdown -r now
```

<b>Note:</b>

This guide is based on a specific issue and might not apply to all problems encountered with the MinION Mk1C device. For different issues, please contact the Nanopore official technical support or check the Nanopore community forums for assistance.

# How to Update the MinION Mk1C Software via SSH

While the MinION Mk1C software can generally be updated via the user interface, some scenarios require manual software updates. If you prefer to update through the command line or need feedback during updates, this guide is for you.

Please always refer to the release notes for any version-specific commands that may be supplied.

<b>Prerequisites:</b>

1. Ensure that the MinION Mk1C is connected to a stable internet connection.
2. Have SSH access to the MinION Mk1C.

<b>Steps to Update Only the MinKNOW Software:</b>

1. Connect to the MinION Mk1C using SSH.
2. Open a terminal window.
3. Run the following commands one at a time to update the list of available software:

```
sudo apt autoclean
sudo apt clean
sudo apt autoremove
sudo apt-key adv --fetch-keys https://cdn.oxfordnanoportal.com/apt/ont-repo.pub
sudo apt update
sudo apt install ont-mk1c-release
```

These commands should complete successfully, without any errors or warnings.

4. If you encounter any error or warning message, ensure that your MinION Mk1C is connected to the Internet. You can diagnose network connection problems by running the following command:

```
sudo apt --fix-broken install
sudo apt update
sudo apt install ont-mk1c-release
```

5. Once the commands have completed successfully without errors, reboot the device using the command:

```
sudo reboot
```

<b>Steps to Update the Operating System and MinKNOW Software Simultaneously:</b>

If you want to update both the operating system and MinKNOW software at the same time, follow these steps:

1. Connect to the MinION Mk1C using SSH.
2. Open a terminal window.
3. Run the following commands one at a time:

```
export DEBIAN_FRONTEND=noninteractive 
sudo apt upgrade
```

4. Reboot the device using the command:

```
sudo reboot
```

Note: This will update all packages to their latest versions, so it may take longer than the MinKNOW-only update.

Remember to check the 'Current version' under software version number after the updates to ensure that your software has been updated successfully.

# Install LunarVim in Fish Shell

## Prerequisites

Before starting, ensure that you have the following installed:

- Latest version of Neovim (v0.9.0+)
- Git
- Make
- Pip
- Python
- NPM
- Node.js
- Cargo

## Installation Steps

1. Enter your mamba environment:

   ```bash
   ma
   ```

2. Execute the LunarVim installation script using bash:

   ```bash
   env LV_BRANCH='release-1.3/neovim-0.9' bash -c "bash <(curl -s https://raw.githubusercontent.com/LunarVim/LunarVim/release-1.3/neovim-0.9/utils/installer/install.sh)"
   ```

   During this step, you may encounter an error related to the compilation of `nvim-treesitter[python]`. This occurs when the `cc1plus` compiler is not found.

   if you are neovim with > 10 version,use:

   ```bash
   bash <(curl -s https://raw.githubusercontent.com/lunarvim/lunarvim/master/utils/installer/install.sh) 
   ```

3. To resolve this issue, you can install a C++ compiler within your Mamba environment using the following command:

   ```bash
   mamba install -c conda-forge cxx-compiler
   ```

4. Afterwards, set the `CXX` environment variable to point to your environment's `g++`:

   ```bash
   export CXX=$CONDA_PREFIX/bin/g++
   ```

```bash
set -U fish_user_paths $fish_user_paths ~/.local/bin/ .  # Set the environmental variable
make sure to manually run ':Lazy sync' when starting lvim for the first time.
```

# Installing Plugins in LunarVim

This guide shows how to install and use the Dracula theme in LunarVim. Similar steps can be followed for other plugins.

## Step 1: Open the LunarVim Configuration File

First, you need to open your LunarVim configuration file. You can do this with the following command:

```bash
nvim ~/.config/lvim/config.lua
```

## Step 2: Add the Plugin Information

In the configuration file, locate the `lvim.plugins` section. Here you add the plugin information. For example, to add the Dracula theme, you can use:

```lua
lvim.plugins = {
    {
        "Mofiqul/dracula.nvim",
    },
    -- more plugins can go here
}
```

## Step 3: Set the Theme

Next, set your chosen theme. For Dracula, you can use:

```lua
lvim.colorscheme = "dracula"
```

## Step 4: Save and Exit

Save and exit the file.

# Install and using Zellij to manage terminal windows

```bash
cargo install --locked zellij

zellij setup --dump-config > ~/.config/zellij/config.kdl # Create config file

zellij options --theme dracula # change theme

ctrl + O and d to leave the session, zellij a to attach a session, zellij -s to create a new session

zellij action close-pane # Close the current pane.
```

# Accelerating Your Command Line Operations using Zellij and SLURM (original by yangyang Li, https://yangyangli.top/)

This tutorial will guide you through the steps to expedite your command line operations on a compute node from a login node, using `zellij` to maintain an interactive `srun` session in a SLURM-managed cluster. This process can greatly improve your command execution speed, especially if your login node is heavily loaded.

## Prerequisites

- Basic understanding of command line interface (CLI)
- Access to a SLURM-managed cluster
- `zellij` and `srun` command line tools installed

## Steps

### Step 1: Start a Zellij Session

To begin, open your command line interface and start a new Zellij session by typing:

```bash
zellij -s q
```

The `-s q` option gives a name (in this case, "q") to your Zellij session.

### Step 2: Check Your Current Node

You can verify the node you're currently on by executing the `hostname` command:

```bash
hostname
```

This command will return the name of the node you're currently logged into. At this point, you should be on the login node.

### Step 3: Start an Interactive Srun Session

Next, you can start an interactive `srun` session with specified resources. This command will move you from the login node to a compute node, enhancing your command execution speed:

```bash
srun -n 12 -p b1171 --account=b1171 -t 10:00:00 --mem 16g --pty bash
```

In this example, `-n 12` specifies the number of CPUs, `-p b1171` specifies the partition, `--account=b1171` specifies the account, `-t 10:00:00` specifies the time limit, and `--mem 16g` specifies the amount of memory.

### Step 4: Verify the Compute Node

Once the srun session starts, you can verify that you're on a compute node by running the `hostname` command again.

### Step 5: Exit the Zellij Session

When you're done with your tasks on the compute node, you can detach from the Zellij session without ending the `srun` session by using the shortcut `Ctrl + o`, followed by `d`.

## Conclusion

This tutorial walked you through the process of using Zellij and srun to expedite your command line operations on a compute node. Remember, when you exit Zellij, any tasks running within the `srun` session will continue to run until they're completed or until the time limit is reached. Thus, you can significantly improve your productivity by leveraging the computational power of the compute nodes.

Of course! Here's a markdown-formatted guide to solve the problem you described:

---

# Solving Name Resolution Issue in Windows WSL Ubuntu

If you're experiencing issues with SSH and name resolution within the Windows Subsystem for Linux (WSL) Ubuntu system, the following steps can help you address it.

## Symptoms

- Failed SSH attempts due to "failure in name resolution."
- Using `ping` results in an error, e.g., `ping: quest.northwestern.edu temporary failure in name resolution`.
- Accessing websites via a browser on Windows works without issues.

## Solution

### 1. Configuring WSL to Stop Automatically Generating `resolv.conf`

WSL can automatically generate the `resolv.conf` file, which may cause DNS resolution problems. To prevent WSL from auto-generating this file:

1. Open WSL terminal.
2. Create or edit the `wsl.conf` file by typing:
   ```
   sudo vim /etc/wsl.conf
   ```
3. Add the following content to the file:
   ```
   [network]
   generateResolvConf = false
   ```
4. Save the file and exit.

### 2. Manually Setting Up DNS Configuration

You'll need to manually set up the DNS configuration in the `resolv.conf` file:

1. Still within the WSL terminal, create a new `resolv.conf` file:
   ```
   sudo vim /etc/resolv.conf
   ```

2. Add the following content:
```
   # This file was automatically generated by WSL. To stop automatic generation of this file, 
   # add the following entry to /etc/wsl.conf:
   # [network]
   # generateResolvConf = false
   # nameserver 172.23.0.1
   nameserver 8.8.8.8
   ```
3. Save the file and exit.


### 3. Test the Solution

After performing the above steps, test the solution:

1. In the WSL terminal, try pinging a known website:
   ```
   ping google.com
   ```

# Precision-Recall Diagrams using `precision-recall-diagrams.py`

This guide provides a quick walkthrough on how to visualize the trade-off between precision and recall using the provided Python script.

## How to Use

1. **Prepare Your Data**:
   
   The script uses predefined data. If you want to use your own data, modify the `prs` array inside the script:
   
   ```python
   prs = np.array([
       [precision1, recall1],
       [precision2, recall2],
       ...
   ])
   labels = ["Label 1", "Label 2", ...]
   ```

2. **Run the Script**:

   Navigate to the script's directory and execute:

   ```bash
   python precision-recall-diagrams.py
   ```

# Understanding Poetry and Pip in Python Package Management**

## Poetry vs. Pip in Python

**Pip Install**
- **Tool**: Standard package installer for Python.
- **Use Case**: Used by end-users for installing Python packages from PyPI or other sources.
- **Focus**: Manages the installation of individual packages without handling project-wide dependencies.

**Poetry Install**
- **Tool**: Modern Python dependency management and packaging tool.
- **Use Case**: Suited for developers managing dependencies of entire projects.
- **Functionality**: Installs all dependencies based on `pyproject.toml` and may create virtual environments for isolation.
- **Advantages**: Manages project dependencies, ensuring consistent development environments.

## Command Line Tools Installation
- Both `pip` and `poetry` install command line tools of packages to `Scripts` (Windows) or `bin` (Unix-like systems) in the Python environment.
- These directories are typically included in the system's PATH, allowing command line recognition and execution of these tools.

## Entry Points
- **Mechanism**: Allows packages to specify one or more command line script entry points.
- **Implementation**: Scripts created during package installation, placed in executable directories.
- **Configuration**:
  - `setup.py` for pip: `entry_points={'console_scripts': ['mycommand = mypackage.mymodule:main_func']}`
  - `pyproject.toml` for poetry: `[tool.poetry.scripts] mycommand = 'mypackage.mymodule:main_func'`
- **Function**: These scripts are lightweight wrappers calling functions or classes in the package.

## Summary
- `poetry install` enhances Python's environment by managing and installing a project's packages and dependencies, updating module paths (`sys.path`), and setting up command line tools. This process ensures that packages and their dependencies are accessible and functional in the specified Python environment.

# Understanding Poetry, Python Environment, and Entry Points

### Poetry Installation Process
- **Dependency Resolution**: Reads `pyproject.toml` to determine required packages and versions.
- **Virtual Environment**: Installs packages in the active or a new virtual environment.
- **Package Installation**: Downloads and installs packages from PyPI or other sources into the Python environment.
- **Module Path Update**: Adds installed packages' paths to Python's `sys.path`.
- **Entry Point Setup**: Creates scripts for command line tools as defined in `pyproject.toml`.
- **Dependency Locking**: Updates `poetry.lock` to ensure consistent versions across environments.

### Python Environment
- **Python Path (`sys.path`)**: A list where Python looks for modules to import.
- **Site-packages Directory**: Common location for installed packages (e.g., `/Users/qgn1237/mambaforge/envs/test_octopusv/lib/python3.9/site-packages`).

### Entry Points and Package Discovery
- **Executable Scripts**: Located in `bin` or `Scripts` directory of the Python environment.
- **Example Script Path**: `/Users/qgn1237/mambaforge/envs/test_octopusv/bin/octopusv`.
- **Import Logic in Scripts**: Can successfully import modules (e.g., `from octopusv.cli.cli import app`) due to their presence in the `site-packages` directory.

### Locating Installed Packages
- **`octopusv` Package**: Typically found in the `site-packages` directory.
- **Metadata Directory (`octopusv-0.1.0.dist-info`)**: Contains metadata like `direct_url.json`, `entry_points.txt`, `INSTALLER`, `METADATA`, and `RECORD`.

# Automated Software Testing with Nox, Pytest, and GitHub Actions

This tutorial explains how to automate your software testing process using Nox, Pytest, and GitHub Actions.

<div align=center>
<img src="img/action1.png">
</div>

## Step 1: Configure Nox

Create a `noxfile.py` with sessions for different Python versions:

```python
import nox

@nox.session(python="3.9")
def tests_39(session):
    session.run("poetry", "install", external=True)
    session.run("pytest", "--cov", "--cov-report=html")

@nox.session(python="3.10")
def tests_310(session):
    session.run("poetry", "install", external=True)
    session.run("pytest", "--cov", "--cov-report=html")
```

## Step 2: Write Pytest Test Cases

Create a `tests` folder with your test cases. Here's an example test script:

```python
import subprocess
import filecmp
import pytest
import os

# The test dir
TEST_DATA_DIR = "tests/data/VCF_for_testing_correct"
OUTPUT_DIR = "tests/output"
os.makedirs(OUTPUT_DIR, exist_ok=True)


def run_octopusv(input_file, output_file):
    """
    Run octopusv.
    """
    cmd = f"octopusv correct {input_file} {output_file}"
    subprocess.run(cmd, shell=True, check=True)


def test_mate_bnd_independent():
    input_vcf = os.path.join(TEST_DATA_DIR, "test_mate_bnd_independent.vcf")
    output_vcf = os.path.join(OUTPUT_DIR, "output_mate_bnd_independent.vcf")
    expected_output_vcf = os.path.join(TEST_DATA_DIR, "example_mate_bnd_independent.vcf")

    run_octopusv(input_vcf, output_vcf)

    # Compare output and expected file
    assert filecmp.cmp(output_vcf, expected_output_vcf)

```

## Step 3: Setup GitHub Actions

Create a workflow file in `.github/workflows/tests.yml`. This file defines the CI pipeline.

```yaml
name: Tests

on:
  push:
    branches:
      - main
    paths:
      - "**.py"
      - "**/Dockerfile"
      - ".github/workflows/*.yml"
      - poetry.lock
      - pyproject.toml

  pull_request:
    branches:
      - main
    paths:
      - "**.py"
      - "**/Dockerfile"
      - ".github/workflows/*.yml"
      - poetry.lock
      - pyproject.toml

jobs:
  tests:
    name: ${{ matrix.session }} ${{ matrix.python }} / ${{ matrix.os }}
    runs-on: ${{ matrix.os }}
    strategy:
      fail-fast: false
      matrix:
        include:
          - { python: "3.10", os: "ubuntu-latest", session: "pre-commit" }
          - { python: "3.10", os: "ubuntu-latest", session: "safety" }
          - { python: "3.10", os: "ubuntu-latest", session: "mypy" }
          - { python: "3.9", os: "ubuntu-latest", session: "mypy" }
          - { python: "3.10", os: "ubuntu-latest", session: "tests" }
          - { python: "3.9", os: "ubuntu-latest", session: "tests" }
          - { python: "3.10", os: "macos-latest", session: "tests" }
          - { python: "3.9", os: "macos-latest", session: "tests" }
          - { python: "3.10", os: "ubuntu-latest", session: "docs-build" }

    env:
      NOXSESSION: ${{ matrix.session }}
      FORCE_COLOR: "1"
      PRE_COMMIT_COLOR: "always"

    steps:
      - name: Check out the repository
        uses: actions/checkout@v4
        with:
          fetch-depth: 1

      - name: Cache pip packages
        uses: actions/cache@v3
        with:
          path: ~/.cache/pip
          key: ${{ runner.os }}-pip-${{ matrix.python }}-${{ hashFiles('**/requirements.txt') }}
          restore-keys: |
            ${{ runner.os }}-pip-${{ matrix.python }}-
            ${{ runner.os }}-pip-

      - name: Cache Poetry virtual environment
        uses: actions/cache@v3
        with:
          path: .venv
          key: ${{ runner.os }}-poetry-${{ matrix.python }}-${{ hashFiles('pyproject.toml') }}
          restore-keys: |
            ${{ runner.os }}-poetry-${{ matrix.python }}-
            ${{ runner.os }}-poetry-

      - name: Set up Python ${{ matrix.python }}
        uses: actions/setup-python@v4
        with:
          python-version: ${{ matrix.python }}

      - name: Upgrade pip
        run: |
          pip install --constraint=.github/workflows/constraints.txt pip
          pip --version

      - name: Upgrade pip in virtual environments
        shell: python
        run: |
          import os
          import pip

          with open(os.environ["GITHUB_ENV"], mode="a") as io:
              print(f"VIRTUALENV_PIP={pip.__version__}", file=io)

      - name: Install Dependencies
        if: matrix.os == 'ubuntu-latest'
        run: |
          sudo apt-get update
          sudo apt-get install -q --no-install-recommends --no-install-suggests libhts-dev  libssl-dev

      - name: Cache macOS dependencies
        if: matrix.os == 'macos-latest'
        uses: actions/cache@v3
        with:
          path: ~/brew-cache
          key: macos-brew-${{ hashFiles('**/Dockerfile') }}
          restore-keys: |
            macos-brew-

      - name: Configure Homebrew
        if: matrix.os == 'macos-latest'
        uses: Homebrew/actions/setup-homebrew@master

      - name: Install Dependencies
        if: matrix.os == 'macos-latest'
        env:
          HOMEBREW_NO_AUTO_UPDATE: 1
          HOMEBREW_NO_INSTALL_CLEANUP: 1 # Do not run brew cleanup automatically.
          HOMEBREW_NO_INSTALLED_DEPENDENTS_CHECK: 1 # Do not automatically update packages.

        run: |
          brew install htslib
          brew install openssl@3
          ln -sf  $(brew --prefix openssl)/include/openssl /usr/local/include/openssl
          ln -sf  $(brew --prefix openssl)/lib/*a /usr/local/lib/
          ln -sf  $(brew --prefix openssl)/lib/*dylib /usr/local/lib/

      - name: Install Poetry
        run: |
          pipx install --pip-args=--constraint=.github/workflows/constraints.txt poetry
          poetry --version

      - name: debug
        if: matrix.os == 'macos-latest'
        run: |
          which python

      - name: Install Nox
        run: |
          pipx install --pip-args=--constraint=.github/workflows/constraints.txt nox
          pipx inject --pip-args=--constraint=.github/workflows/constraints.txt nox nox-poetry
          nox --version

      - name: Compute pre-commit cache key
        if: matrix.session == 'pre-commit'
        id: pre-commit-cache
        shell: python
        run: |
          import hashlib
          import sys

          python = "py{}.{}".format(*sys.version_info[:2])
          payload = sys.version.encode() + sys.executable.encode()
          digest = hashlib.sha256(payload).hexdigest()
          result = "${{ runner.os }}-{}-{}-pre-commit".format(python, digest[:8])

          print("::set-output name=result::{}".format(result))

      - name: Restore pre-commit cache
        uses: actions/cache@v3
        if: matrix.session == 'pre-commit'
        with:
          path: ~/.cache/pre-commit
          key: ${{ steps.pre-commit-cache.outputs.result }}-${{ hashFiles('.pre-commit-config.yaml') }}
          restore-keys: |
            ${{ steps.pre-commit-cache.outputs.result }}-

      - name: Run Nox for Python 3.9 tests session
        if: matrix.session == 'tests' && matrix.python == '3.9'
        run: nox --session tests_39 --force-color

      - name: Run Nox for Python 3.10 tests session
        if: matrix.session == 'tests' && matrix.python == '3.10'
        run: nox --session tests_310 --force-color

      - name: Upload coverage data
        if: always() && matrix.session == 'tests'
        uses: "actions/upload-artifact@v3"
        with:
          name: coverage-data-${{ matrix.python }}-${{ matrix.os }}
          path: "htmlcov/index.html"

      - name: Upload documentation
        if: matrix.session == 'docs-build'
        uses: actions/upload-artifact@v3
        with:
          name: docs
          path: docs/_build
```

## How It Works

1. **Triggering GitHub Actions**: 
   - Defined in `.github/workflows/tests.yml`.
   - Triggered on push or pull request to the main branch.
   - Monitors changes in Python files, Dockerfiles, workflow files, and dependency lock files.

2. **Running Nox Sessions via GitHub Actions**:
   - GitHub Actions uses the `nox` command, as specified in `noxfile.py`.
   - Runs different Python versions (3.9, 3.10) using separate Nox sessions (`tests_39`, `tests_310`).

3. **Nox Session Execution**:
   - Installs dependencies via `poetry`.
   - Executes Pytest with coverage options: `session.run("pytest", "--cov", "--cov-report=html")`.

4. **Pytest Testing Process**:
   - Located in the `/tests` directory.
   - Test scripts like `test_correct_end_to_end.py` are run.
   - Tests compare the output of `octopusv` commands with expected results using `filecmp.cmp`.

5. **Coverage Report Generation and Upload**:
   - Pytest generates an HTML coverage report.
   - GitHub Actions uploads this report as an artifact for review.

This setup forms a comprehensive CI/CD pipeline, ensuring that each code change is automatically tested for reliability and coverage.

Absolutely, your understanding is correct. Here's a simple Markdown tutorial in English to guide you through the process of working with Overleaf and GitHub without direct integration:

# Working with Overleaf and GitHub

This guide will help you maintain a LaTeX project on Overleaf and synchronize it with a GitHub repository.

## Initial Setup

- **Create a repository on GitHub** for your LaTeX project.
- **Clone the repository to your local machine** using `git clone <repository-url>`.
- **Create your LaTeX project on Overleaf** and begin your writing.

## Regular Workflow

### Editing on Overleaf

1. **Edit your document on Overleaf** as usual.
2. **Download the project** from Overleaf by clicking on `Menu` and selecting `Source` under `Download`.

### Synchronizing with GitHub

3. **Extract the downloaded ZIP file** to your local repository folder.
4. Open your terminal or command prompt and **navigate to your repository folder**.
5. Run the following commands to update your GitHub repository:
   ```sh
   git status                  # Check the changed files
   git add .                   # Add all changes to git
   git commit -m "Description of the changes" # Commit changes
   git push                    # Push changes to GitHub
   ```

### Continuing Work

6. Before you start a new session, **pull the latest version from GitHub**:
   ```sh
   git pull                    # Update your local repository
   ```
7. **Upload the updated files to Overleaf**:
   - Go to your Overleaf project.
   - Click on `Menu` and under `Upload` select `ZIP` to upload your local repository files.

### Collaboration

- Inform your collaborators to **pull from GitHub** before they start editing.
- After they have made changes, they should also **push their updates to GitHub**.
- Always **communicate with your team** when you push or pull changes to avoid conflicts.

## Best Practices

- **Commit frequently** with clear, descriptive messages.
- **Pull from GitHub** before you start working to get the latest version.
- **Push to GitHub** after each significant change or at the end of your working session.
- Use `.gitignore` to exclude files that don't need to be version controlled (e.g., `.aux`, `.log` files).

By following these steps, you can ensure that your LaTeX project is well-maintained and up-to-date across both Overleaf and GitHub.

Remember to replace `<repository-url>` with the actual URL of your GitHub repository. This guide assumes a basic knowledge of git commands and the functionality of Overleaf. If you or your collaborators are not familiar with git, you may need a more detailed tutorial on git usage.

# How to make sure the two files are the same

```bash
cmp --silent sim.srt.bam visor.ba || echo "files are different"
```

# Understanding the NODES(A/I/O/T) Column in `sinfo --summarize`, check nodes availability in SLURM

The `NODES(A/I/O/T)` column in the output of the `sinfo --summarize` command provides a quick summary of the state of the nodes in each partition of a Slurm-managed cluster. Here's what each component of this column represents:

- **A (Active nodes)**: The number of nodes that are currently running jobs. These nodes are fully occupied and are executing tasks.
  
- **I (Idle nodes)**: The number of nodes that are idle and available to accept new jobs. These nodes are ready to be utilized but are not currently in use.

- **O (Offline nodes)**: The number of nodes that are offline or otherwise unavailable. These nodes might be down for maintenance or temporarily out of service.

- **T (Total nodes)**: The total number of nodes in the partition. This is the sum of active, idle, and offline nodes.

### Example Interpretation

If the `NODES(A/I/O/T)` column shows `25/0/3/28`, it means:

- **25** nodes are currently active and running jobs.
- **0** nodes are idle and available for new jobs.
- **3** nodes are offline or unavailable.
- **28** nodes in total are part of this partition.

# VCF Contig and Chromosome Name Corrector Tutorial

When working with VCF files from different sources, chromosome naming inconsistencies can cause issues with tools like bcftools. This Python script helps standardize VCF files by correcting chromosome names and adding proper contig headers based on a reference genome.

The script handles several common issues:
- Missing contig headers in VCF files
- Inconsistent chromosome naming (e.g., "1" vs "chr1")
- Non-standard scaffold names
- Missing chromosome length information

Key features:
- Reads chromosome names and lengths from reference FASTA
- Creates mapping for different chromosome naming conventions
- Adds or updates contig information in VCF headers
- Filters out variants on chromosomes not present in reference
- Maintains VCF format compatibility

Usage example:
```bash
python VCF_contig_CHROM_correcter.py input.vcf output.vcf reference.fasta
```

Common workflow:
```bash
# 1. Fix chromosome names and add contig headers
python VCF_contig_CHROM_correcter.py input.vcf corrected.vcf GRCh38.fa

# 2. Sort the corrected VCF
bcftools sort corrected.vcf -o sorted.vcf

# 3. Normalize variants (optional)
bcftools norm -f reference.fa sorted.vcf -o normalized.vcf

# 4. Check the results
bcftools stats normalized.vcf > stats.txt
```

Tips for usage:
- Always backup your original VCF files before processing
- Verify chromosome names in your reference FASTA first
- Check the output VCF header to ensure contig lines are added correctly
- Use bcftools stats to validate the output file

Common chromosome naming variations handled:
- Numerical (1, 2, 3...)
- Chr-prefixed (chr1, chr2, chr3...)
- Sex chromosomes (X/chrX, Y/chrY)
- Mitochondrial (M/MT/chrM)
- Alternative contigs and scaffolds

The script creates a comprehensive mapping between different chromosome naming conventions and ensures that only variants on chromosomes present in your reference genome are retained. This standardization is particularly important when:
- Merging VCFs from different sources
- Running variant callers or analysis tools
- Comparing variants across datasets
- Preparing files for submission to databases

Error handling:
- Invalid VCF format lines are skipped
- Chromosomes not in reference are excluded
- Original headers are preserved (except contig lines)

Remember to check your reference genome's chromosome naming convention and adjust the script's mapping if needed for special cases in your data.

# ONT Direct RNA004 Pipeline Generator Documentation

This Python script [ONT_Direct_RNA004_steps_generator.py](/scripts/ONT_Direct_RNA004_steps_generator.py) generates a SLURM job for analyzing Oxford Nanopore Direct RNA004 sequencing data. The pipeline includes basecalled BAM file conversion, quality filtering, genome alignment, and quality control analysis.

Basic usage:
```bash
python ONT_Direct_RNA004_steps_generator.py \
    --bam /path/to/calls.bam \
    --reference /path/to/genome.mmi \
    --sample_name sample1
```

All available arguments:
```bash
--bam           : Input BAM file from Dorado basecalling (required)
--reference     : Reference genome minimap2 index path (required)
--min_quality   : Minimum read quality score (default: 8)
--min_length    : Minimum read length in bp (default: 500)
--threads       : Number of threads (default: 8)
--mem           : Memory allocation (default: 70G)
--time          : Wall time limit (default: 48:00:00)
--sample_name   : Sample name for output files (default: sample)
```

Example with all parameters:
```bash
python ONT_Direct_RNA004_steps_generator.py \
    --bam /projects/data/calls.bam \
    --reference /database/GRCh38.mmi \
    --min_quality 10 \
    --min_length 500 \
    --threads 8 \
    --mem 70G \
    --time 48:00:00 \
    --sample_name PC310cells
```

The script will generate a SLURM job script named `run_ont_rna_[sample_name].slurm` that performs:
- BAM to FASTQ conversion using samtools
- Read filtering with chopper
- Genome alignment using minimap2
- BAM file sorting and indexing
- Quality control with NanoPlot and Qualimap

To submit the generated SLURM job:
```bash
sbatch run_ont_rna_PC310cells.slurm
```

Output files will include:
- clean_reads.fastq: Filtered reads
- [sample_name]_direct_RNA.bam: Aligned reads
- [sample_name]_direct_RNA.bam.bai: BAM index
- nanoplot_qc_results/: NanoPlot quality reports
- qualimap/: Qualimap quality reports

# Quick guide for VCF file sort, zip, and indexing

## Installation
No additional installation is needed if you have standard bioinformatics tools:
- awk (built-in with Linux/Unix)
- bgzip
- tabix

## Usage
./[VCF_sort_bgzip_tabix.py](/scripts/VCF_sort_bgzip_tabix.py) -i input.vcf

## Input
- A standard VCF file (uncompressed)

## Output
The script will generate three files using the input file's prefix:
1. `{prefix}_sorted.vcf` - Sorted VCF file
2. `{prefix}_sorted.vcf.gz` - Compressed VCF file
3. `{prefix}_sorted.vcf.gz.tbi` - Tabix index file

## What it does
The script performs three sequential operations:
1. Sorts the VCF file while preserving headers
2. Compresses the sorted file using bgzip
3. Creates an index using tabix

## Example
```bash
./VCF_sort_bgzip_tabix.py -i na12878_truth.vcf
```
This will create:
- na12878_truth_sorted.vcf
- na12878_truth_sorted.vcf.gz
- na12878_truth_sorted.vcf.gz.tbi

# VCF Savior tool guide

[VCF_savior.py](/scripts/VCF_savior.py) is a comprehensive tool designed to fix common issues in VCF files and make them compatible with truvari for benchmarking structural variants. The tool automatically detects and fixes various problems including header definitions, genotype fields, and chromosome naming conventions.

Basic usage:
```bash
python VCF_savior.py -i input.vcf -o output.vcf
```

For standardizing chromosome names according to reference genome version:
```bash
python VCF_savior.py -i input.vcf -o output.vcf -g 38  # for GRCh38 (chr1, chrX format)
python VCF_savior.py -i input.vcf -o output.vcf -g 37  # for GRCh37 (1, X format)
```

The tool performs the following automatic fixes:
- Adds missing INFO and FORMAT definitions in the header
- Standardizes chromosome names for regular chromosomes (1-22, X, Y, M)
- Fixes genotype (GT) fields and ensures GT is the first FORMAT field
- Sets all FILTER fields to PASS
- Fixes SVLEN values based on SVTYPE and END positions
- Sorts and indexes the output VCF file

Output files:
- Fixed VCF file (.vcf)
- Sorted VCF file (_sorted.vcf)
- Compressed VCF file (_sorted.vcf.gz)
- Index file (_sorted.vcf.gz.tbi)

The tool is designed to be robust and tolerant of various VCF format issues, making it particularly useful when preparing VCF files for truvari benchmarking.

# Promethion SSH Connection Guide

This guide provides instructions for establishing an SSH connection to your Oxford Nanopore Promethion sequencing device.

## Prerequisites
- The Promethion device is powered on and connected to your network
- You know either the IP address or Device ID of your Promethion
- You have SSH client software installed on your computer

## Connection Steps

1. **Enable Remote Connection on the Promethion**
   - On the Promethion touchscreen, navigate to **Host Settings > System**
   - Find the **Remote Connection** toggle switch and enable it

2. **Connect via SSH**
   - Open a terminal and use one of the following commands:
     ```bash
     ssh prom@192.168.1.xxx  # Using IP Address
     ```
     or
     ```bash
     ssh prom@P2-######      # Using Device ID (preferred)
     ```
   - Replace `P2-######` with your actual Device ID

3. **Enter Password**
   - When prompted for a password, enter:
     ```
     prom
     ```