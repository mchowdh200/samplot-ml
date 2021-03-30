#!/bin/bash
# set up an ec2 instance with samplot and nvme mounted storage

set exuo pipefail

## some common dependencies
sudo apt-get update
sudo apt-get install -y \
    build-essential curl git libcurl4-openssl-dev wget \
    software-properties-common automake libtool pkg-config libssl-dev \
    ncurses-dev awscli python-pip libbz2-dev liblzma-dev unzip

## mount storage
export TMPDIR=/mnt/local/temp
sudo mkdir /mnt/local

# setup raid 0 if more than one drive specified
# the nvme drive naming convention is not consistent enough
# so I have just resorted to filtering out nvme disks with
# the expected size (ex smallest c5d has a 50GB so > 40 is my threshold)
num_drives=$(lsblk -o NAME,SIZE | grep 'nvme'| awk '$2 ~ /G$/ && $2+0 > 40' | wc -l)
if [[ $num_drives > 1 ]]; then
    drive_list=$(lsblk -o NAME,SIZE | grep 'nvme' |
                 awk '$2 ~ /G$/ && $2+0 > 40' |
                 awk 'BEGIN{ORS=" "}{print "/dev/"$1 }')
    sudo mdadm --create --verbose \
         /dev/md0 \
         --level=0 \
         --raid-devices=$num_drives $drive_list
    sudo mkfs -t ext4 /dev/md0
    sudo mount /dev/md0 /mnt/local
else
    sudo mkfs -t ext4 /dev/nvme0n1 
    sudo mount /dev/nvme0n1 /mnt/local
fi

sudo chown ubuntu /mnt/local
mkdir /mnt/local/data
mkdir /mnt/local/temp
echo "export TMPDIR=/mnt/local/temp" >> ~/.profile

## tmux/neovim setup
echo 'source-file ~/.tmux.d/.tmux.conf' > ~/.tmux.conf
git clone https://github.com/mchowdh200/.tmux.d.git ~/.tmux.d

sudo add-apt-repository ppa:neovim-ppa/stable -y
sudo apt-get update -y
sudo apt-get install -y neovim
git clone https://github.com/mchowdh200/.vim.git ~/.vim
mkdir ~/.config
mkdir ~/.config/nvim
printf 'set runtimepath^=~/.vim runtimepath+=~/.vim/after\nlet &packpath=&runtimepath\nsource ~/.vim/vimrc' > ~/.config/nvim/init.vim
pip install jedi neovim
echo 'alias vim=nvim' >> ~/.profile
echo 'export EDITOR=nvim' >> ~/.profile

## setup path
mkdir /mnt/local/bin
echo 'PATH=$PATH:/mnt/local/bin' >> ~/.profile

## install gargs
wget https://github.com/brentp/gargs/releases/download/v0.3.9/gargs_linux -O /mnt/local/bin/gargs
chmod +x /mnt/local/bin/gargs

## conda setup
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /mnt/local/miniconda.sh
bash /mnt/local/miniconda.sh -b -p /mnt/local/miniconda
eval "$(/mnt/local/miniconda/bin/conda shell.bash hook)"
conda init
conda config --add channels bioconda
conda install -y -c conda-forge mamba

## setup snakemake
mamba create -c conda-forge -c bioconda -y -n snakemake snakemake boto3

## setup htslib
cd /mnt/local
wget https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2
bunzip2 htslib-1.9.tar.bz2
tar -xvf htslib-1.9.tar
cd htslib-1.9
autoheader
autoconf
./configure --enable-libcurl --enable-s3
make
sudo make install
echo "export LD_LIBRARY_PATH=/usr/local/lib" >> ~/.profile
echo "export HTSLIB_LIBRARY_DIR=/usr/local/lib" >> ~/.profile
echo "export HTSLIB_INCLUDE_DIR=/usr/local/include" >> ~/.profile
cd -

## setup samtools
wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2
bunzip2 samtools-1.9.tar.bz2
tar -xvf samtools-1.9.tar
cd samtools-1.9
autoheader
autoconf -Wno-syntax
./configure --with-htslib=system --enable-configure-htslib
make
sudo make install
cd -

## setup bcftools
wget https://github.com/samtools/bcftools/releases/download/1.9/bcftools-1.9.tar.bz2
bunzip2 bcftools-1.9.tar.bz2
tar -xvf bcftools-1.9.tar
cd bcftools-1.9
autoheader
autoconf -Wno-syntax
./configure --with-htslib=system --enable-configure-htslib
make
sudo make install
cd -
