#!/bin/bash
# set up an ec2 instance with samplot and nvme mounted storage

set exuo pipefail

# some common dependencies ------------------------------------------------------------------------
sudo apt-get update
sudo apt-get install -y build-essential git libfuse-dev libcurl4-openssl-dev libxml2-dev mime-support automake libtool pkg-config libssl-dev ncurses-dev awscli python-pip libbz2-dev liblzma-dev unzip openjdk-8-jre-headless

# mount storage -----------------------------------------------------------------------------------
export TMPDIR=/mnt/local
sudo mkfs -t ext4 /dev/nvme0n1 # the device name changes from time to time (???)
# sudo mkfs -t ext4 /dev/nvme1n1
sudo mkdir /mnt/local
# sudo mount /dev/nvme1n1 /mnt/local
sudo mount /dev/nvme0n1 /mnt/local
sudo chown ubuntu /mnt/local

# conda setup -------------------------------------------------------------------------------------
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh
bash ~/miniconda.sh -b -p $HOME/miniconda
eval "$($HOME/miniconda/bin/conda shell.bash hook)"
conda init
conda config --add channels bioconda

# tmux/neovim setup -------------------------------------------------------------------------------
echo "source-file ~/.tmux.d/.tmux.conf" > ~/.tmux.conf
git clone https://github.com/mchowdh200/.tmux.d.git ~/.tmux.d

sudo add-apt-repository ppa:neovim-ppa/stable -y
sudo apt-get update -y
sudo apt-get install -y neovim
git clone https://github.com/mchowdh200/.vim.git ~/.vim
mkdir ~/.config
mkdir ~/.config/nvim
printf "set runtimepath^=~/.vim runtimepath+=~/.vim/after\nlet &packpath=&runtimepath\nsource ~/.vim/vimrc" > ~/.config/nvim/init.vim
pip install jedi neovim
echo "alias vim=nvim" >> ~/.profile
echo "export EDITOR=nvim" >> ~/.profile

# setup path --------------------------------------------------------------------------------------
mkdir /mnt/local/bin
echo "PATH=$PATH:/mnt/local/bin" >> ~/.profile

cd /mnt/local
# install htslib ----------------------------------------------------------------------------------
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

# install samtools --------------------------------------------------------------------------------
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

# install bcftools --------------------------------------------------------------------------------
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

# install duphold ---------------------------------------------------------------------------------
wget https://github.com/brentp/duphold/releases/download/v0.2.1/duphold -O /mnt/local/bin/duphold
chmod +x /mnt/local/bin/duphold

# install gargs -----------------------------------------------------------------------------------
wget https://github.com/brentp/gargs/releases/download/v0.3.9/gargs_linux -O /mnt/local/bin/gargs
chmod +x /mnt/local/bin/gargs

# setup conda environments
# samplot -----------------------------------------------------------------------------------------
conda create -y --name samplot
conda activate samplot
conda install -y python=3.6 matplotlib numpy cython
pip install pysam
git clone https://github.com/ryanlayer/samplot.git ~/samplot
conda deactivate

# smoove ------------------------------------------------------------------------------------------
# conda create -y --name smoove
# conda activate smoove
# conda install -y smoove

