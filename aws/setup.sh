#!/bin/bash
# set up an i3.large/xlarge instance with samplot and nvme mounted ssd

set exuo pipefail

# some common dependencies
sudo apt-get update
sudo apt-get install -y build-essential git libfuse-dev libcurl4-openssl-dev libxml2-dev mime-support automake libtool pkg-config libssl-dev ncurses-dev awscli python-pip libbz2-dev liblzma-dev unzip openjdk-8-jre-headless

# mount storage
export TMPDIR=/mnt/local
sudo mkfs -t ext4 /dev/nvme0n1
sudo mkdir /mnt/local
sudo mkdir /mnt/local
sudo mount /dev/nvme0n1 /mnt/local
sudo chown ubuntu /mnt/local

# anaconda
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh
bash ~/miniconda.sh -b -p $HOME/miniconda
eval "$($HOME/miniconda/bin/conda shell.bash hook)"
conda init
conda config --add channels bioconda

# tmux/neovim setup
echo "source-file ~/.tmux.d/.tmux.conf" > ~/.tmux.conf
git clone https://github.com/mchowdh200/.tmux.d.git ~/.tmux.d

sudo add-apt-repository ppa:neovim-ppa/stable -y
sudo apt-get update -y
sudo apt-get install -y neovim
git clone https://github.com/mchowdh200/.vim.git ~/.vim
mkdir ~/.config
mkdir ~/.config/nvim
printf "set runtimepath^=~/.vim runtimepath+=~/.vim/after\nlet &packpath=&runtimepath\nsource ~/.vimrc" > ~/.config/nvim/init.vim
pip install jedi neovim
echo "alias vim=nvim" >> ~/.profile
echo "export EDITOR=nvim" >> ~/.profile

# setup path
mkdir /mnt/local/bin
echo "PATH=$PATH:/mnt/local/bin" >> ~/.profile

# setup conda environments
conda create -y --name samplot
conda activate samplot
conda install -y samplot

conda create -y --name smoove
conda activate smoove
conda install -y smoove

# other misc utilities/config
wget https://github.com/brentp/gargs/releases/download/v0.3.9/gargs_linux -O /mnt/local/bin/gargs
chmod +x /mnt/local/bin/gargs
