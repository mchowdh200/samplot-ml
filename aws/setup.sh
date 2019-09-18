#!/bin/bash
# set up an i3.large/xlarge instance with samplot and nvme mounted ssd

set exuo pipefail

# mount storage
export TMPDIR=/mnt/local
sudo apt-get update
sudo apt-get install -y build-essential git libfuse-dev libcurl4-openssl-dev libxml2-dev mime-support automake libtool pkg-config libssl-dev ncurses-dev awscli python-pip libbz2-dev liblzma-dev unzip openjdk-8-jre-headless
sudo mkfs -t ext4 /dev/nvme0n1
sudo mkdir /mnt/local
sudo mkdir /mnt/local
sudo mount /dev/nvme0n1 /mnt/local
sudo chown ubuntu /mnt/local

# setup path
mkdir /mnt/local/bin
export PATH="$PATH:/mnt/local/bin"

# setup samplot with anaconda
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh
bash ~/miniconda.sh -b -p $HOME/miniconda
eval "$($HOME/miniconda/bin/conda shell.bash hook)"
conda init
conda config --add channels bioconda
conda create -y --name samplot
conda activate samplot
conda install -y samplot

# other misc utilities/config
wget https://github.com/brentp/gargs/releases/download/v0.3.9/gargs_linux -O /mnt/local/bin/gargs
chmod +x /mnt/local/bin/gargs
