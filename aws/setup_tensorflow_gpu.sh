#!/bin/bash
# set up an ec2 instance with samplot and nvme mounted storage

set exuo pipefail

# some common dependencies ------------------------------------------------------------------------
sudo apt-get update
sudo apt-get install -y build-essential git libfuse-dev libcurl4-openssl-dev libxml2-dev mime-support automake libtool pkg-config libssl-dev ncurses-dev awscli  libbz2-dev liblzma-dev unzip

# mount storage -----------------------------------------------------------------------------------
sudo mkfs -t ext4 /dev/xvdf
sudo mkdir /mnt/local
sudo mount /dev/xvdf /mnt/local
sudo chown ubuntu /mnt/local

# nvidia stuff ------------------------------------------------------------------------------------
# Add NVIDIA package repositories
wget https://developer.download.nvidia.com/compute/cuda/repos/ubuntu1804/x86_64/cuda-repo-ubuntu1804_10.0.130-1_amd64.deb
sudo dpkg -i cuda-repo-ubuntu1804_10.0.130-1_amd64.deb
sudo apt-key adv --fetch-keys https://developer.download.nvidia.com/compute/cuda/repos/ubuntu1804/x86_64/7fa2af80.pub
sudo apt-get update
wget http://developer.download.nvidia.com/compute/machine-learning/repos/ubuntu1804/x86_64/nvidia-machine-learning-repo-ubuntu1804_1.0.0-1_amd64.deb
sudo apt install ./nvidia-machine-learning-repo-ubuntu1804_1.0.0-1_amd64.deb
sudo apt-get update

# Install NVIDIA driver
sudo apt-get install -y --no-install-recommends nvidia-driver-418
# Reboot. Check that GPUs are visible using the command: nvidia-smi

# Install development and runtime libraries (~4GB)
# sudo apt-get install -y --no-install-recommends \
#     cuda-10-0 \
#     libcudnn7=7.6.2.24-1+cuda10.0  \
#     libcudnn7-dev=7.6.2.24-1+cuda10.0

# Install TensorRT. Requires that libcudnn7 is installed above.
# sudo apt-get install -y --no-install-recommends libnvinfer5=5.1.5-1+cuda10.0 \
#     libnvinfer-dev=5.1.5-1+cuda10.0


# conda setup -------------------------------------------------------------------------------------
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh
bash ~/miniconda.sh -b -p /mnt/local/miniconda
eval "$(/mnt/local/miniconda/bin/conda shell.bash hook)"
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


# install gargs -----------------------------------------------------------------------------------
wget https://github.com/brentp/gargs/releases/download/v0.3.9/gargs_linux -O /mnt/local/bin/gargs
chmod +x /mnt/local/bin/gargs

# setup conda environments
# samplot -----------------------------------------------------------------------------------------
conda create -y --name tensorflow
conda activate tensorflow
conda install -y python=3.7 pip matplotlib numpy cudatoolkit=10.0 cudnn joblib scikit-learn
pip install tensorflow-gpu tensorflow-addons
conda deactivate

# setup environment variables that will enable use of s3 buckets with tensorflow datasets ---------
echo "export AWS_REGION=us-east-1" >> ~/.profile
echo "export S3_USE_HTTPS=1" >> ~/.profile
echo "export S3_VERIFY_SSL=1" >> ~/.profile


