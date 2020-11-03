#!/bin/bash

# install cmake
sudo apt install cmake -y

# setup boost
cd /mnt/local
wget http://downloads.sourceforge.net/project/boost/boost/1.65.0/boost_1_65_0.tar.bz2
tar xf boost_1_65_0.tar.bz2
cd boost_1_65_0
./bootstrap.sh
./b2 --prefix=/mnt/local/boost_1_65_0_install link=static install

export BOOST_ROOT=/mnt/local/boost_1_65_0_install

cd /mnt/local
git clone https://github.com/Illumina/paragraph.git
cd paragraph
mkdir build
cd build
cmake /mnt/local/paragraph
make

# create conda env
conda create --name paragraph -y
source activate paragraph
conda install -y pysam intervaltree


