#!/bin/bash

# install cmake
sudo apt install cmake -y

# setup boost
cd ~
wget http://downloads.sourceforge.net/project/boost/boost/1.65.0/boost_1_65_0.tar.bz2
tar xf boost_1_65_0.tar.bz2
cd boost_1_65_0
./bootstrap.sh
./b2 --prefix=$HOME/boost_1_65_0_install link=static install

export BOOST_ROOT=$HOME/boost_1_65_0_install


cd ~
git clone https://github.com/Illumina/paragraph.git
cd paragraph
mkdir build
cd build
cmake ~/paragraph
make

# TODO create conda env
# install pysam
# intervaltree
conda create --name paragraph -y
conda activate paragrph
conda install pysam intervaltree


