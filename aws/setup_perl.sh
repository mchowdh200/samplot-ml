#!/bin/bash

# setup bio perl and other dependencies for cnv-jacg package

# 1) perl stuff
sudo apt install bioperl -y
perl -MCPAN -e 'install Bio::DB::HTS'
perl -MCPAN -e 'Statistics::Basic'

# 2) R stuff
sudo apt install r-base -y
sudo su - -c "R -e \"install.packages('randomForest', repos='http://cran.rstudio.com/')\""

# 3) cnv-jacg
wget https://github.com/sunnyzxh/CNV-JACG/archive/v1.1.zip
unzip v1.1.zip
