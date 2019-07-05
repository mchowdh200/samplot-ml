#!/bin/bash

set -e

MODEL_PATH=$1
DATA_LIST=data/giab/giab_list.txt

VCF_DIR=data/giab/VCF/$(basename $MODEL_PATH .h5)
VCF_FILE=data/giab/VCF/HG002-smoove.genotyped.vcf.gz
# BED_FILE=data/giab/BED/HG002_SVs_Tier1_v0.6.bed
BED_FILE=predictions/$(basename $MODEL_PATH .h5).bed
if [ ! -d $VCF_DIR ]; then
    mkdir $VCF_DIR
fi

source activate tensorflow2-beta
python3 run.py predict -mp $MODEL_PATH -h5 -i $DATA_LIST \
    | python3 pred2bed.py > predictions/$(basename $MODEL_PATH .h5).bed
source deactivate

source activate cyvcf2
python3 annotate.py $VCF_FILE $BED_FILE | bgzip -c > $VCF_DIR/ml.vcf.gz
bash filter_vcf.sh $VCF_DIR
source deactivate




