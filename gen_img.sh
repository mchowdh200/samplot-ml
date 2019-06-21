#!/bin/bash

# TODO don't forget to 'source activate samplot'

# ideally you will pass in these args by using xargs/gargs/parallel
# along with a bed file with corresponding tab delimitted columns.
CHROM=chr$1 # our reference genome has chr as a prefix to the contig name
START=$2
END=$3
SAMPLE=$4
GENOTYPE=$5

START_DIR=$PWD
CRAI=data/cram-indices
OUT_DIR=/scratch/Shares/layer/projects/samplot/ml/data/1kg/high_cov/imgs
cd $CRAI

# reference path
FASTA=/scratch/Shares/layer/ref/GRCh38_full_analysis_set_plus_decoy_hla/GRCh38_full_analysis_set_plus_decoy_hla.fa

# bam file names
# b=$( grep $n sample_bam_paths.txt )
BAM_LIST=data/cram-indices/cram-index-paths.txt # yeah I know they're crams...
BAMS=$(grep $n $BAM_LIST)

# output file
OUT=$OUT_DIR/${CHROM}_${START}_${END}_${SAMPLE}_${GENOTYPE}.png
echo $OUT

# TODO do I need to add reference (-r) for the crams

source activate samplot
~/Repositories/samplot/src/samplot.py -c $CHROM -s $START -e $END -t DEL -b $BAMS -o $OUT
source deactivate

cd $START_DIR


# cat data/del.sample.bed | gargs -p 1 "get_img {0} {1} {2} {3} {4}"


