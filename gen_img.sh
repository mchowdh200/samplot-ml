#!/bin/bash
set -eu

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
BAM_LIST=~/Repositories/samplot-ml/data/cram-indices/cram-index-paths.txt # yeah I know they're crams...
BAMS=$(grep $SAMPLE $BAM_LIST)

# output file
OUT=$OUT_DIR/${CHROM}_${START}_${END}_${SAMPLE}_${GENOTYPE}.png
echo $OUT

if [ ! -f $OUT ]; then
    ~/Repositories/samplot/src/samplot.py -c $CHROM -s $START -e $END -t DEL -b $BAMS -o $OUT -r $FASTA
fi
cd $START_DIR


#cat data/del.sample.bed | head | gargs 'rj -l log/ -n {0}_{1}_{2}_{3}_{4} -c "bash gen_img.sh {0} {1} {2} {3} {4}"'

