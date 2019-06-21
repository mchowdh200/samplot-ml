#!/bin/bash
# don't forget to 'source activate samplot'

# c=$1 # CHROM
# s=$2 # START
# e=$3 # END
# n=$4 # SAMPLE
# g=$5 # GENOTYPE
CHROM=$1
START=$2
END=$3
SAMPLE=$4
GENOTYPE=$5

# bam file names (CSV)
# b=$( grep $n sample_bam_paths.txt )
BAMS=$( grep $n sample_bam_paths.txt )

# output file
OUT=imgs/${CHROM}_${START}_${END}_${SAMPLE}_${GENOTYPE}.png
echo $OUT

# TODO do I need to add reference (-r) for the crams
~/Repositories/samplot/src/samplot.py -c $CHROM -s $START -e $END -t DEL -b $BAMS -o $OUT
