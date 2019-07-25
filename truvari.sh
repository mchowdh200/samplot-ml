#!/bin/bash
# VCF_DIR=data/giab/VCF
# BED_DIR=data/giab/BED
# lumpy_vcf=$VCF_DIR/HG002-smoove.genotyped.vcf.gz
# fasta=/scratch/Shares/layer/ref/hs37-1kg/human_g1k_v37.20.fasta
# fasta=data/FASTA/human_g1k_v37.20.fasta

out_dir=duphold_comp
rm -rf $out_dir

comp_vcf=$1
truth_set=$2
fasta=$3

truvari \
    -b $truth_set \
    -c $comp_vcf \
    -o $out_dir \
    --pctsim 0 \
    --sizemax 15000000 \
    --sizemin 300 --sizefilt 270 \
    --pctov 0.6 \
    --refdist 20 \
    --noprog \
    --no-ref a
    # --debug
    #--giabreport \

    #--reference $fasta \


    # -b $VCF_DIR/HG002_SVs_Tier1_v0.6.DEL.vcf.gz \
    # --includebed $BED_DIR/HG002_SVs_Tier1_v0.6.bed \
