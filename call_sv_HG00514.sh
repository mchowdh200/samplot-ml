#!/bin/bash

set -u

FASTA=~/data/ref/GRCh38_full_analysis_set_plus_decoy_hla/GRCh38_full_analysis_set_plus_decoy_hla.fa
CRAM=~/data/HG00514/CRAM/HG00514.alt_bwamem_GRCh38DH.20150715.CHS.high_coverage.cram
OUT_DIR=~/data/HG00514/VCF/smoove
SAMPLE=HG00514
EXCLUDE_BED=~/data/BED/exclude.cnvnator_100bp.GRCh38.20170403.bed

smoove call -x \
    --outdir $OUT_DIR \
    --name $SAMPLE \
    --fasta $FASTA \
    --exclude $EXCLUDE_BED \
    -p 1 \
    --genotype -d \
    $CRAM
