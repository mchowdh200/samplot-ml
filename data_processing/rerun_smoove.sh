#!/bin/bash

function call_svs(){
    set -eu
    local bam=$1
    local sample=$(basename $bam .cram | cut -d'.' -f1)
    local fasta=$2
    local exclude=$3
    local outdir=$4

    if [[ ! -d $outdir ]]; then
        mkdir $outdir
        mkdir $outdir/$sample
    fi

    smoove call \
        --name $sample \
        --fasta $fasta \
        --exclude $exclude \
        --genotype \
        --duphold \
        --removepr \
        --outdir $outdir/$sample \
        $bam
}
export -f call_svs

bam_dir="/mnt/local/BAM"
fasta="/mnt/local/ref/GRCh38_full_analysis_set_plus_decoy_hla/GRCh38_full_analysis_set_plus_decoy_hla.fa"
exclude="/mnt/local/BED/exclude.cnvnator_100bp.GRCh38.20170403.bed"
outdir="/mnt/local/VCF"

ls $bam_dir/*.cram | gargs -p 3 "call_svs {} $fasta $exclude $outdir"
