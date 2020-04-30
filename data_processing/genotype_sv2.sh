#!/bin/bash

# 1. call variants with lumpy script (found in the smoove directory)
# 2. using the output vcf + snv vcf + ped, genotype using sv2

set -euo pipefail

function call_and_genotype {
    local data_dir=$1
    local sample=$2

    # paths for reference genomes
    hg19_fasta=$3
    hg38_fasta=$4
    hg19_exclude=$5
    hg38_exclude=$6

    [[ ! -d $data_dir/$sample/VCF/SV2 ]] && mkdir $data_dir/$sample/VCF/SV2

    if [[ $sample == "HG002" ]]; then
        genome="hg19"
        fasta=$hg19_fasta
        snv_vcf="$data_dir/$sample/VCF/HG002_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-22_v.3.3.2_highconf_triophased.vcf.gz"
        exclude=$hg19_exclude
    else 
        genome="hg38"
        fasta=$hg38_fasta
        exclude=$hg38_exclude
        if [[ $sample == "HG00514" ]]; then
            snv_vcf="$data_dir/$sample/VCF/20170309_han.gatk.vcf.gz"
        elif [[ $sample == "HG00733" ]]; then
            snv_vcf="$data_dir/$sample/VCF/20170309_pur.gatk.vcf.gz"
        elif [[ $sample == "NA19240" ]]; then
            snv_vcf="$data_dir/$sample/VCF/20170309_yri.gatk.vcf.gz"
        fi
    fi

    bam=$(find $data_dir/$sample/BAM/ -name '*.bam' -or -name '*.cram')
    
    # call svs with smoove
    # smoove call $bam \
    #     -x \
    #     -n $sample \
    #     -f $fasta \
    #     -e $exclude \
    #     -o $data_dir/$sample/VCF/SV2 \

    # genotype with sv2
    sv2 -$genome $fasta
    sv2 -g $genome \
        -i $bam \
        -v $data_dir/$sample/VCF/SV2/$sample-smoove.vcf.gz \
        -snv $snv_vcf \
        -p $data_dir/$sample/VCF/$sample.ped \
        -odir $data_dir/$sample/VCF/SV2 \
        -o ${sample}_sv2
}
export -f call_and_genotype



data_dir=$1
# samples="HG002\nHG00514\nHG00733\nNA19240"
samples="HG00514\nHG00733\nNA19240"

# call svs with lumpy (and existing split/discordant reads from previous smoove run)
printf "$samples" | gargs -p 4 -e \
    "call_and_genotype $data_dir {} \\
        ~/data/ref/g1k_v37_decoy/g1k_v37_decoy.fa \\
        ~/data/ref/GRCh38_full_analysis_set_plus_decoy_hla/GRCh38_full_analysis_set_plus_decoy_hla.fa \\
        ~/data/BED/ceph18.b37.lumpy.exclude.2014-01-15.bed \\
        ~/data/BED/exclude.cnvnator_100bp.GRCh38.20170403.bed"

# genotype with sv2
