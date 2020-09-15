#!/bin/bash

# set -e

data_dir=$1
samples=(HG00514 HG00733 NA19240)

# get vcfs
[[ ! -d $data_dir/VCF ]] && mkdir $data_dir/VCF
aws s3 cp --recursive s3://layerlabcu/samplot-ml/chaisson/VCF/smoove/ $data_dir/VCF
vcfs=$(find $data_dir/VCF -name '*.vcf.gz')

# get reference
[[ ! -d $data_dir/ref ]] && mkdir $data_dir/ref
aws s3 cp --recursive s3://layerlabcu/samplot-ml/ref/GRCh38_full_analysis_set_plus_decoy_hla/ $data_dir/ref/
fasta=$data_dir/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa

# get alignments
[[ ! -d $data_dir/BAM ]] && mkdir $data_dir/BAM
aws s3 cp --recursive s3://layerlabcu/samplot-ml/chaisson/BAM/ $data_dir/BAM/
bams=$(find $data_dir/BAM -name '*.cram')

# generate images
for sample in ${samples[@]}; do
    echo $sample -------------------------------
    if [[ ! -d $data_dir/$sample ]]; then
        mkdir $data_dir/$sample
        mkdir $data_dir/$sample/imgs
        mkdir $data_dir/$sample/crop
    fi

    vcf=$(echo "$vcfs" | grep $sample)
    bam=$(echo "$bams" | grep $sample)

    bcftools query -f '%CHROM\t%POS\t%INFO/END\n' $vcf |
        gargs -p 36 \
            "bash gen_img.sh \\
                --chrom {0} --start {1} --end {2} \\
                --sample $sample --genotype DEL \\
                --min-mqual 10 \\
                --fasta $fasta \\
                --bam-file $bam \\
                --out-dir $data_dir/$sample/imgs"

    bash crop.sh -p 36 -d $data_dir/$sample/imgs -o $data_dir/$sample/crop
done
