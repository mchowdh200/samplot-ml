#!/bin/bash

data_dir=/mnt/local/data

# get reference
aws s3 cp --recursive s3://layerlabcu/ref/genomes/GRCh38_full_analysis_set_plus_decoy_hla/ \
    $data_dir/GRCh38_full_analysis_set_plus_decoy_hla/

# TODO get bam
aws s3 cp s3://layerlabcu/samplot-ml/HG002/BAM/hg002.bam $data_dir/BAM/
aws s3 cp s3://layerlabcu/samplot-ml/HG002/BAM/hg002.bam.bai $data_dir/BAM/

# download vcfs
aws s3 cp s3://layerlabcu/samplot-ml/manta_varied_breakpoints/variants-1.vcf $data_dir/VCF/
aws s3 cp s3://layerlabcu/samplot-ml/manta_varied_breakpoints/variants-2.vcf $data_dir/VCF/
aws s3 cp s3://layerlabcu/samplot-ml/manta_varied_breakpoints/variants-3.vcf $data_dir/VCF/
aws s3 cp s3://layerlabcu/samplot-ml/manta_varied_breakpoints/variants-4.vcf $data_dir/VCF/
aws s3 cp s3://layerlabcu/samplot-ml/manta_varied_breakpoints/variants-5.vcf $data_dir/VCF/
aws s3 cp s3://layerlabcu/samplot-ml/manta_varied_breakpoints/variants0.vcf $data_dir/VCF/
aws s3 cp s3://layerlabcu/samplot-ml/manta_varied_breakpoints/variants1.vcf $data_dir/VCF/
aws s3 cp s3://layerlabcu/samplot-ml/manta_varied_breakpoints/variants2.vcf $data_dir/VCF/
aws s3 cp s3://layerlabcu/samplot-ml/manta_varied_breakpoints/variants3.vcf $data_dir/VCF/
aws s3 cp s3://layerlabcu/samplot-ml/manta_varied_breakpoints/variants4.vcf $data_dir/VCF/
aws s3 cp s3://layerlabcu/samplot-ml/manta_varied_breakpoints/variants5.vcf $data_dir/VCF/

[[ ! -d $data_dir/out ]] && mkdir $data_dir/out
for vcf in $(ls $data_dir/VCF); do
    if [[ ! -d $data_dir/out/$(basename $vcf .vcf) ]]; then
        mkdir $data_dir/out/$(basename $vcf .vcf)
        mkdir $data_dir/out/$(basename $vcf .vcf)/imgs
        mkdir $data_dir/out/$(basename $vcf .vcf)/crop
    fi

    # generate images
    bcftools query -f '%CHROM\t%POS\t%INFO/END\n' $data_dir/$vcf |
        gargs -p 36 -e \
            "bash gen_img.sh \\
                --chrom {0} --start {1} --end {2} \\
                --sample HG002 --genotype DEL \\
                --min-mqual 10 \\
                --bam-file $data_dir/BAM/hg002.bam \\
                --out-dir $data_dir/out/$(basename $vcf .vcf)/imgs"
    # crop images
    bash crop.sh -p 36 -d $data_dir/out/$(basename $vcf .vcf)/imgs -o $data_dir/out/$(basename $vcf .vcf)/crop
done
