#!/bin/bash

data_dir=/mnt/local/data
[[ ! -d $data_dir ]] && mkdir $data_dir

# get reference
aws s3 cp --recursive \
    s3://layerlabcu/ref/genomes/g1k_v37_decoy/ \
    $data_dir/g1k_v37_decoy/

# get BAM
[[ ! -d $data_dir/BAM ]] && mkdir $data_dir/BAM
aws s3 cp s3://layerlabcu/samplot-ml/HG002/BAM/hg002.bam $data_dir/BAM/
aws s3 cp s3://layerlabcu/samplot-ml/HG002/BAM/hg002.bam.bai $data_dir/BAM/

# get VCFs
[[ ! -d $data_dir/VCF ]] && mkdir $data_dir/VCF
aws s3 cp s3://layerlabcu/samplot-ml/manta_varied_breakpoints/variants-1.vcf $data_dir/VCF
aws s3 cp s3://layerlabcu/samplot-ml/manta_varied_breakpoints/variants-2.vcf $data_dir/VCF
aws s3 cp s3://layerlabcu/samplot-ml/manta_varied_breakpoints/variants-3.vcf $data_dir/VCF
aws s3 cp s3://layerlabcu/samplot-ml/manta_varied_breakpoints/variants-4.vcf $data_dir/VCF
aws s3 cp s3://layerlabcu/samplot-ml/manta_varied_breakpoints/variants-5.vcf $data_dir/VCF
aws s3 cp s3://layerlabcu/samplot-ml/manta_varied_breakpoints/variants0.vcf $data_dir/VCF
aws s3 cp s3://layerlabcu/samplot-ml/manta_varied_breakpoints/variants1.vcf $data_dir/VCF
aws s3 cp s3://layerlabcu/samplot-ml/manta_varied_breakpoints/variants2.vcf $data_dir/VCF
aws s3 cp s3://layerlabcu/samplot-ml/manta_varied_breakpoints/variants3.vcf $data_dir/VCF
aws s3 cp s3://layerlabcu/samplot-ml/manta_varied_breakpoints/variants4.vcf $data_dir/VCF
aws s3 cp s3://layerlabcu/samplot-ml/manta_varied_breakpoints/variants5.vcf $data_dir/VCF

# run duphold
[[ ! -d $data_dir/out ]] && mkdir $data_dir/out
for vcf in $(ls $data_dir/VCF); do
    duphold -v $data_dir/VCF/$vcf \
            -b $data_dir/BAM/hg002.bam \
            -f $data_dir/g1k_v37_decoy/g1k_v37_decoy.fa |
    bgzip -c > $data_dir/out/$(basename $vcf .vcf).duphold.vcf.gz
    tabix -f $data_dir/out/$(basename $vcf .vcf).duphold.vcf.gz
done
