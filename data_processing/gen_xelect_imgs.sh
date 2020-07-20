#!/bin/bash

data_dir="../data_listings"
samples=( XA01 XA61 XA79 XE09 XE32 XH82 XH96 )

# get the regions
aws s3 cp s3://layerlabcu/samplot-ml/salmon/SVplaudit_Xelect_full.del.bed $data_dir

# get bam indices
for sample in ${samples[@]};do
    aws s3 cp s3://layerlabcu/xelect_bams/$sample.bam.bai .
done

# generate images from regions
sed -E -e 's/0\/1|1\/0/het/g' -e 's/1\/1/alt/g' -e 's/0\/0/ref/g' \
    $data_dir/SVplaudit_Xelect_full.del.bed |
    gargs --dry-run -p 4 "bash gen_img.sh \\
        --chrom {0} --start {1} --end {2} --sample {3} --genotype {4}-{5} \\
        --min-mqual 10 --bam-file s3://layerlabcu/xelect_bams/{3}.bam \\
        -o /mnt/local/data"


# for sample in ${samples[@]};do
#     # get the regions
#     aws s3 cp s3://layerlabcu/samplot-ml/salmon/$sample.bed $data_dir

#     # get the bam index
#     aws s3 cp s3://layerlabcu/xelect_bams/$sample.bam.bai .

#     # feed to gen_img
#     cat $data_dir/$sample.bed | gargs --dry-run -p 4 \
#         "bash gen_img.sh -c {0} -s {1} -e {2} -n $sample -g {3} \\
#                          -m 10 -b $sample.bam -o ~/data/"
# done
