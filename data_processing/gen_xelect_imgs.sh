#!/bin/bash

data_dir="../data_listings"
out_dir="/mnt/local/data"
samples=( XA01 XA61 XA79 XE09 XE32 XH82 XH96 )

for sample in ${samples[@]};do
    # get the regions
    aws s3 cp s3://layerlabcu/samplot-ml/salmon/$sample.bed $data_dir

    # get the bam index
    aws s3 cp s3://layerlabcu/xelect_bams/$sample.bam.bai .

    # feed to gen_img
    cat $data_dir/$sample.bed | gargs -p 4 \
        "bash gen_img.sh -c {0} -s {1} -e {2} -n $sample -g {3} \\
                         -m 10 -b s3://layerlabcu/xelect_bams/$sample.bam -o $out_dir/imgs"

    bash crop.sh -p 4 -d $out_dir/imgs -o $out_dir/crop
done

aws s3 cp --recursive $out_dir/imgs/ s3://layerlabcu/samplot-ml/salmon/imgs/
aws s3 cp --recursive $out_dir/crop/ s3://layerlabcu/samplot-ml/salmon/crop/
