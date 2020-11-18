#!/bin/bash

data_dir=/mnt/local/data
[[ ! -d $data_dir ]] && mkdir $data_dir

# get ref
if [[ ! -d $data_dir/ref ]]; then
    aws s3 cp --recursive \
        s3://layerlabcu/ref/genomes/GRCh38_full_analysis_set_plus_decoy_hla/ \
        $data_dir/ref/
fi
ref=$data_dir/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa

# get BAM
if [[ ! -d $data_dir/BAM/ ]]; then
    aws s3 cp s3://layerlabcu/samplot-ml/chaisson/BAM/HG00514.alt_bwamem_GRCh38DH.20150715.CHS.high_coverage.cram $data_dir/BAM/
    aws s3 cp s3://layerlabcu/samplot-ml/chaisson/BAM/HG00514.alt_bwamem_GRCh38DH.20150715.CHS.high_coverage.cram.crai $data_dir/BAM/
fi
bam=$data_dir/BAM/HG00514.alt_bwamem_GRCh38DH.20150715.CHS.high_coverage.cram

# get BEDs
[[ ! -d $data_dir/BED/ ]] &&
    aws s3 cp --recursive s3://layerlabcu/samplot-ml/peter_SV/ $data_dir/BED/

# gen images for each bed
for bed in $(ls $data_dir/BED); do
    out_dir=$(basename $bed .bed)
    img_type=$(cut -d'.' -f1-2 <<<$out_dir)
    [[ ! -d $data_dir/$out_dir ]] && mkdir $out_dir

    # cat bed into gargs
    # pipe into gen_img.sh
    cat $bed | gargs -p 36 \
        "bash gen_img.sh \\
            -c {0} -s {1} -e {2} -n HG00514 -g $img_type \\
            -m 10 -f $ref -b $bam -o $out_dir"
done
