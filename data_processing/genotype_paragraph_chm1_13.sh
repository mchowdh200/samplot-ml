#!/bin/bash

set -eu

data_dir=/mnt/local/data
[[ ! -d $data_dir ]] && mkdir $data_dir

out_dir=$data_dir/out
[[ ! -d $out_dir ]] && mkdir $out_dir

# =============================================================================
# get data
# =============================================================================
# reference
if [[ ! -d $data_dir/hg38 ]]; then
    aws s3 cp --recursive s3://layerlabcu/ref/genomes/hg38/ $data_dir/hg38/
fi
ref=$data_dir/hg38/hg38.fa

# BAMs
# TODO test with just one BAM
# aws s3 cp s3://layerlabcu/samplot-ml/CHM1_13/BAM/CHM13.5-CHM1.5.bam $data_dir/BAM/
# aws s3 cp s3://layerlabcu/samplot-ml/CHM1_13/BAM/CHM13.5-CHM1.5.bam.bai $data_dir/BAM/

# if [[ ! -d $data_dir/BAM ]]; then
#     aws s3 cp --recursive s3://layerlabcu/samplot-ml/CHM1_13/BAM/ $data_dir/BAM/
# fi

# VCFs
# if [[ ! -d $data_dir/VCF ]]; then
    # mkdir $data_dir/VCF
    aws s3 cp s3://layerlabcu/samplot-ml/CHM1_13/VCF/CHM1.del.unique.replaced_sample.vcf.gz $data_dir/VCF/
    aws s3 cp s3://layerlabcu/samplot-ml/CHM1_13/VCF/CHM13.del.unique.replaced_sample.vcf.gz $data_dir/VCF/

    tabix -f $data_dir/VCF/CHM1.del.unique.replaced_sample.vcf.gz
    tabix -f $data_dir/VCF/CHM13.del.unique.replaced_sample.vcf.gz
# fi

# =============================================================================
# get read depth/read length for bams and create manifests for paragraph runs.
# =============================================================================

# for each mixture BAM
# for bam in $(find $data_dir/BAM/ -name '*.bam'); do
#     read_depth=$(/mnt/local/paragraph/build/bin/idxdepth -b $bam -r $ref |
#         grep -m 1 'depth' | sed -E 's/\s+//g' | cut -d':' -f2)
#     read_length=$(samtools view $bam |
#         awk '{print length($10)}' | head -100000|
#         python3 -c 'import sys, statistics; print(statistics.mean([float(i.rstrip()) for i in sys.stdin]))')
#     printf "id\tpath\tdepth\tread length\nCHM1_13\t$bam\t$read_depth\t$read_length\n" > $data_dir/$(basename $bam .bam).txt

# done


# =============================================================================
# Paragraph runs
# =============================================================================

for vcf in $(find $data_dir/VCF/ -name '*.vcf.gz'); do
    sample=$(cut -d'.' -f1 <<<$(basename $vcf))
    echo $sample =============================================================
    [[ ! -d $out_dir/$sample ]] && mkdir $out_dir/$sample

    for manifest in $(find $data_dir/ -name '*.txt'); do
        python3 /mnt/local/paragraph/build/bin/multigrmpy.py \
            -i $vcf \
            -m $manifest \
            -r $ref \
            -o $out_dir/$sample/$(basename $manifest .txt)
    done
done
