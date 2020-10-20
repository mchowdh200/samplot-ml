#!/bin/bash

set -eu

[[ ! -d /mnt/local/data ]] && mkdir /mnt/local/data

[[ ! -d /mnt/local/data/out ]] && mkdir /mnt/local/data/out

# get reference
aws s3 cp --recursive \
    s3://layerlabcu/ref/genomes/g1k_v37_decoy/ \
    /mnt/local/data/g1k_v37_decoy/

# get BAM
aws s3 cp s3://layerlabcu/samplot-ml/HG002/BAM/hg002.bam \
    /mnt/local/data/BAM/
aws s3 cp s3://layerlabcu/samplot-ml/HG002/BAM/hg002.bam.bai \
    /mnt/local/data/BAM/

# get vcfs
aws s3 cp --recursive \
    s3://layerlabcu/samplot-ml/manta_varied_breakpoints/ \
    /mnt/local/data/VCF/

# get read_depth/read_length
read_depth=$(~/paragraph/build/bin/idxdepth \
    -b /mnt/local/data/BAM/hg002.bam \
    -r /mnt/local/data/g1k_v37_decoy/g1k_v37_decoy.fa |
    grep -m 1 'depth' | sed -E 's/\s+//g' | cut -d':' -f2)
read_length=$(samtools view /mnt/local/data/BAM/hg002.bam |
    awk '{print length($10)}' | head -100000 |
    python3 -c 'import sys, statistics; print(statistics.mean([float(i.rstrip()) for i in sys.stdin]))')
printf "id\tpath\tdepth\tread length\nHG002\t/mnt/local/data/BAM/hg002.bam\t$read_depth\t$read_length\n" > /mnt/local/data/samples.txt

# ===
# TEST ON THE ORIG VCF
# [[ ! -d /mnt/local/data/out/variants0 ]] && mkdir /mnt/local/data/out/variants0
# python3 ~/paragraph/build/bin/multigrmpy.py \
#     -i /mnt/local/data/VCF/variants0.vcf \
#     -m /mnt/local/data/samples.txt \
#     -r /mnt/local/data/g1k_v37_decoy/g1k_v37_decoy.fa \
#     -o /mnt/local/data/out/variants0 &> log.txt
# ===


for vcf in $(ls /mnt/local/data/VCF); do
    echo $vcf ---------------------------------------------
    [[ ! -d /mnt/local/data/out/$(basename $vcf .vcf) ]] && mkdir /mnt/local/data/out/$(basename $vcf .vcf)
    python3 ~/paragraph/build/bin/multigrmpy.py \
        -i /mnt/local/data/VCF/$vcf \
        -m /mnt/local/data/samples.txt \
        -r /mnt/local/data/g1k_v37_decoy/g1k_v37_decoy.fa \
        -o /mnt/local/data/out/$(basename $vcf .vcf) &> log.txt
done
