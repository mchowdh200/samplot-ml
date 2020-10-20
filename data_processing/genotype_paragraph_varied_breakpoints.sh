#!/bin/bash

set -eu

[[ ! -d /mnt/local/data ]] && mkdir /mnt/local/data

out_dir=/mnt/local/data/out
mkdir $out_dir

# get reference
aws s3 cp --recursive \
    s3://layerlabcu/ref/genomes/g1k_v37_decoy/ \
    /mnt/local/data/g1k_v37_decoy/

# get BAM
aws s3 cp --recursive \
    s3://layerlabcu/samplot-ml/HG002/BAM/ \
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

# TODO For each vcf: run paragraph

for vcf in $(ls /mnt/local/data/VCF); do
    [[ ! -d /mnt/local/data/out/$(basename $vcf) ]] && mkdir /mnt/loca/data/out/$(basename $vcf)
    python3 ~/paragraph/build/bin/multigrmpy.py \
        -i $vcf \
        -m /mnt/local/data/samples.txt \
        -r /mnt/local/data/g1k_v37_decoy/g1k_v37_decoy.fa \
        -o /mnt/local/data/out/$(basename $vcf)
done
