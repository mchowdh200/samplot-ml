#!/bin/bash

set -eu

exclude_regions_bed=$1
true_pos_1kg=$2
out_dir=$3

# genomes related to test set
removed_genomes=(NA12878 NA12891 NA12892 
                 HG00512 HG00513 HG00514 
                 HG00731 HG00732 HG00733
                 NA19238 NA19239 NA19240)
removed_genomes=$(echo ${removed_genomes[@]} | tr " " "|") # construct regex from array

# # sample false positives --------------------------------------------------------------------------
grep -Ev "$removed_genomes" $exclude_regions_bed |
    grep -E '0\/1|1\/1' |
    python3 sample_exclude_regions.py 65 1.0 > $out_dir/ref.bed

# # sample true negatives ---------------------------------------------------------------------------
grep -Ev "$removed_genomes" $exclude_regions_bed |
    grep -E '0\/0' |
    python3 sample_exclude_regions.py 1 0.05 >> $out_dir/ref.bed

# sample true positives ---------------------------------------------------------------------------
# heterozygous deletions
grep -Ev "$removed_genomes" $true_pos_1kg |
    grep -E '0\|1|1\|0' |
    python3 sample_del.py 6 > $out_dir/het.bed

# heterozygous deletions
grep -Ev "$removed_genomes" $true_pos_1kg |
    grep -E '1\|1' |
    python3 sample_del.py 7 > $out_dir/alt.bed
