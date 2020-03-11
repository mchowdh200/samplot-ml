#!/bin/bash

set -eu

exclude_regions=$1
noexclude_regions=$2
true_pos_1kg=$3
out_dir=$4

# genomes related to test set
removed_genomes=(NA12878 NA12891 NA12892 
                 HG00512 HG00513 HG00514 
                 HG00731 HG00732 HG00733
                 NA19238 NA19239 NA19240)
removed_genomes=$(echo ${removed_genomes[@]} | tr " " "|") # construct regex from array

# sample false positives --------------------------------------------------------------------------
grep -Ev "$removed_genomes" $exclude_regions |
    grep -E '0\/1|1\/1' |
    awk '{ if ($6 > 0.7) print $0 }' | 
    python3 sample_ref_regions.py 105 1.0 "ref-fp" > $out_dir/ref-fp.bed
printf "ref-fp: "
wc -l $out_dir/ref-fp.bed

# sample true negatives ---------------------------------------------------------------------------
# regions where duphold disagrees with svtyper
grep -Ev "$removed_genomes" $noexclude_regions |
    awk '{ if ($5 <= 0.1 && $6 <= 0.7) print $0 }' | # AB=$5, DHFFC=$6
    python3 sample_ref_regions.py 1 0.25 "ref-tn" > $out_dir/ref-tn.bed
printf "ref-tn: "
wc -l $out_dir/ref-tn.bed

# regions where duphold agrees with svtyper
grep -Ev "$removed_genomes" $noexclude_regions |
    awk '{ if ($5 <= 0.1 && $6 > 0.7) print $0 }' |
    python3 sample_ref_regions.py 1 0.011 "ref-tn" >> $out_dir/ref-tn.bed
printf "ref-tn: "
wc -l $out_dir/ref-tn.bed


# sample true positives ---------------------------------------------------------------------------
# heterozygous deletions
grep -Ev "$removed_genomes" $true_pos_1kg |
    grep -E '0\|1|1\|0' |
    python3 sample_del.py 6 > $out_dir/het.bed
printf "tp-het: "
wc -l $out_dir/het.bed

# heterozygous deletions
grep -Ev "$removed_genomes" $true_pos_1kg |
    grep -E '1\|1' |
    python3 sample_del.py 7 > $out_dir/alt.bed
printf "tp-alt: "
wc -l $out_dir/alt.bed

