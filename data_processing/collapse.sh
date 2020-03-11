#!/bin/bash

# take a sorted(!) bed file with the following format
# chr  start  end  sample  AB  DHFFC GT
# and intersect it with itself.
# prints two sets of intervals side by side
# if the chr start end in the first interval are the 
# same as the previous line's, then only change the end 
# of the second interval.
set -eu

function self_overlap() {
    bed=$1

    bedtools intersect \
        -a $bed \
        -b $bed \
        -f 0.7 -r -wa -wb |
        awk 'BEGIN{OFS="\t"} {
            if ($1==chrA && $2==startA && $3==endA) {
                endB=$10
            } else {
                print chrA,startA,endA,sampleA,ABA,DHFFCA,GTA,chrB,startB,endB;
                chrA=$1; startA=$2; endA=$3;
                sampleA=$4; ABA=$5; DHFFCA=$6; GTA=$7
                chrB=$8; startB=$9; endB=$10;
            }
        }; 
        END{print chrA,startA,endA,sampleA,ABA,DHFFCA,GTA,chrB,startB,endB;}'
}
export -f self_overlap

bed=$1
out_dir=$2
fai=$3 # ref genome index

[[ ! -d $out_dir/temp ]] && mkdir $out_dir/temp && mkdir $out_dir/temp/self_overlap

# separate bed by chromosome
cut -f1 $fai | gargs -p 62 "grep {} $bed | bedtools sort > $out_dir/temp/'{}'.bed"

# remove 0 byte files
find $out_dir/temp/ -type f -size 0 -delete

# self intersection by chromosome
find $out_dir/temp/ -name '*.bed' |
    gargs -p 62 "self_overlap {} > $out_dir/temp/self_overlap/\$(basename {})"
