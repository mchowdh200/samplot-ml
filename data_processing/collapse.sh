#!/bin/bash

# take a sorted(!) bed file with the following format
# chr  start  end  sample  AB  DHFFC GT
# and intersect it with itself.
# prints two sets of intervals side by side
# if the chr start end in the first interval are the 
# same as the previous line's, then only change the end 
# of the second interval.
set -eu

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
