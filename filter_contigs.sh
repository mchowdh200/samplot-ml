#!/bin/bash
VCF=$1
list="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y"
contigs=""
for i in $list; do
    contigs="$contigs,chr$i"
done
contigs=${contigs:1} # remove the leading comma

bcftools view -r $contigs $1
