#!/bin/bash

set -eu

vcf=$1
# exclude_regions=$2

# bcftools view -i 'SVTYPE="DEL" & DHFFC>0.7 & GT[*]!="RR"' $vcf | # this will give us 0/1 & 1/1
# bcftools view -i 'SVTYPE="DEL" & DHFFC>0.7' $vcf | # this will give us everything -- including 0/0
bcftools view -i 'SVTYPE="DEL"' $vcf | # no DHFFC filtering
bcftools query -f '%CHROM\t%POS\t%INFO/END\t[%SAMPLE\t%AB\t%DHFFC\t%GT]\n' |
uniq
# bedtools intersect -wa \
#     -a stdin \
#     -b $exclude_regions |
