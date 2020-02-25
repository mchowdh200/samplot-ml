#!/bin/bash

set -eu

vcf=$1

# bcftools view -i 'SVTYPE="DEL" & DHFFC>0.7 & GT[*]!="RR"' $vcf | # this will give us 0/1 & 1/1
bcftools view -i 'SVTYPE="DEL" & DHFFC>0.7' $vcf | # this will give us everything -- including 0/0
bcftools query -f '%CHROM\t%POS\t%INFO/END\t[%SAMPLE\t%AB\t%DHFFC\t%GT]\n' |
bedtools intersect -wa -a stdin -b BED/exclude.cnvnator_100bp.GRCh38.20170403.bed |
uniq
