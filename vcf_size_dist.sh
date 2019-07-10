#!/bin/bash

ORIG=data/giab/VCF/HG002-smoove.genotyped.ml.vcf.gz
FP=duphold_comp/fp.vcf
FN=duphold_comp/fn.vcf

bcftools view -i 'SVTYPE="DEL"' $ORIG \
    | bcftools query -f '%CHROM\t%POS\t%INFO/END\t%SVLEN\n' > duphold_comp/orig.summary.bed

bcftools view -i 'SVTYPE="DEL"' $FP \
    | bcftools query -f '%CHROM\t%POS\t%INFO/END\t%SVLEN\n' > duphold_comp/fp.summary.bed

bcftools view -i 'SVTYPE="DEL"' $FN \
    | bcftools query -f '%CHROM\t%POS\t%INFO/END\t%SVLEN\n' > duphold_comp/fn.summary.bed

cut -f4 duphold_comp/orig.summary.bed \
    | (stats -a -b 50 -th 500 -t "Original VCF")

cut -f4 duphold_comp/fp.summary.bed \
    | (stats -a -b 50 -th 500 -t "False Positives")

cut -f4 duphold_comp/fn.summary.bed \
    | (stats -a -b 50 -th 500 -t "False Negatives")



