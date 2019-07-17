#!/bin/bash

VCF=$1
bcftools view -i 'SVTYPE="DEL"' $VCF \
    | bcftools query -f '%CHROM\t%POS\t%INFO/END\t%SVLEN\n' \
    | cut -f4 \
    | (stats -b 50 -th 500 -t "SVLEN size distribution")
