#!/bin/env bash
# given:
#     * input vcf
#     * sample name
# output:
#     * sample's DEL regions in bed format to stdout

# 1. get sample's part of vcf
# 2. get SVTYPE = DEL (only het/alt genotypes)
# 3. bcftools query into bed format with format:
#    %CHROM\t%POS\t%INFO/END\t%SVTYPE\t[%SAMPLE]\n

vcf=$1
sample=$2

# note use single & for within sample logic
bcftools view -s $sample -i 'SVTYPE="DEL"' $vcf |
    bcftools query -i 'GT!="0/0" & GT!="./."' \
        -f '%CHROM\t%POS\t%INFO/END\t%SVTYPE\t[%SAMPLE]\n'
