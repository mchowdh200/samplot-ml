#!/bin/bash

VCF=data/VCF/ALL.wgs.integrated_sv_map_v1_GRCh38.20130502.svs.genotypes.vcf.gz
OUT=data/BED


bcftools view -i 'SVTYPE="DEL"' $VCF \
    | bcftools query -f '%CHROM\t%POS\t%INFO/END[\t%SAMPLE,%GT]\n'\
    | awk '$3-$2 > 500' \
    | python sample_del.py > $OUT/del.gt500.bed

bcftools view -i 'SVTYPE="DEL"' $VCF \
    | bcftools query -f '%CHROM\t%POS\t%INFO/END[\t%SAMPLE,%GT]\n'\
    | awk '$3-$2 < 500' \
    | python sample_del.py > $OUT/del.lt500.bed

bcftools view -i 'SVTYPE="DEL"' $VCF \
    | bcftools query -f '%CHROM\t%POS\t%INFO/END[\t%SAMPLE,%GT]\n'\
    | python sample_del.py > $OUT/del.sample.bed

