#!/bin/bash

set -e

vcf=$1
out_dir=$2

bcftools view -i 'SVTYPE="DEL"' $vcf \
    | bcftools query -f '%CHROM\t%POS\t%INFO/END[,%SAMPLE\t%GT]\n' \
    | python3 sample_del.py > $out_dir/1kg_regions.bed

rm *.out *.err

    # TODO
    # get all 0/1 and 1/1 genotypes --> del.bed
    # get sample of 0/0 genotypes --> ref.bed
    # combine both to get

