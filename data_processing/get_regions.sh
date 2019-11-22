#!/bin/bash

set -e

vcf=$1
out_dir=$2

bcftools view -i 'SVTYPE="DEL"' $vcf \
    | bcftools query -f '%CHROM\t%POS\t%INFO/END[,%SAMPLE\t%GT]\n' \
    | python3 sample_del.py > $out_dir/1kg_regions.bed
 
# genomes related to the test set genomes
# also excluding chromosome X because I noticed a lot of false posotives
exclude=(NA12878 NA12891 NA12892 
         HG00512 HG00513 HG00514 
         HG00731 HG00732 HG00733
         NA19238 NA19239 NA19240
         X)

exclude=$(echo ${exclude[@]} | tr " " "|") # construct regex from array

# sample 100k from each category while excludeing test set genomes
grep -E -v "$exclude" $out_dir/1kg_regions.bed | grep "ref" | shuf -n 150000 > $out_dir/ref.sampled.bed
grep -E -v "$exclude" $out_dir/1kg_regions.bed | grep "het" | shuf -n 150000 > $out_dir/het.sampled.bed
grep -E -v "$exclude" $out_dir/1kg_regions.bed | grep "alt" | shuf -n 150000 > $out_dir/alt.sampled.bed

cat $out_dir/ref.sampled.bed \
    $out_dir/het.sampled.bed \
    $out_dir/alt.sampled.bed > $out_dir/training_regions.bed

