#!/bin/bash

set -e

in_dir=$1
out_dir=$2

ls $in_dir/*.vcf.gz | gargs -p 4 "
    bcftools view -i 'SVTYPE=\"DEL\"' {} \\
        | bcftools query -f '%CHROM\t%POS\t%INFO/END\t[%SAMPLE\t%GT\t%DHFFC]\n'\\
        | sed -e 's/0|0/ref/' -e 's/1|0/het/' -e 's/0|1/het/' -e 's/1|1/alt/'" > $out_dir/1kg_regions.bed

# genomes related to the test set genome
# also excluding chromosome X because I noticed a lot of false positives
exclude=(NA12878 NA12891 NA12892
         HG00512 HG00513 HG00514 
         HG00731 HG00732 HG00733 
         NA19238 NA19239 NA19240
         X)

exclude=$(echo ${exclude[@]} | tr " " "|") # construct regex from array

# sample 100k from each category while excludeing test set genomes 
grep -E -v "$exclude" $out_dir/1kg_regions.bed | grep "ref" | shuf -n 150000 > $out_dir/ref.bed
grep -E -v "$exclude" $out_dir/1kg_regions.bed \
    | grep "het" \
    | awk '{ if($6 <= 0.7){print $0;} }' \
    | shuf -n 150000 > $out_dir/het.bed
grep -E -v "$exclude" $out_dir/1kg_regions.bed \
    | grep "alt" \
    | awk '{ if($6 <= 0.7){print $0;} }' \
    | shuf -n 150000 > $out_dir/alt.bed

cat $out_dir/ref.bed \
    $out_dir/het.bed \
    $out_dir/alt.bed > $out_dir/training_regions.bed

