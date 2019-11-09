
input_bed=$1
out_dir=$2

# genomes related to the test set genomes
exclude=(NA12878 NA12891 NA12892 
         HG00512 HG00513 HG00514 
         HG00731 HG00732 HG00733
         NA19238 NA19239 NA19240) 
exclude=$(echo ${exclude[@]} | tr " " "|") # construct regex from array

# sample 100k from each category while excludeing test set genomes
grep -E -v "$exclude" $input_bed | grep "ref" | shuf -n 150000 > $out_dir/ref.sampled.bed
grep -E -v "$exclude" $input_bed | grep "het" | shuf -n 150000 > $out_dir/het.sampled.bed
grep -E -v "$exclude" $input_bed | grep "alt" | shuf -n 150000 > $out_dir/alt.sampled.bed

