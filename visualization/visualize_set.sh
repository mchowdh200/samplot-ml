#!/bin/bash
set -eu

while (( "$#" )); do
    case "$1" in
        -v|--vcf)
            vcf=$2
            shift 2;;
        -s|--sample)
            sample=$2
            shift 2;;
        -p|--pred-bed)
            pred_bed=$2
            shift 2;;
        -i|--image-dir)
            image_dir=$2
            shift 2;;
        -o|--out-dir)
            out_dir=$2
            shift 2;;
        -n|--set-name)
            set_name=$2
            shift 2;;
    esac
done
[[ -z $vcf ]] && echo Missing argument -v/--vcf && exit 1
[[ -z $pred_bed ]] && echo Missing argument -p/--pred-bed && exit 1
[[ -z $image_dir ]] && echo Missing argument -i/--image-dir && exit 1
[[ -z $out_dir ]] && echo Missing argument -o/--out-dir && exit 1

# get a BED like format that has the region, prediction and dhffc
# 4, 5, 6 -> CNN prediction | 10 -> DHFFC
prediction_set=$(bcftools query -f '%CHROM\t%POS\t%INFO/END\t[%DHFFC]\n' $vcf \
    | intersectBed -wa -wb -f 1.0 -r -a $pred_bed -b stdin \
    | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$10}')

image_files=$(echo "$prediction_set" \
    | awk -v image_dir=$image_dir -v sample=$sample '{print image_dir"/"$1"_"$2"_"$3"_"sample"_DEL.png"}')

# construct an imagemagik montage command that collates images along
# with predictions/dhffc into single pdf.

# create the -label part of the command
montage_command=$(echo "$prediction_set" \
    | awk '{print "-pointsize 30 -label "$1"_"$2"_"$3":prediction=["$4","$5","$6"]:DHFFC="$7}')

# combine with the actual image file
montage_command=$(paste -d '\t' <(echo "$montage_command") <(echo "$image_files"))
montage_command="${montage_command} -mode concatenate -tile 1x1 $out_dir/${sample}_${set_name}.pdf"
# echo "$montage_command"
montage $montage_command

