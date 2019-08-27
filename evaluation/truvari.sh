#!/bin/bash

while (( "$#" )); do
    case "$1" in
        -c|--comp-vcf)
            comp_vcf=$2
            shift 2;;
        -b|--base-vcf) # ie truth set
            base_vcf=$2
            shift 2;;
        -r|--reference)
            reference=$2
            shift 2;;
        -o|--out-dir)
            out_dir=$2
            shift 2;;
    esac
done
[[ -z $comp_vcf ]] && echo Missing argument --comp-vcf && exit 1
[[ -z $base_vcf ]] && echo Missing argument --base-vcf && exit 1
[[ -z $reference ]] && echo Missing argument --reference && exit 1
[[ -z $out_dir ]] && echo Missing argument --out-dir && exit 1

# if the out dir already exists, truvari fails
[[ -d $out_dir ]] && rm -rf $out_dir

truvari \
    -b $base_vcf \
    -c $comp_vcf \
    -o $out_dir \
    --reference $reference \
    --pctsim 0 \
    --sizemax 1000000 \
    --sizemin 300 --sizefilt 270 \
    --pctov 0.6 \
    --refdist 20 \
    --noprog \
    --no-ref a \
    # --includebed ~/Repositories/samplot-ml/data/giab/BED/HG002_SVs_Tier1_v0.6.bed \
    # --sizemax 15000000 \
