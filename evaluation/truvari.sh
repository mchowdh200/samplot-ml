#!/bin/bash

while (( "$#" )); do
    case "$1" in
        -h|--help)
            HELP=true
            shift 1;;
        -c|--comp-vcf)
            comp_vcf=$2
            shift 2;;
        -b|--base-vcf) # ie truth set
            base_vcf=$2
            shift 2;;
        -o|--out-dir)
            out_dir=$2
            shift 2;;
        --) # end argument parsing
            shift
            break;;
        -*|--*=) # unsupported flags
            echo "Error: Unsupported flag $1" >&2
            exit 1;;
    esac
done

HELP_MESSAGE="
-h/--help:              Display help message
-c/--comp-vcf:          Comparision (call set) vcf
-b/--base-vcf:          Base (truth set) vcf
-o/--out-dir:           Output results directory
"
[[ ! -z $HELP ]] && echo "$HELP_MESSAGE" && exit
[[ -z $comp_vcf ]] && printf "\nMissing argument --comp-vcf\n" && echo "$HELP_MESSAGE" && exit 1
[[ -z $base_vcf ]] && printf "\nMissing argument --base-vcf\n" && echo "$HELP_MESSAGE" && exit 1
[[ -z $out_dir ]] && printf "\nMissing argument --out-dir\n" && echo "$HELP_MESSAGE" && exit 1

# if the out dir already exists, truvari fails
[[ -d $out_dir ]] && rm -rf $out_dir

truvari \
    -b $base_vcf \
    -c $comp_vcf \
    -o $out_dir \
    --pctsim 0 \
    --sizemax 1000000 \
    --pctov 0.6 \
    --refdist 20 \
    --no-ref a \
    --sizemin 300 --sizefilt 270 \
    # --includebed ~/Repositories/samplot-ml/data/giab/BED/HG002_SVs_Tier1_v0.6.bed \
