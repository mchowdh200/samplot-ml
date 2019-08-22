#!/bin/bash

set -e

while (( "$#" )); do
    case "$1" in
        -s|--smoove-vcf)
            smoove_vcf=$2
            shift 2;;
        -p|--predictions-bed)
            predictions_bed=$2
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

[[ -z $smoove_vcf ]] && echo Missing argument -s/--smoove-vcf && exit 1
[[ -z $predictions_bed ]] && echo Missing argument -p/--predictions-bed && exit 1
[[ -z $base_vcf ]] && echo Missing argument -b/--base-vcf && exit 1
[[ -z $reference ]] && echo Missing argument -b/--reference-vcf && exit 1
[[ -z $out_dir ]] && echo Missing argument -o/--out-dir && exit 1


[[ -f $out_dir/pr_curve.txt ]] && rm $out_dir/pr_curve.txt

for threshold in $(seq 0 0.05 1); do

    # create vcf for each threshold
    python3 filter_threshold.py $smoove_vcf $predictions_bed $threshold \
        | bgzip  -c > $out_dir/$threshold.vcf.gz
    tabix $out_dir/$threshold.vcf.gz

    # evaluate generated vcf using truvari
    bash truvari.sh \
        --base-vcf $base_vcf \
        --comp-vcf $out_dir/$threshold.vcf.gz \
        --reference $reference \
        --out-dir $out_dir/$threshold.truvari

    # get precision, recall, F1 from the log file output
    P=$(echo $(grep '"precision":' $out_dir/$threshold.truvari/log.txt) \
        | cut -d ' ' -f2 | cut -d ',' -f1)
    R=$(echo $(grep '"recall":' $out_dir/$threshold.truvari/log.txt) \
        | cut -d ' ' -f2 | cut -d ',' -f1)
    F1=$(echo $(grep '"f1":' $out_dir/$threshold.truvari/log.txt) \
        | cut -d ' ' -f2 | cut -d ',' -f1)

    echo $threshold $P $R $F1 >> $out_dir/pr_curve.txt
done
