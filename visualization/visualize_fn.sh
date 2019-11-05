#!/bin/bash
set -e

while (( "$#" )); do
    case "$1" in
        -s|--sample)
            sample=$2
            shift 2;;
        -c|--CNN-vcf)
            cnn_vcf=$2
            shift 2;;
        -b|--baseline-vcf) # ie. the smoove, manta, etc. baseline callset
            baseline_vcf=$2
            shift 2;;
        -f|--fn-vcf)
            fn_vcf=$2
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
    esac
done
[[ -z $sample ]] && echo Missing argument -s/--sample && ERR=1
[[ -z $cnn_vcf ]] && echo Missing argument -c/--cnn_vcf && ERR=1
[[ -z $baseline_vcf ]] && echo Missing argument -b/--baseline_vcf && ERR=1
[[ -z $fn_vcf ]] && echo Missing argument -f/--fn_vcf && ERR=1
[[ -z $pred_bed ]] && echo Missing argument -p/--pred-bed && ERR=1
[[ -z $image_dir ]] && echo Missing argument -i/--image-dir && ERR=1
[[ -z $out_dir ]] && echo Missing argument -o/--out-dir && ERR=1
[[ ! -z $ERR ]] && exit 1

# new_false_negatives <- subtract the CNN tp calls from the base calls 
# prediction_set <- intersect $new_false_negatives with $pred_bed
#bedtools subtract -A -header -a $baseline_vcf -b $cnn_vcf \


prediction_set=$(subtractBed -A -header -a $baseline_vcf -b $cnn_vcf \
    | bcftools query -f '%CHROM\t%POS\t%INFO/END\t[%DHFFC]\n' \
    | intersectBed -wa -wb -f 1.0 -r -a $pred_bed -b stdin \
    | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$10}' \
    | intersectBed -wa -a stdin -b $fn_vcf \
    | sort -u # remove duplicates
)

# generate images with mapq threshold
FASTA="~/data/ref/GRCh38_full_analysis_set_plus_decoy_hla/GRCh38_full_analysis_set_plus_decoy_hla.fa"
echo "$prediction_set" | gargs -p 4 \
    "bash ../data_processing/gen_img.sh \\
        --chrom {0} --start {1} --end {2} --sample $sample --genotype del --min-mqual 60  \\
        --fasta $FASTA \\
        --bam-list ~/data/$sample/cram_list.txt \\
        --bam-dir ~/data/$sample/CRAM/ \\
        --out-dir ~/data/$sample/VCF/manta/CNN_10_21/visualization/"

# image_files=$(echo "$prediction_set" \
#     | awk -v image_dir=$image_dir '{print image_dir$1"_"$2"_"$3"*"}' \
#     | gargs "find {}")
image_files=$(ls ~/data/$sample/VCF/manta/CNN_10_21/visualization/*.png)

montage_command=$(echo "$prediction_set" \
    | awk '{print "-pointsize 30 -label "$1"_"$2"_"$3":prediction=["$4","$5","$6"]:DHFFC="$7}')

montage_command=$(paste -d '\t' <(echo "$montage_command") <(echo "$image_files"))
montage_command="${montage_command} -mode concatenate -tile 1x1 $out_dir/${sample}_pred_vis.pdf"

# echo $montage_command # hack: left out the "" so that the string is all on one line
montage $montage_command


