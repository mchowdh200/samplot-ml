#!/bin/bash
set -e

# arg parser -----------------------------------------------------------------
while (( "$#" )); do
    case "$1" in
        -f|--fasta)
            fasta=$2
            shift 2;;
        -r|--regions-bed) # bed with format: $chr $start $end $sample $genotype
            regions_bed=$2
            shift 2;;
        -l|--bam-list) # should be s3 paths
            bam_list=$2
            shift 2;;
        -d|--bam-dir) # place to download bams
            bam_dir=$2
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
[[ -z $regions_bed ]] && echo Missing argument --regions-bed && exit 1
[[ -z $bam_list ]] && echo Missing argument --bam-list && exit 1
[[ -z $bam_dir ]] && echo Missing argument --bam-dir && exit 1
[[ -z $out_dir ]] && echo Missing argument --out-dir && exit 1


# -----------------------------------------------------------------------------
# pass parameters to samplot command.  If the region is > 1Mb then use the zoom flag
# -----------------------------------------------------------------------------
function samplot_cmd {
    local ch=$1
    local st=$2
    local en=$3
    local sa=$4
    local gt=$5
    local bam=$6
    local out_dir=$7
    out=$out_dir/${ch}_${st}_${en}_${sa}_${gt}.png
    echo $out
    if [[ ! -f $OUT ]]; then
        if [[ $(( $en - $st )) -gt 1000000 ]]; then
            samplot.py -c $ch -s $st -e $en -t $gt -b $bam -o $out --zoom 500
        else
            samplot.py -c $ch -s $st -e $en -t $gt -b $bam -o $out
        fi
    fi
}
export -f samplot_cmd

# -----------------------------------------------------------------------------
# given a sample and chromosome, download (if necessary) the bam/bai
# and generate images from regions_bed subset containing the sample/chromosome
# -----------------------------------------------------------------------------
function get_img {
    local chr=$1
    local sample=$2
    local regions_bed=$3
    local bam_list=$4
    local bam_dir=$5
    local out_dir=$6

    # echo ----------------------------------------------------------------------
    # echo $chr $sample $bam_list $out_dir

    # get bam/bai if they're not present
    bam_url=$(grep -P "(?=.*${chr}\.)(?=.*$sample)" $bam_list)
    bam_file=$(basename $bam_url)
    bai_url=$(dirname $bam_url)/$(basename $bam_url .bam).bai
    bai_file=$(basename $bai_url)

    # TODO make sure to export AWS keys before running this script
    # echo "aws s3 cp $bam_url $bam_dir"
    # exit
    [[ ! -f $bam_file ]] && aws s3 cp $bam_url $bam_dir
    [[ ! -f $bai_file ]] && aws s3 cp $bai_url $bam_dir

    # feed the regions from sample & chr to samplot
    grep -P "(?=.*$chr\t)(?=.*$sample)" $regions_bed | gargs \
        "samplot_cmd {0} {1} {2} {3} {4} $bam_file $out_dir"
    # rm $bam_file
    # rm $bai_file
}
export -f get_img






cd $bam_dir

samples=(HG00514 HG00733 NA19240)
chromosomes=($(echo chr{{1..22},X}))

for sample in ${samples[@]}; do
    printf '%s\n' "${chromosomes[@]}" | head -4 | gargs -p 4 "get_img {} $sample $regions_bed $bam_list $bam_dir $out_dir"

    exit
    # TODO instead of for loop just printf '%s\n' "${chromosome[@]}" array into gargs
    # TODO write a function containing the loop body
    # for chr in $chromosomes; do
    #     # get bams
    #     bam_url=$(grep -P "(?=.*${chr}\.)(?=.*$sample)" $bam_list)
    #     bam_file=$(basename $bam_url)
    #     bai_url=$(dirname $bam_url)/$(basename $bam_url .bam).bai
    #     bai_file=$(basename $bai_url)

    #     # TODO make sure to export AWS keys before running this script
    #     [[ ! -f $bam_file ]] && aws s3 cp $bam_url $bam_dir
    #     [[ ! -f $bai_file ]] && aws s3 cp $bai_url $bam_dir


    #     # feed the regions from sample & chr to samplot
    #     # grep -P "(?=.*$chr\t)(?=.*$sample)" $regions_bed | gargs -p 1 "echo {}"
    #     grep -P "(?=.*$chr\t)(?=.*$sample)" $regions_bed | gargs \
    #         "samplot_cmd {0} {1} {2} {3} {4} $bam_file $out_dir"
    #     exit
    # done
done

cd -
