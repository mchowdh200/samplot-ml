#!/bin/bash
set -e

# TODO combine beds

# arg parser -----------------------------------------------------------------
while (( "$#" )); do
    case "$1" in
        -f|--fasta)
            fasta=$2
            shift 2;;
        -r|--regions-bed)
            regions_bed=$2
            shift 2;;
        -l|--bam-list)
            bam_list=$2
            shift 2;;
        -d|--bam-dir)
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

cd $bam_dir

samples=(HG00514 HG00733 NA19240)
chromosomes=$(echo chr{{1..22},X})

for sample in ${samples[@]}; do
    for chr in $chromosomes; do
        # get bams
        bam_url=$(grep -P "(?=.*${chr}\.)(?=.*$sample)" $bam_list)
        bam_file=$(basename $bam_url)
        bai_url=$(dirname $bam_url)/$(basename $bam_url .bam).bai
        bai_file=$(basename $bai_url)

        [[ ! -f $bam_file ]] && wget $bam_url
        [[ ! -f $bai_file ]] && wget $bai_url

        # TODO should I upload this to s3 as well? (use & to keep in background)
        # make sure to export AWS keys before running this script

        # feed the regions from sample & chr to samplot
        # grep -P "(?=.*$chr\t)(?=.*$sample)" $regions_bed

        exit
    done
done

cd -
