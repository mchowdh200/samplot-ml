#!/bin/env bash

function crop
{
    input_img=$1
    outdir=$2
    output_img=$outdir/$(basename $input_img)

    convert \
        -crop 2090x575+175+200 \
        -fill white \
        -draw "rectangle 0,30 500,50" \
        -draw "rectangle 600,0 700,50" \
        $input_img $output_img
    echo $output_img
}
export -f crop

while (( "$#" )); do
    case "$1" in
        -i|--imgdir)
            imgdir=$2
            shift 2;;
        -o|--outdir)
            outdir=$2
            shift 2;;
        -p|--processes)
            processes=$2
            shift 2;;
        -g|--gargs-bin)
            gargs_bin=$2
            shift 2;;
        --) # end argument parsing
            shift
            break;;
        -*|--*=) # unsupported flags
            echo "Error: Unsupported flag $1" >&2
            exit 1;;
    esac
done

[[ ! -d $outdir ]] && mkdir $outdir
find $imgdir -name '*.png' | $gargs_bin -e -p $processes "crop {0} $outdir"
