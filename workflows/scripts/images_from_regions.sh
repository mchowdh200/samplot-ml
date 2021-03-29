#!/bin/env bash
# given a set of bed regions, reference fasta, and bam/cram
# this script will output the regions to gargs which will generate
# individual samplot images using gen_img.sh
set -eu
while (( "$#" )); do
    case "$1" in
        -g|--gargs-bin)
            gargs_bin=$2
            shift 2;;
        -f|--fasta)
            fasta=$2
            shift 2;;
        -r|--regions)
            regions=$2
            shift 2;;
        -b|--bam)
            bam=$2
            shift 2;;
        -d|--delimiter) # delimiter between chrm, start, end, etc in filename
            delimiter=$2
            shift 2;;
        -o|--outdir)
            outdir=$2
            shift 2;;
        -p|--processes)
            processes=$2
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

# format of regions bed is:
# chrom start end svtype sample
echo $PWD
cat $regions | $gargs_bin -e -p $processes "bash scripts/gen_img.sh \\
    --chrom {0} --start {1} --end {2} --genotype {3} --sample {4} \\
    --min-mqual 10 --fasta $fasta --bam $bam \\
    --delimiter $delimiter --outdir $outdir"
