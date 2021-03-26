#!/bin/env bash
set -e

while (( "$#" )); do
    case "$1" in
        -c|--chrom)
            chrom=$2
            shift 2;;
        -s|--start)
            start=$2
            shift 2;;
        -e|--end)
            end=$2
            shift 2;;
        -n|--sample)
            sample=$2
            shift 2;;
        -g|--genotype)
            genotype=$2
            shift 2;;
        -m|--min-mqual)
            min_mq=$2
            shift 2;;
        -f|--fasta)
            fasta=$2
            shift 2;;
        -b|--bam)
            bam=$2
            shift 2;;
        -d|--delimiter)
            delimiter=$2
            shift 2;;
        -o|--outdir)
            outdir=$2
            shift 2;;
        --) # end argument parsing
            shift
            break;;
        -*|--*=) # unsupported flags
            echo "Error: Unsupported flag $1" >&2
            exit 1;;
    esac
done

# out=$outdir/${chrom}_${end}_${sample}_${genotype}.png
out=$outdir/$(echo "$chrom $end $sample $genotype.png" | tr ' ' $delimiter)
svlen=$(($end-$start))

if [[ $svlen -gt 5000 ]]; then
    samplot plot \
        --zoom 1000 \
        --chrom $chrom --start $start --end $end \
        --min_mqual  $min_mq \
        --sv_type DEL \
        --bams $bam \
        --reference $fasta \
        --output_file $out
else
    samplot plot \
        --chrom $chrom --start $start --end $end \
        --min_mqual  $min_mq \
        --sv_type DEL \
        --bams $bam \
        --reference $fasta \
        --output_file $out
fi



