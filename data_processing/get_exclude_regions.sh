#!/bin/bash

# smoove_vcf=$1
# exclude_bed=$2
while (( "$#" )); do
    case "$1" in
        -c|--call-vcf)
            call_vcf=$2
            shift 2;;
        -e|--exclude-bed)
            exclude_bed=$2
            shift 2;;
        -f|--fasta)
            fasta=$2
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
            echo "Error: unsupported flag $1" >&2
            exit 1;;
    esac
done
[[ -z $call_vcf  ]] && echo Missing argument --call-vcf && err=1
[[ -z $exclude_bed  ]] && echo Missing argument --exclude-bed && err=1
[[ -z $fasta  ]] && echo Missing argument --fasta && err=1
[[ -z $bam_dir  ]] && echo Missing argument --bam-dir && err=1
[[ -z $out_dir  ]] && echo Missing argument --out-dir && err=1
[[ ! -z $err ]] && exit 1


bcftools view -i 'SVTYPE="DEL"' $call_vcf \
    | bcftools query -f '%CHROM\t%POS\t%INFO/END\t[%SAMPLE]\t[%GT]\n' \
    | intersectBed -a stdin -b $exclude_bed -wa \
    | sort -u \
    | sed -e 's/0\/1/het/' -e 's/1\/1/alt/' \
    | gargs -p 2 \
        "bash gen_img.sh \\
            --chrom {0} --start {1} --end {2} --sample {3} --genotype {4} \\
            --min-mqual 5 \\
            --fasta $fasta \\
            --bam-dir $bam_dir \\
            --out-dir $out_dir"
