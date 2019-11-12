#!/bin/bash

set -e

while (( "$#" )); do
    case "$1" in
        --fasta)
            fasta=$2
            shift 2;;
        --cram-url)
            cram_url=$2
            shift 2;;
        --crai-url)
            crai_url=$2
            shift 2;;
        --download-dir)
            download_dir=$2
            shift 2;;
        --out-dir)
            out_dir=$2
            shift 2;;
        --s3-dist)
            s3_dist=$2
            shift 2;;
        --) # end argument parsing
            shift
            break;;
        -*|--*=) # unsupported flags
            echo "Error: Unsupported flag $1" >&2
            exit 1;;
    esac
done

# download the alignments
cram=$(basename $cram_url)
crai=$(basename $crai_url)
sample=$(cut -d'.' -f1 $cram)
aws s3 cp $cram_url $download_dir # don't forget to set environment variables

# genotype with smoove. annotate with duphold
mkdir $out_dir/$sample
smoove call -p 2 \
    --outdir $out_dir/$sample \
    --name $sample \
    --fasta $fasta \
    --genotype \
    --duphold \
    $download_dir/$cram

# copy results to s3
aws s3 cp $out_dir/$sample/$sample-smoove.genotyped.vcf.gz $s3_dest
aws s3 cp $out_dir/$sample/$sample-smoove.genotyped.vcf.gz.csi $s3_dest

# delete local data
rm -r $out_dir/$sample/
rm $download_dir/$cram
rm $download_dir/$crai

