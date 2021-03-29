#!/bin/env bash

while (( "$#" )); do
    case "$1" in
        -b|--bam-bucket)
            bam_bucket=$2
            shift 2;;
        -f|--fasta)
            fasta=$2
            shift 2;;
        -i|--fai)
            fai=$2
            shift 2;;
        -v|--vcf)
            vcf=$2
            shift 2;;
        -d|--delimiter)
            vcf=$2
            shift 2;;
        -o|--outdir)
            vcf=$2
            shift 2;;
        --) # end argument parsing
            shift
            break;;
        -*|--*=) # unsupported flags
            echo "Error: Unsupported flag $1" >&2
            exit 1;;
    esac
done

echo "samples:"

# filenames
#bams=$(aws s3 ls $bam_bucket | grep -E '*.bam$' | sed -E 's/\s+/\t/g' | cut -f4)
# OR do the same thing with a local directory of bams
# TODO for now I'm just going to download bams locally and create this config
bams=$(ls $bam_bucket | grep -E '*.bam$')

### Sample bam pairs
# get the sample name from the bam header
for bam in $bams; do
    sample=$(samtools view -H $bam_bucket/$bam | grep SM |
             awk 'BEGIN{RS="\t"; FS=":"} /SM/ {print}')
    printf "  $sample: \"$bam_bucket/$bam\"\n"
done

### reference fasta and index
echo "fasta:"
printf "  data_source: \"s3\"\n"
printf "  file: \"$fasta\"\n"

echo "fai:"
printf "  data_source: \"s3\"\n"
printf "  file: \"$fai\"\n"

echo "vcf:"
printf "  data_source: \"s3\"\n"
printf "  file: \"$vcf\"\n"

echo "image_filename_delimiter: \"$delimiter\""
echo "image_filename_delimiter: \"$delimiter\""
echo "outdir: \"$outdir\""
