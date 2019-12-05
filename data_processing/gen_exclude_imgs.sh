#!/bin/bash
# Used on c5d ec2 instance to generate images from the exclude regions in a subset
# of 1000 genomes samples with sv's called by smoove and filtered by duphold
mkdir /mnt/local/data
mkdir /mnt/local/data/CRAM
mkdir /mnt/local/data/ref
mkdir /mnt/local/data/imgs

n_tasks=$1

# get ref
echo "Downloading reference genome..."
aws s3 cp s3://1000genomes/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa /mnt/local/data/ref/
aws s3 cp s3://1000genomes/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai /mnt/local/data/ref/

# get the 1kg CRAIs
echo "Downloading cram indices..."
cat ../data_listings/1kg_high_cov_crams_s3.txt \
    | gargs -p 48 "aws s3 cp {}.crai /mnt/local/data/CRAM/"

cat "../data_listings/1kg_smoove_genotyped_exclude.bed" \
    | gargs -p $n_tasks "bash ../data_processing/gen_img.sh \\
    --chrom {0} --start {1} --end {2} --sample {3} --genotype REF \\
    --min-mqual 5 \\
    --fasta /mnt/local/data/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa \\
    --bam-list ~/samplot-ml/data_listings/1kg_high_cov_crams_s3.txt \\
    --bam-dir /mnt/local/data/CRAM \\
    --out-dir /mnt/local/data/imgs"

aws s3 cp --recursive /mnt/local/data/imgs s3://layerlab/samplot-ml/1kg/exclude_ref_regions/

