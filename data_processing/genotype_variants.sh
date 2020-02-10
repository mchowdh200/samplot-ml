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
        --header-contigs) # text file containing vcf header contigs
            header_contigs=$2
            shift 2;;
        --out-dir)
            out_dir=$2
            shift 2;;
        # --s3-dest)
        #     s3_dist=$2
        #     shift 2;;
        --) # end argument parsing
            shift
            break;;
        -*|--*=) # unsupported flags
            echo "Error: Unsupported flag $1" >&2
            exit 1;;
    esac
done

s3_dest="s3://layerlab/samplot-ml/1kg/lumpy_sensitive/"

# download the alignments
cram=$(basename $cram_url)
crai=$cram.crai
sample=$(cut -d'.' -f1 <<<$cram)

[[ ! -d $out_dir/$sample ]] && mkdir $out_dir/$sample

# don't forget to set AWS environment variables
[[ ! -f $out_dir/$sample/$cram ]] && aws s3 cp $cram_url $out_dir/$sample 
[[ ! -f $out_dir/$sample/$crai ]] && aws s3 cp $cram_url.crai $out_dir/$sample 

# call sv's with smoove (really just want the split/disc. reads)
smoove call \
    --outdir $out_dir/$sample \
    --name $sample \
    --fasta $fasta \
    $out_dir/$sample/$cram

# modify the original lumpy command output by smoove to be more sensitive
tail -1 $out_dir/$sample/$sample-lumpy-cmd.sh \
    | awk '{$6=2; $8=2; print $0}' > $out_dir/$sample/lumpy_sensitive.sh

# execute the lumpy command and genotype with svtyper
bash $out_dir/$sample/lumpy_sensitive.sh \
    | bcftools annotate -h $header_contigs \
    | svtyper -B $out_dir/$sample/$cram \
              -T $fasta \
    | bcftools sort \
    | bgzip -c > $out_dir/$sample/$sample.genotyped.vcf.gz
tabix $out_dir/$sample/$sample.genotyped.vcf.gz
    
# copy vcf results to s3
aws s3 cp $out_dir/$sample/$sample.genotyped.vcf.gz $s3_dest
aws s3 cp $out_dir/$sample/$sample.genotyped.vcf.gz.tbi $s3_dest

# cleanup bam/crams
rm $out_dir/$sample/*.bam*
rm $out_dir/$sample/*.cram*

