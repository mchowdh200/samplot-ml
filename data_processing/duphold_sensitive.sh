#!/bin/bash
set -eu

function svtyper_sample {
    while (( "$#" )); do
        case "$1" in
            --vcf-dir)
                local vcf_dir=$2
                shift 2;;
            --cram-dir)
                local cram_dir=$2
                shift 2;;
            --source-vcf)
                local source_vcf=$2
                shift 2;;
            --fasta)
                local fasta=$2
                shift 2;;
            --out-dir)
                local out_dir=$2
                shift 2;;
            --s3-dest)
                local s3_dest=$2
                shift 2;;
            --) # end argument parsing
                shift
                break;;
            -*|--*=) # unsupported flags
                echo "Error: Unsupported flag $1" >&2
                exit 1;;
        esac
    done

    local sample=$(basename $source_vcf | cut -d'.' -f1)
    echo $sample

    # check if this sample has already been processed
    if ! grep -q $sample <<<$(aws s3 ls $s3_dest); then
        # get cram/index
        local cram_url=$(grep $sample ../data_listings/1kg_high_cov_crams_s3.txt)
        local cram=$cram_dir/$(basename $cram_url)
        if [[ ! -f $cram ]]; then
            aws s3 cp $cram_url $cram_dir
            aws s3 cp $cram_url.crai $cram_dir
        fi

        # get vcf/index
        local vcf=$vcf_dir/$(basename $source_vcf)
        if [[ ! -f $vcf ]]; then
            aws s3 cp $source_vcf $vcf_dir
            aws s3 cp $source_vcf.tbi $vcf_dir
        fi
        
        # run svtyper on the vcf/cram
        bcftools view -i 'SVTYPE="DEL"' $vcf |
        duphold -v $vcf -f $fasta -b $cram |
        bgzip -c > $out_dir/$sample.sensitive.duphold.vcf.gz
        tabix $out_dir/$sample.sensitive.duphold.vcf.gz

        # upload results 
        # aws s3 cp $out_dir/$sample.sensitive.duphold.vcf.gz $s3_dest
        # aws s3 cp $out_dir/$sample.sensitive.duphold.vcf.gz.tbi $s3_dest

        # and clean up
        # rm $cram $cram.crai
        # rm $vcf $vcf.tbi
        # rm $out_dir/$sample.svtyper.duphold.vcf.gz
        # rm $out_dir/$sample.svtyper.duphold.vcf.gz.tbi

    fi
}
export -f svtyper_sample

if [ ! -d /mnt/local/data ]; then
    mkdir /mnt/local/data
    mkdir /mnt/local/data/CRAM
    mkdir /mnt/local/data/ref
    mkdir /mnt/local/data/VCF
    mkdir /mnt/local/data/output
fi

cram_dir=/mnt/local/data/CRAM
ref_dir=/mnt/local/data/ref
vcf_dir=/mnt/local/data/VCF
out_dir=/mnt/local/data/output

n_tasks=$1

# echo "Downloading reference genome..."
if [[ ! -f $ref_dir/GRCh38_full_analysis_set_plus_decoy_hla.fa ]]; then
    aws s3 cp s3://1000genomes/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa $ref_dir
    aws s3 cp s3://1000genomes/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai $ref_dir
fi
fasta=$ref_dir/GRCh38_full_analysis_set_plus_decoy_hla.fa

# for each vcf in lumpy_sensitive set
s3_loc="s3://layerlab/samplot-ml/1kg/lumpy_sensitive"
s3_dest="s3://layerlab/samplot-ml/1kg/lumpy_sensitive/duphold/"
vcf_list=$(aws s3 ls "$s3_loc/" |
    grep -E 'vcf.gz$' |
    sed -E 's/ +/\t/g' |
    cut -f4)

echo "$vcf_list" | head -1 |
    gargs -p $n_tasks \
        "svtyper_sample --source-vcf $s3_loc/{} \\
                        --vcf-dir $vcf_dir \\
                        --cram-dir $cram_dir \\
                        --fasta $fasta \\
                        --out-dir $out_dir \\
                        --s3-dest $s3_dest"

