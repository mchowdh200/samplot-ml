#!/bin/bash

set -eu

data_dir=$1
samples=(HG002 HG00514 HG00733 NA19240)

for sample in ${samples[@]}; do
    echo $sample ================================
    if [[ $sample == "HG002"  ]]; then
        truth_set=HG002_SVs_Tier1_v0.6.DEL.vcf.gz
    else
        truth_set=$sample.BIP-unified.del.fixed.vcf.gz
    fi

    echo Lumpy/SVTyper --------------------------------------------------------

    # * The SV2 vcf header does not contain a contig if it is not
    #   present in the VCF itself which causes problems with truvari.
    #   So we must first annotate the vcf with the full set of contigs
    #   from the reference genome.
    # * Also the INFO field for GENES incorrectly says it has just 1 item,
    #   but it could have many.
    bcftools view -i 'GT!="0/0" & SVTYPE="DEL"'\
        $data_dir/$sample/VCF/SV2-filter/sv2_genotypes/${sample}_smoove_sv2.vcf |
    grep -v contig |
    bcftools annotate -h $data_dir/$sample/VCF/header_contigs.txt |
    sed 's/##INFO=<ID=GENES,Number=1/##INFO=<ID=GENES,Number=./g' |
    bgzip -c > $data_dir/$sample/VCF/SV2-filter/sv2_genotypes/temp.vcf.gz
    tabix $data_dir/$sample/VCF/SV2-filter/sv2_genotypes/temp.vcf.gz

    bash truvari.sh \
        -c $data_dir/$sample/VCF/SV2-filter/sv2_genotypes/temp.vcf.gz \
        -b $data_dir/$sample/VCF/$truth_set \
        -o $data_dir/$sample/VCF/SV2-filter/truvari

    echo manta ----------------------------------------------------------------
    bcftools view -i 'GT!="0/0" & SVTYPE="DEL"'\
        $data_dir/$sample/VCF/SV2-filter-manta/sv2_genotypes/${sample}_manta_sv2.vcf |
    grep -v contig |
    bcftools annotate -h $data_dir/$sample/VCF/header_contigs.txt |
    sed 's/##INFO=<ID=GENES,Number=1/##INFO=<ID=GENES,Number=./g' |
    bgzip -c > $data_dir/$sample/VCF/SV2-filter-manta/sv2_genotypes/temp.vcf.gz
    tabix $data_dir/$sample/VCF/SV2-filter-manta/sv2_genotypes/temp.vcf.gz

    bash truvari.sh \
        -c $data_dir/$sample/VCF/SV2-filter-manta/sv2_genotypes/temp.vcf.gz \
        -b $data_dir/$sample/VCF/$truth_set \
        -o $data_dir/$sample/VCF/SV2-filter-manta/truvari

done
