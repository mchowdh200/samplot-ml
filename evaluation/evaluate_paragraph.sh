#!/bin/bash

set -eu

samples=(HG002 HG002.2x250 HG00514 HG00733 NA19240)

for sample in ${samples[@]}; do
    echo ===
    echo === $sample ==========================================================
    echo ===

    if [[ $sample == "HG002" || $sample == "HG002.2x250"  ]]; then
        truth_set=~/data/HG002/VCF/HG002_SVs_Tier1_v0.6.DEL.vcf.gz
    else
        truth_set=~/data/$sample/VCF/$sample.BIP-unified.del.fixed.vcf.gz
    fi

    tabix -f ~/data/$sample/VCF/paragraph/smoove/genotypes.vcf.gz
    echo === smoove ===========================================================
    bash truvari.sh \
        -c ~/data/$sample/VCF/paragraph/smoove/genotypes.vcf.gz \
        -b $truth_set \
        -o ~/data/$sample/VCF/paragraph/smoove/truvari

    tabix -f ~/data/$sample/VCF/paragraph/manta/genotypes.vcf.gz
    echo === manta ===========================================================
    bash truvari.sh \
        -c ~/data/$sample/VCF/paragraph/manta/genotypes.vcf.gz \
        -b $truth_set \
        -o ~/data/$sample/VCF/paragraph/manta/truvari

done
