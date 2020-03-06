#!/bin/bash
# given a directory of vcfs, get all the DELs and turn them into a bed containing:
# chrom  start  end  sample AB DHFFC GT

set -eu

function query_dels() {
    bcftools view $1 -i 'SVTYPE="DEL"' |
    bcftools query -f '%CHROM\t%POS\t%INFO/END\t[%SAMPLE\t%AB\t%DHFFC\t%GT]\n' |
    uniq
}
export -f query_dels

vcf_dir=$1
n_tasks=$2

find "$vcf_dir" -name '*.vcf.gz' | gargs -p $n_tasks "query_dels {}"
