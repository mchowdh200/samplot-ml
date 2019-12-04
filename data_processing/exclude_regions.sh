#!/bin/bash

function filter_dhffc
{
    vcf=$1
    out_dir=$2
    out_file=$(basename $vcf .vcf.gz).del.dhffc.gt0.7.vcf.gz
    echo $out_dir/$out_file
    
    if [[ ! -f $out_dir/$out_file  ]]; then
        bcftools view -i 'SVTYPE="DEL" & DHFFC > 0.7' $vcf \
            | bgzip -c > $out_dir/$out_file
        tabix $out_dir/$out_file
    fi
}
export -f filter_dhffc

function get_exclude_regions
{
    vcf=$1
    out_dir=$2
    out_file=$(basename $vcf .vcf.gz).bed
    exclude_bed=$3

    bcftools query -f '%CHROM\t%POS\t%INFO/END\t[%SAMPLE]\t [%DHFFC]\n' $vcf \
        | intersectBed -wa -a stdin -b $exclude_bed | sort | uniq > $out_dir/$out_file
}
export -f get_exclude_regions

# base_dir=$1
# duphold_dir=$2
# exclude_dir=$3
# exclude_bed=$4
# lazy...
base_dir="/home/murad/data/1kg/1kg_smoove/base_vcf"
duphold_dir="/home/murad/data/1kg/1kg_smoove/del.gt0.7"
exclude_dir="/home/murad/data/1kg/1kg_smoove/1kg_excluded"
exclude_bed="/home/murad/data/1kg/1kg_smoove/exclude.cnvnator_100bp.GRCh38.20170403.bed"

ls $base_dir/*.vcf.gz | gargs -p 8 "filter_dhffc {} $duphold_dir"
ls $duphold_dir/*.vcf.gz | gargs -p 8 "get_exclude_regions {} $exclude_dir $exclude_bed"
