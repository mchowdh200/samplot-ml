#!/bin/bash
set -eu

OUT_DIR=analysis/comparison

# directory of A and B vcf
A=$1
B=$2
PRED_FILE=$3


bedtools intersect \
    -a <(bcftools view -i 'SVTYPE="DEL"' $A/fp.vcf \
       | bcftools query -f '%CHROM\t%POS\t%INFO/END[\t%SAMPLE,%GT]\n') \
    -b <(bcftools view -i 'SVTYPE="DEL"' $B/fp.vcf \
       | bcftools query -f '%CHROM\t%POS\t%INFO/END[\t%SAMPLE,%GT]\n') \
       | awk '{print "data/giab/crop/"$1"_"$2"_"$3"_DEL.png"}' > $OUT_DIR/fp.intersect.images.txt
bedtools intersect \
    -a <(bcftools view -i 'SVTYPE="DEL"' $A/fp.vcf \
       | bcftools query -f '%CHROM\t%POS\t%INFO/END[\t%SAMPLE,%GT]\n') \
    -b <(bcftools view -i 'SVTYPE="DEL"' $B/fp.vcf \
       | bcftools query -f '%CHROM\t%POS\t%INFO/END[\t%SAMPLE,%GT]\n') > $OUT_DIR/fp.intersect.bed

bedtools intersect -wa -r -f 1.0 -a $PRED_FILE -b $OUT_DIR/fp.intersect.bed \
    | sort -k1,1V -k2,2n -k3,3n > $OUT_DIR/fp.pred.intersect.bed

MONTAGE_COMMAND=$(paste -d '\t' $OUT_DIR/fp.pred.intersect.bed $OUT_DIR/fp.intersect.images.txt \
    | awk '{print "-pointsize 30 -label "$1"_"$2"_"$3":["$4","$5","$6"] "$7" "}')
MONTAGE_COMMAND="${MONTAGE_COMMAND} -mode concatenate -tile 1x1 $OUT_DIR/fp.intersect.pdf"
montage $MONTAGE_COMMAND


bedtools subtract -A \
    -a <(bcftools view -i 'SVTYPE="DEL"' $A/fp.vcf \
       | bcftools query -f '%CHROM\t%POS\t%INFO/END[\t%SAMPLE,%GT]\n') \
    -b <(bcftools view -i 'SVTYPE="DEL"' $B/fp.vcf \
       | bcftools query -f '%CHROM\t%POS\t%INFO/END[\t%SAMPLE,%GT]\n') \
       | awk '{print "data/giab/crop/"$1"_"$2"_"$3"_DEL.png"}' > $OUT_DIR/fp.difference.images.txt
bedtools subtract -A\
    -a <(bcftools view -i 'SVTYPE="DEL"' $A/fp.vcf \
       | bcftools query -f '%CHROM\t%POS\t%INFO/END[\t%SAMPLE,%GT]\n') \
    -b <(bcftools view -i 'SVTYPE="DEL"' $B/fp.vcf \
       | bcftools query -f '%CHROM\t%POS\t%INFO/END[\t%SAMPLE,%GT]\n') > $OUT_DIR/fp.difference.bed

if [ ! -s $OUT_DIR/fp.difference.bed ]; then
    bedtools subtract -A \
        -a <(bcftools view -i 'SVTYPE="DEL"' $B/fp.vcf \
           | bcftools query -f '%CHROM\t%POS\t%INFO/END[\t%SAMPLE,%GT]\n') \
        -b <(bcftools view -i 'SVTYPE="DEL"' $A/fp.vcf \
           | bcftools query -f '%CHROM\t%POS\t%INFO/END[\t%SAMPLE,%GT]\n') \
           | awk '{print "data/giab/crop/"$1"_"$2"_"$3"_DEL.png"}' > $OUT_DIR/fp.difference.images.txt
    bedtools subtract -A \
        -a <(bcftools view -i 'SVTYPE="DEL"' $B/fp.vcf \
           | bcftools query -f '%CHROM\t%POS\t%INFO/END[\t%SAMPLE,%GT]\n') \
        -b <(bcftools view -i 'SVTYPE="DEL"' $A/fp.vcf \
           | bcftools query -f '%CHROM\t%POS\t%INFO/END[\t%SAMPLE,%GT]\n') > $OUT_DIR/fp.difference.bed
fi

bedtools intersect -r -f 1.0 -wa -a $PRED_FILE -b $OUT_DIR/fp.difference.bed \
    | sort -k1,1V -k2,2n -k3,3n > $OUT_DIR/fp.pred.difference.bed

MONTAGE_COMMAND=$(paste -d '\t' $OUT_DIR/fp.pred.difference.bed $OUT_DIR/fp.difference.images.txt \
    | awk '{print "-pointsize 30 -label "$1"_"$2"_"$3":["$4","$5","$6"] "$7" "}')
MONTAGE_COMMAND="${MONTAGE_COMMAND} -mode concatenate -tile 1x1 $OUT_DIR/fp.difference.pdf"
montage $MONTAGE_COMMAND


