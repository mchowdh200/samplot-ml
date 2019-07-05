#!/bin/bash

VCF_DIR=$1

bcftools view -i 'ml_het>0.5 || ml_alt>0.5' $VCF_DIR/ml.vcf.gz \
    | bgzip -c > $VCF_DIR/ml.gt_0.5.vcf.gz
tabix $VCF_DIR/ml.gt_0.5.vcf.gz

bcftools view -i 'ml_het>0.8 || ml_alt>0.8' $VCF_DIR/ml.vcf.gz \
    | bgzip -c > $VCF_DIR/ml.gt_0.8.vcf.gz
tabix $VCF_DIR/ml.gt_0.8.vcf.gz

bcftools view -i 'ml_het>0.9 || ml_alt>0.9' $VCF_DIR/ml.vcf.gz \
    | bgzip -c > $VCF_DIR/ml.gt_0.9.vcf.gz
tabix $VCF_DIR/ml.gt_0.9.vcf.gz

bcftools view -i 'ml_argmax > 0' $VCF_DIR/ml.vcf.gz \
    | bgzip -c > $VCF_DIR/ml.argmax.vcf.gz
tabix $VCF_DIR/ml.argmax.vcf.gz

bcftools view -i 'ml_het + ml_alt > ml_ref' $VCF_DIR/ml.vcf.gz \
    | bgzip -c > $VCF_DIR/ml.majority.vcf.gz
tabix $VCF_DIR/ml.majority.vcf.gz

