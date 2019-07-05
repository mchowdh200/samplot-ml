#!/bin/bash
VCF_DIR=data/giab/VCF
BED_DIR=data/giab/BED
lumpy_vcf=$VCF_DIR/HG002-smoove.genotyped.vcf.gz
# fasta=/scratch/Shares/layer/ref/hs37-1kg/human_g1k_v37.20.fasta
fasta=data/FASTA/human_g1k_v37.20.fasta

#truvari –s 300 –S 270 –b HG002_SVs_Tier1_v0.6.DEL.vcf.gz –c $lumpy_vcf –o eval-no-support –passonly –pctsim = 0 –r 20 –giabreport –f $fasta –no-ref –includebed HG002_SVs_Tier1_v0.6.bed –O 0.6


out_dir=duphold_comp
rm -rf $out_dir

#truvari \
#~/src/truvari.brentp/truvari.py \

# vcf=$VCF_DIR/HG002-smoove.genotyped.ml.gt_0.9.vcf.gz
#vcf=$VCF_DIR/HG002-smoove.genotyped.ml.gt_0.8.vcf.gz
# vcf=$VCF_DIR/HG002-smoove.genotyped.DEL.DHFFC_lt_0.7.vcf.gz 
# vcf=$VCF_DIR/HG002-smoove.genotyped.ml.gt_0.5.vcf.gz
# vcf=$VCF_DIR/HG002-smoove.genotyped.ml.argmax.vcf.gz
# vcf=$VCF_DIR/HG002-smoove.genotyped.ml.majority.vcf.gz
# vcf=$VCF_DIR/HG002-smoove.genotyped.ml.vcf.gz
vcf=$1

    # -b $lumpy_vcf \
truvari \
    -b $VCF_DIR/HG002_SVs_Tier1_v0.6.DEL.vcf.gz \
    -c $vcf \
    -o $out_dir \
    --reference $fasta \
    --pctsim 0 \
    --sizemax 15000000 \
    --sizemin 300 --sizefilt 270 \
    --includebed $BED_DIR/HG002_SVs_Tier1_v0.6.bed \
    --pctov 0.6 \
    --giabreport \
    --refdist 20 \
    --no-ref a

    #--passonly \


#truvari.py --sizemax 15000000 -s 300 -S 270 -b HG002_SVs_Tier1_v0.6.DEL.vcf.gz -c $dupholded_vcf -o $out \
   #--passonly --pctsim=0  -r 20 --giabreport -f $fasta --no-ref --includebed HG002_SVs_Tier1_v0.6.bed -O 0.6
#

