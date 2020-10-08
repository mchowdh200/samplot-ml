#!/bin/bash

# setup bio perl and other dependencies for cnv-jacg package
# 1) perl stuff
sudo apt install bioperl -y
perl -MCPAN -e 'install Bio::DB::HTS'
perl -MCPAN -e 'Statistics::Basic'

# 2) R stuff
sudo apt install r-base -y
sudo su - -c "R -e \"install.packages('randomForest', repos='http://cran.rstudio.com/')\""

# 3) cnv-jacg
wget https://github.com/sunnyzxh/CNV-JACG/archive/v1.1.zip
unzip v1.1.zip


# pull data from s3 -----------------------------------------------------------
data_dir=/mnt/local/data
mkdir $data_dir

mkdir $data_dir/ref/
aws s3 cp --recursive \
    s3://layerlabcu/ref/genomes/GRCh38_full_analysis_set_plus_decoy_hla/ \
    $data_dir/ref/GRCh38_full_analysis_set_plus_decoy_hla/
aws s3 cp --recursive \
    s3://layerlabcu/ref/genomes/g1k_v37_decoy/ \
    $data_dir/ref/g1k_v37_decoy/
aws s3 cp --recursive \
    s3://layerlabcu/ref/genomes/g1k_v37_decoy/ \
    $data_dir/ref/g1k_v37_decoy/


mkdir $data_dir/HG002
aws s3 cp --recursive s3://layerlabcu/samplot-ml/HG002/smoove/ $data_dir/HG002/smoove/
aws s3 cp --recursive s3://layerlabcu/samplot-ml/HG002/manta/ $data_dir/HG002/manta/
aws s3 cp --recursive s3://layerlabcu/samplot-ml/HG002/smoove.2x250/ $data_dir/HG002/smoove.2x250/
aws s3 cp --recursive s3://layerlabcu/samplot-ml/HG002/manta.2x250/ $data_dir/HG002/manta.2x250/
aws s3 cp --recursive s3://layerlabcu/samplot-ml/HG002/BAM/ $data_dir/HG002/BAM/

mkdir $data_dir/chaisson
aws s3 cp --recursive s3://layerlabcu/samplot-ml/chaisson/VCF/smoove/ $data_dir/chaisson/smoove/
aws s3 cp --recursive s3://layerlabcu/samplot-ml/chaisson/VCF/manta/ $data_dir/chaisson/manta/
aws s3 cp --recursive s3://layerlabcu/samplot-ml/chaisson/BAM/ $data_dir/chaisson/BAM/


# create the cnv bed files
# hg002
bcftools query -f '%CHROM\t%POS\t%INFO/END\tDEL\n' \
    $data_dir/HG002/smoove/HG002-smoove.genotyped.del.vcf.gz \
    > $data_dir/HG002/smoove/HG002.smoove.bed 
bcftools query -f '%CHROM\t%POS\t%INFO/END\tDEL\n' \
    $data_dir/HG002/smoove.2x250/HG002-smoove.genotyped.del.vcf.gz \
    > $data_dir/HG002/smoove.2x250/HG002.smoove.2x250.bed 
bcftools query -f '%CHROM\t%POS\t%INFO/END\tDEL\n' \
    $data_dir/HG002/manta/diploidSV.del.duphold.vcf.gz \
    > $data_dir/HG002/manta/HG002.manta.bed 
bcftools query -f '%CHROM\t%POS\t%INFO/END\tDEL\n' \
    $data_dir/HG002/manta.2x250/diploidSV.del.duphold.vcf.gz \
    > $data_dir/HG002/manta.2x250/HG002.manta.2x250.bed 

# hg00514
bcftools query -f '%CHROM\t%POS\t%INFO/END\tDEL\n' \
    $data_dir/chaisson/smoove/HG00514-smoove.genotyped.del.vcf.gz \
    > $data_dir/chaisson/smoove/HG00514.smoove.bed 
bcftools query -f '%CHROM\t%POS\t%INFO/END\tDEL\n' \
    $data_dir/chaisson/manta/HG00514.manta.vcf.gz \
    > $data_dir/chaisson/manta/HG00514.manta.bed 

# hg00733
bcftools query -f '%CHROM\t%POS\t%INFO/END\tDEL\n' \
    $data_dir/chaisson/smoove/HG00733-smoove.genotyped.del.vcf.gz \
    > $data_dir/chaisson/smoove/HG00733.smoove.bed 
bcftools query -f '%CHROM\t%POS\t%INFO/END\tDEL\n' \
    $data_dir/chaisson/manta/HG00733.manta.vcf.gz \
    > $data_dir/chaisson/manta/HG00733.manta.bed 

# na19240
bcftools query -f '%CHROM\t%POS\t%INFO/END\tDEL\n' \
    $data_dir/chaisson/smoove/NA19240-smoove.genotyped.del.vcf.gz \
    > $data_dir/chaisson/smoove/NA19240.smoove.bed 
bcftools query -f '%CHROM\t%POS\t%INFO/END\tDEL\n' \
    $data_dir/chaisson/manta/NA19240.manta.vcf.gz \
    > $data_dir/chaisson/manta/NA19240.manta.bed 

# run cnv-jacg ----------------------------------------------------------------

# HG002
perl ~/CNV-JACG-1.1/bin/CNV-JACG.pl \
    -p $data_dir/HG002/smoove/HG002.smoove.bed \
    -b $data_dir/HG002/BAM/hg002.bam \
    -r $data_dir/ref/g1k_v37_decoy/g1k_v37_decoy.fa \
    -o $data_dir/HG002/smoove/results/
perl ~/CNV-JACG-1.1/bin/CNV-JACG.pl \
    -p $data_dir/HG002/manta/HG002.manta.bed \
    -b $data_dir/HG002/BAM/hg002.bam \
    -r $data_dir/ref/g1k_v37_decoy/g1k_v37_decoy.fa\
    -o $data_dir/HG002/manta/results/
perl ~/CNV-JACG-1.1/bin/CNV-JACG.pl \
    -p $data_dir/HG002/smoove.2x250/HG002.smoove.2x250.bed \
    -b $data_dir/HG002/BAM/HG002.2x250.sorted.bam \
    -r $data_dir/ref/hs37d5/hs37d5.fa \
    -o $data_dir/HG002/smoove.2x250/results/
perl ~/CNV-JACG-1.1/bin/CNV-JACG.pl \
    -p $data_dir/HG002/manta.2x250/HG002.manta.2x250.bed \
    -b $data_dir/HG002/BAM/HG002.2x250.sorted.bam \
    -r $data_dir/ref/hs37d5/hs37d5.fa \
    -o $data_dir/HG002/manta.2x250/results/

# HG00514
perl ~/CNV-JACG-1.1/bin/CNV-JACG.pl \
    -p $data_dir/chaisson/smoove/HG00514.smoove.bed \
    -b $data_dir/chaisson/BAM/HG00514.alt_bwamem_GRCh38DH.20150715.CHS.high_coverage.cram \
    -r $data_dir/ref/GRCh38_full_analysis_set_plus_decoy_hla/GRCh38_full_analysis_set_plus_decoy_hla.fa \
    -o $data_dir/chaisson/HG00514/smoove/results/
perl ~/CNV-JACG-1.1/bin/CNV-JACG.pl \
    -p $data_dir/chaisson/manta/HG00514.manta.bed \
    -b $data_dir/chaisson/BAM/HG00514.alt_bwamem_GRCh38DH.20150715.CHS.high_coverage.cram \
    -r $data_dir/ref/GRCh38_full_analysis_set_plus_decoy_hla/GRCh38_full_analysis_set_plus_decoy_hla.fa \
    -o $data_dir/chaisson/HG00514/manta/results/

# HG00733
perl ~/CNV-JACG-1.1/bin/CNV-JACG.pl \
    -p $data_dir/chaisson/smoove/HG00733.smoove.bed \
    -b $data_dir/chaisson/BAM/HG00733.alt_bwamem_GRCh38DH.20150715.PUR.high_coverage.cram \
    -r $data_dir/ref/GRCh38_full_analysis_set_plus_decoy_hla/GRCh38_full_analysis_set_plus_decoy_hla.fa \
    -o $data_dir/chaisson/HG00733/smoove/results/
perl ~/CNV-JACG-1.1/bin/CNV-JACG.pl \
    -p $data_dir/chaisson/manta/HG00514.manta.bed \
    -b $data_dir/chaisson/BAM/HG00733.alt_bwamem_GRCh38DH.20150715.PUR.high_coverage.cram \
    -r $data_dir/ref/GRCh38_full_analysis_set_plus_decoy_hla/GRCh38_full_analysis_set_plus_decoy_hla.fa \
    -o $data_dir/chaisson/HG00733/manta/results/

# NA19240
perl ~/CNV-JACG-1.1/bin/CNV-JACG.pl \
    -p $data_dir/chaisson/smoove/NA19240.smoove.bed \
    -b $data_dir/chaisson/BAM/NA19240.alt_bwamem_GRCh38DH.20150715.YRI.high_coverage.cram \
    -r $data_dir/ref/GRCh38_full_analysis_set_plus_decoy_hla/GRCh38_full_analysis_set_plus_decoy_hla.fa \
    -o $data_dir/chaisson/NA19240/smoove/results/
perl ~/CNV-JACG-1.1/bin/CNV-JACG.pl \
    -p $data_dir/chaisson/manta/NA19240.manta.bed \
    -b $data_dir/chaisson/BAM/NA19240.alt_bwamem_GRCh38DH.20150715.YRI.high_coverage.cram \
    -r $data_dir/ref/GRCh38_full_analysis_set_plus_decoy_hla/GRCh38_full_analysis_set_plus_decoy_hla.fa \
    -o $data_dir/chaisson/NA19240/manta/results/
