#!/bin/bash

set -eu

# # === get references ===========================================================
# [[ ! -d /mnt/local/data/ref ]] && mkdir /mnt/local/data/ref
# aws s3 cp --recursive \
#     s3://layerlabcu/ref/genomes/GRCh38_full_analysis_set_plus_decoy_hla/ \
#     /mnt/local/data/GRCh38_full_analysis_set_plus_decoy_hla/
# aws s3 cp --recursive \
#     s3://layerlabcu/ref/genomes/g1k_v37_decoy/ \
#     /mnt/local/data/g1k_v37_decoy/
# aws s3 cp --recursive \
#     s3://layerlabcu/ref/genomes/hs37d5/ \
#     /mnt/local/data/hs37d5/


# # === HG002 ====================================================================
# [[ ! -d /mnt/local/data/HG002 ]] && mkdir /mnt/local/data/HG002
# # === bam ===
# [[ ! -d /mnt/local/data/HG002/BAM ]] && mkdir /mnt/local/data/HG002/BAM
# aws s3 cp s3://layerlabcu/samplot-ml/HG002/BAM/hg002.bam \
#           /mnt/local/data/HG002/BAM/
# aws s3 cp s3://layerlabcu/samplot-ml/HG002/BAM/hg002.bam.bai \
#           /mnt/local/data/HG002/BAM/
# # === vcf ===
# [[ ! -d /mnt/local/data/HG002/VCF ]] && mkdir /mnt/local/data/HG002/VCF
# aws s3 cp s3://layerlabcu/samplot-ml/HG002/smoove/VCF/HG002-smoove.genotyped.del.vcf.gz \
#           /mnt/local/data/HG002/VCF/HG002-smoove.vcf.gz
# aws s3 cp s3://layerlabcu/samplot-ml/HG002/manta/VCF/diploidSV.del.duphold.vcf.gz \
#           /mnt/local/data/HG002/VCF/HG002-manta.vcf.gz
# === paragraph ===
# read_depth=$(~/paragraph/build/bin/idxdepth \
#     -b /mnt/local/data/HG002/BAM/hg002.bam \
#     -r /mnt/local/data/g1k_v37_decoy/g1k_v37_decoy.fa |
#     grep -m 1 'depth' | sed -E 's/\s+//g' | cut -d':' -f2)
# read_length=$(samtools view /mnt/local/data/HG002/BAM/hg002.bam |
#     awk '{print length($10)}' | head -100000 |
#     python3 -c 'import sys, statistics; print(statistics.mean([float(i.rstrip()) for i in sys.stdin]))')
# printf "id\tpath\tdepth\tread length\nHG002\t/mnt/local/data/HG002/BAM/hg002.bam\t$read_depth\t$read_length\n" > /mnt/local/data/HG002/samples.txt

# smoove
# [[ ! -d /mnt/local/data/HG002/smoove ]] && mkdir /mnt/local/data/HG002/smoove
# python3 ~/paragraph/build/bin/multigrmpy.py \
#     -i /mnt/local/data/HG002/VCF/HG002-smoove.vcf.gz \
#     -m /mnt/local/data/HG002/samples.txt \
#     -r /mnt/local/data/g1k_v37_decoy/g1k_v37_decoy.fa \
#     -o /mnt/local/data/HG002/smoove/

# manta
# [[ ! -d /mnt/local/data/HG002/manta ]] && mkdir /mnt/local/data/HG002/manta
# python3 ~/paragraph/build/bin/multigrmpy.py \
#     -i /mnt/local/data/HG002/VCF/HG002-manta.vcf.gz \
#     -m /mnt/local/data/HG002/samples.txt \
#     -r /mnt/local/data/g1k_v37_decoy/g1k_v37_decoy.fa \
#     -o /mnt/local/data/HG002/manta/

# === HG002.2x250 ==============================================================
[[ ! -d /mnt/local/data/HG002.2x250 ]] && mkdir /mnt/local/data/HG002.2x250
# === bam ===
[[ ! -d /mnt/local/data/HG002.2x250/BAM ]] && mkdir /mnt/local/data/HG002.2x250/BAM
# aws s3 cp s3://layerlabcu/samplot-ml/HG002/BAM/HG002.2x250.sorted.bam \
#           /mnt/local/data/HG002.2x250/BAM/
# aws s3 cp s3://layerlabcu/samplot-ml/HG002/BAM/HG002.2x250.sorted.bam.bai \
#           /mnt/local/data/HG002.2x250/BAM/
# === vcf ===
[[ ! -d /mnt/local/data/HG002.2x250/VCF ]] && mkdir /mnt/local/data/HG002.2x250/VCF
aws s3 cp s3://layerlabcu/samplot-ml/HG002/smoove.2x250/HG002-smoove.genotyped.del.vcf.gz \
          /mnt/local/data/HG002.2x250/VCF/HG002.2x250-smoove.vcf.gz
aws s3 cp s3://layerlabcu/samplot-ml/HG002/manta.2x250/diploidSV.duphold.del.vcf.gz \
          /mnt/local/data/HG002.2x250/VCF/HG002.2x250-manta.vcf.gz
# === paragraph ===
read_depth=$(~/paragraph/build/bin/idxdepth \
    -b /mnt/local/data/HG002.2x250/BAM/HG002.2x250.sorted.bam \
    -r /mnt/local/data/hs37d5/hs37d5.fa |
    grep -m 1 'depth' | sed -E 's/\s+//g' | cut -d':' -f2)
read_length=$(samtools view /mnt/local/data/HG002.2x250/BAM/HG002.2x250.sorted.bam |
    awk '{print length($10)}' | head -100000 |
    python3 -c 'import sys, statistics; print(statistics.mean([float(i.rstrip()) for i in sys.stdin]))')
printf "id\tpath\tdepth\tread length\nHG002\t/mnt/local/data/HG002.2x250/BAM/HG002.2x250.sorted.bam\t$read_depth\t$read_length\n" > /mnt/local/data/HG002.2x250/samples.txt

# smoove
[[ ! -d /mnt/local/data/HG002.2x250/smoove ]] && mkdir /mnt/local/data/HG002.2x250/smoove
python3 ~/paragraph/build/bin/multigrmpy.py \
    -i /mnt/local/data/HG002.2x250/VCF/HG002.2x250-smoove.vcf.gz \
    -m /mnt/local/data/HG002.2x250/samples.txt \
    -r /mnt/local/data/hs37d5/hs37d5.fa \
    -o /mnt/local/data/HG002.2x250/smoove/

# manta
[[ ! -d /mnt/local/data/HG002.2x250/manta ]] && mkdir /mnt/local/data/HG002.2x250/manta
python3 ~/paragraph/build/bin/multigrmpy.py \
    -i /mnt/local/data/HG002.2x250/VCF/HG002.2x250-manta.vcf.gz \
    -m /mnt/local/data/HG002.2x250/samples.txt \
    -r /mnt/local/data/hs37d5/hs37d5.fa \
    -o /mnt/local/data/HG002.2x250/manta/

exit

# === HG00514 =================================================================
[[ ! -d /mnt/local/data/HG00514 ]] && mkdir /mnt/local/data/HG00514
# === bam ===
[[ ! -d /mnt/local/data/HG00514/BAM ]] && mkdir /mnt/local/data/HG00514/BAM
aws s3 cp s3://layerlabcu/samplot-ml/chaisson/BAM/HG00514.alt_bwamem_GRCh38DH.20150715.CHS.high_coverage.cram \
          /mnt/local/data/HG00514/BAM/
aws s3 cp s3://layerlabcu/samplot-ml/chaisson/BAM/HG00514.alt_bwamem_GRCh38DH.20150715.CHS.high_coverage.cram.crai \
          /mnt/local/data/HG00514/BAM/
# === vcf ===
[[ ! -d /mnt/local/data/HG00514/VCF ]] && mkdir /mnt/local/data/HG00514/VCF
aws s3 cp s3://layerlabcu/samplot-ml/chaisson/VCF/smoove/HG00514-smoove.genotyped.del.vcf.gz \
          /mnt/local/data/HG00514/VCF/HG00514-smoove.vcf.gz
aws s3 cp s3://layerlabcu/samplot-ml/chaisson/VCF/manta/HG00514.manta.vcf.gz \
          /mnt/local/data/HG00514/VCF/HG00514-manta.vcf.gz

# === paragraph ===
read_depth=$(~/paragraph/build/bin/idxdepth \
    -b /mnt/local/data/HG00514/BAM/HG00514.alt_bwamem_GRCh38DH.20150715.CHS.high_coverage.cram
    -r /mnt/local/data/GRCh38_full_analysis_set_plus_decoy_hla/GRCh38_full_analysis_set_plus_decoy_hla.fa |
    grep -m 1 'depth' | sed -E 's/\s+//g' | cut -d':' -f2)
read_length=$(samtools view /mnt/local/data/HG00514/BAM/HG00514.alt_bwamem_GRCh38DH.20150715.CHS.high_coverage.cram |
    awk '{print length($10)}' | head -100000 |
    python3 -c 'import sys, statistics; print(statistics.mean([float(i.rstrip()) for i in sys.stdin]))')
printf "id\tpath\tdepth\tread length\nHG00514\t/mnt/local/data/HG00514/BAM/HG00514.alt_bwamem_GRCh38DH.20150715.CHS.high_coverage.cram\t$read_depth\t$read_length\n" > /mnt/local/data/HG00514/samples.txt

# smoove
[[ ! -d /mnt/local/data/HG00514/smoove ]] && mkdir /mnt/local/data/HG00514/smoove
python3 ~/paragraph/build/bin/multigrmpy.py \
    -i /mnt/local/data/HG00514/VCF/HG00514-smoove.vcf.gz \
    -m /mnt/local/data/HG00514/samples.txt \
    -r /mnt/local/data/GRCh38_full_analysis_set_plus_decoy_hla/GRCh38_full_analysis_set_plus_decoy_hla.fa \
    -o /mnt/local/data/HG00514/smoove/

# manta
[[ ! -d /mnt/local/data/HG00514/manta ]] && mkdir /mnt/local/data/HG00514/manta
python3 ~/paragraph/build/bin/multigrmpy.py \
    -i /mnt/local/data/HG00514/VCF/HG00514-manta.vcf.gz \
    -m /mnt/local/data/HG00514/samples.txt \
    -r /mnt/local/data/GRCh38_full_analysis_set_plus_decoy_hla/GRCh38_full_analysis_set_plus_decoy_hla.fa \
    -o /mnt/local/data/HG00514/manta/
