#!/bin/bash
#SBATCH -p tesla-k40
#SBATCH --gres=gpu:1
#SBATCH --job-name=evaluate
#SBATCH --mail-type=ALL
#SBATCH --mail-user=much8161@colorado.edu
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --mem=32gb
#SBATCH --time=24:00:00
#SBATCH --output=/Users/much8161/Repositories/samplot-ml/evaluation/logs/evaluate.out
#SBATCH --error=/Users/much8161/Repositories/samplot-ml/evaluation/logs/evaluate.err

set -eu

model=$1
samples=(HG002 HG00514 HG00733 NA19240)
data_dir=/scratch/Shares/layer/projects/samplot/ml/data

echo "Using $model"

echo "SMOOVE TEST SETS"
for sample in ${samples[@]}; do
    echo "SAMPLE $sample"
    if [[ ! -f $data_dir/$sample/VCF/$(basename $model .h5)/ml.vcf.gz ]]; then
        find $data_dir/$sample/crop/*.png > $data_dir/$sample/$sample-img-list.txt
        bash create_test_vcfs.sh \
            --model-path $model \
            --data-list $data_dir/$sample/$sample-img-list.txt \
            --vcf $data_dir/$sample/VCF/smoove/$sample-smoove.genotyped.del.vcf.gz \
            --out-dir $data_dir/$sample/VCF/$(basename $model .h5)
    fi

    if [[ $sample == "HG002" ]]; then 
        truth_set=HG002_SVs_Tier1_v0.6.del.vcf.gz
    else
        truth_set=$sample.BIP-unified.del.fixed.vcf.gz
    fi

    bash truvari.sh \
        -c $data_dir/$sample/VCF/$(basename $model .h5)/ml.vcf.gz \
        -b $data_dir/$sample/VCF/$truth_set \
        -o $data_dir/$sample/VCF/$(basename $model .h5)/truvari

done


echo "MANTA TEST SETS"
for sample in ${samples[@]}; do
    echo "SAMPLE $sample"
    if [[ ! -f $data_dir/$sample/VCF/manta/$(basename $model .h5)/ml.vcf.gz ]]; then
        find $data_dir/$sample/manta/crop/*.png > $data_dir/$sample/manta/$sample-img-list.txt
        bash create_test_vcfs.sh \
            --model-path $model \
            --data-list $data_dir/$sample/manta/$sample-img-list.txt \
            --vcf $data_dir/$sample/VCF/manta/results/variants/diploidSV.del.vcf.gz \
            --out-dir $data_dir/$sample/VCF/manta/$(basename $model .h5)
    fi

    if [[ $sample == "HG002" ]]; then 
        truth_set=HG002_SVs_Tier1_v0.6.del.vcf.gz
    else
        truth_set=$sample.BIP-unified.del.fixed.vcf.gz
    fi

    bash truvari.sh \
        -c $data_dir/$sample/VCF/manta/$(basename $model .h5)/ml.vcf.gz \
        -b $data_dir/$sample/VCF/$truth_set \
        -o $data_dir/$sample/VCF/manta/$(basename $model .h5)/truvari

done

