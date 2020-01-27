#!/bin/bash
#SBATCH -p tesla-k40
#SBATCH --job-name=samplot-ml
#SBATCH --mail-type=ALL
#SBATCH --mail-user=much8161@colorado.edu
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --mem=128gb
#SBATCH --gres=gpu:1
#SBATCH --time=24:00:00
#SBATCH --output=/Users/much8161/Repositories/samplot-ml/model_code/training-log.out
#SBATCH --error=/Users/much8161/Repositories/samplot-ml/model_code/training-log.err

# fiji
# data_dir=/scratch/Shares/layer/projects/samplot/ml/data/1kg/high_cov/
# processes=16

# home server
# data_dir="/data/1kg_high_cov/"
# processes=4


# ec2
# set up directory structure ----------------------------------------------------------------------
# data_dir=/mnt/local/data
data_dir="/data/1kg_high_cov/1kg_exclude/"
# mkdir $data_dir
# mkdir $data_dir/train
# mkdir $data_dir/val
mkdir $data_dir/saved_models

# download data listings --------------------------------------------------------------------------
# s3_source="s3://layerlab/samplot-ml/1kg/training_sets/1kg_duphold.12.12.18"
# aws s3 cp $s3_source/train.txt $data_dir
# aws s3 cp $s3_source/val.txt $data_dir

# download the tfrecords --------------------------------------------------------------------------
# aws s3 cp --recursive \
#     s3://layerlab/samplot-ml/1kg/training_sets/1kg_duphold.12.12.18/train/ $data_dir/train
# aws s3 cp --recursive \
#     s3://layerlab/samplot-ml/1kg/training_sets/1kg_duphold.12.12.18/val/ $data_dir/val

find $data_dir/train/*.tfrec > $data_dir/train_tfrec_list.txt
find $data_dir/val/*.tfrec > $data_dir/val_tfrec_list.txt


model_name=CNN.$(date +%m.%d.%H).h5

 python3 run.py train \
    --verbose 1 \
    --processes 4 \
    --batch-size 32 \
    --epochs 200 \
    --model-type CNN \
    --num-classes 3 \
    --train-list $data_dir/train.txt \
    --val-list $data_dir/val.txt \
    --train-tfrec-list $data_dir/train_tfrec_list.txt \
    --val-tfrec-list $data_dir/val_tfrec_list.txt \
    --learning-rate 0.2 \
    --momentum 0.0 \
    --weight-decay 5e-5 \
    --label-smoothing 0.05 \
    --save-to $data_dir/saved_models/$model_name

 # upload model to s3
 aws s3 cp $data_dir/saved_models/$model_name s3://layerlab/samplot-ml/saved_models/

