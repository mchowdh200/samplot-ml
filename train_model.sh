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

# TODO create folder structure?
    # or create the folder structure in the setup script
    # and just set the $data_dir here
# TODO download the train/val.txt
# TODO generate the tfrec lists

# ec2
# set up directory structure ----------------------------------------------------------------------
data_dir=~/data
mkdir $data_dir
mkdir $data_dir
mkdir $data_dir/saved_models

# download data listings --------------------------------------------------------------------------
s3_source="s3://layerlab/samplot-ml/1kg/training_sets/1kg_duphold.12.12.18"
aws s3 cp $s3_source/train.txt $data_dir
aws s3 cp $s3_source/val.txt $data_dir

aws s3 ls s3://layerlab/samplot-ml/1kg/training_sets/1kg_duphold.12.12.18/train/ \
    | sed 's/ \+/\t/g' \
    | cut -f4 \
    | awk '{print "s3://layerlab/samplot-ml/1kg/training_sets/1kg_duphold.12.12.18/train/"$0}' > $data_dir/train_tfrec_list.txt
aws s3 ls s3://layerlab/samplot-ml/1kg/training_sets/1kg_duphold.12.12.18/val/ \
    | sed 's/ \+/\t/g' \
    | cut -f4 \
    | awk '{print "s3://layerlab/samplot-ml/1kg/training_sets/1kg_duphold.12.12.18/val/"$0}' > $data_dir/val_tfrec_list.txt


 python3 run.py train \
    --verbose 1 \
    --processes 4 \
    --batch-size 32 \
    --epochs 50 \
    --model-type CNN \
    --num-classes 3 \
    --train-list $data_dir/train.txt \
    --val-list $data_dir/val.txt \
    --train-tfrec-list $data_dir/train_tfrec_list.txt \
    --val-tfrec-list $data_dir/val_tfrec_list.txt \
    --learning-rate 0.05 \
    --momentum 0.9 \
    --label-smoothing 0.05 \
    --save-to $data_dir/saved_models/CNN.$(date +%m.%d.%H).h5
