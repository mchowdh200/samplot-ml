#!/bin/bash

set -eu

training_set_name=$1
ntasks=$2

s3_source="s3://layerlab/samplot-ml/1kg/duphold_callset/crop/"
s3_destination="s3://layerlab/samplot-ml/1kg/training_sets/${training_set_name}.$(date +%m.%d.%H)"

data_dir=/mnt/local/data
mkdir $data_dir/crop
mkdir $data_dir/train
mkdir $data_dir/val

# get the training set images
aws configure set default.s3.max_concurrent_requests $n_tasks
aws s3 cp --recursive $s3_source $data_dir/crop/


# organize images into training and validation sets
data_list=$(ls U $data_dir/crop)
val_chrms=$(echo ^chr{1..3}_ | tr " " "|")
echo "$data_list" | grep -v -E "$val_chrms" | shuf > $data_dir/train.txt
echo "$data_list" | grep -E "$val_chrms" | shuf > $data_dir/val.txt

# create the tfrecords
python3 datasets.py $data_dir

# upload to s3
aws s3 cp --recursive $data_dir/train/ $s3_destination/train/
aws s3 cp --recursive $data_dir/val/ $s3_destination/val/
