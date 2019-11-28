#!/bin/bash

data_dir=$1 # parent directory where cropped image dir can be found
data_list=$(ls -U $data_dir/crop)

# split into train/validation sets (by chromosome)
train_chrms=$(echo ^{{4..22},X,Y}_ | tr " " "|")
val_chrms=$(echo ^{1..3}_ | tr " " "|")

echo "$data_list" | grep -E "$train_chrms" | shuf > $data_dir/train.txt
echo "$data_list" | grep -E "$val_chrms" | shuf > $data_dir/val.txt

# ratio=0.90
# lines=$(echo "$data_list" | wc -l )
# train_lines=$(python -c "from math import ceil; print(int(ceil($ratio*$lines)))")
# let test_lines=$lines-$train_lines
# data_list=$(echo "$data_list" | shuf)
# echo "$data_list" | head --lines=$train_lines > $data_dir/train.txt
# echo "$data_list" | tail --lines=$test_lines > $data_dir/val.txt
