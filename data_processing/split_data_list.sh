#!/bin/bash


DATA=$1 # eg data/high_cov or data/low_cov
ls $DATA/crop > $DATA/data_list.txt
FILE=$DATA/data_list.txt

RATIO=0.95
LINES=$(wc -l $FILE | cut -f1 -d' ')
TRAIN_LINES=$(python -c "from math import ceil; print(ceil($RATIO*$LINES))")
let TEST_LINES=$LINES-$TRAIN_LINES

shuf $FILE > $DATA/shuffled.txt

head $DATA/shuffled.txt --lines=$TRAIN_LINES > $DATA/train.txt
tail $DATA/shuffled.txt --lines=$TEST_LINES > $DATA/val.txt

