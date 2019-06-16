#!/bin/bash

FILE=$1
RATIO=0.8
LINES=$(wc -l $FILE | cut -f1 -d' ')
TRAIN_LINES=$(python -c "from math import ceil; print(ceil($RATIO*$LINES))")
let TEST_LINES=$LINES-$TRAIN_LINES

shuf $FILE > data/shuffled.txt

head data/shuffled.txt --lines=$TRAIN_LINES > data/train.txt
tail data/shuffled.txt --lines=$TEST_LINES > data/test.txt

