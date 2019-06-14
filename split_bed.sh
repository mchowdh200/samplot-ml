#!/bin/bash

BED=$1
RATIO=0.8
LINES=$(wc -l $BED | cut -f1 -d' ')
TRAIN_LINES=$(python -c "from math import ceil; print(ceil($RATIO*$LINES))")
let TEST_LINES=$LINES-$TRAIN_LINES

shuf $BED > data/shuffled.bed

head data/shuffled.bed --lines=$TRAIN_LINES > data/train.bed
tail data/shuffled.bed --lines=$TEST_LINES > data/test.bed

