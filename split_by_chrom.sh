#!/bin/bash

# path to list of image paths
DATA_LIST=$1
OUT_DIR=$(dirname $1)

TRAIN=$(seq 19 | awk '{print "chr"$1"_"}')
TEST=$(printf "20\n21\n22\nX\nY" | awk '{print "chr"$1"_"}')

grep -f <(printf "$TRAIN") $1 > $OUT_DIR/train.txt
grep -f <(printf "$TEST") $1 > $OUT_DIR/val.txt


