#!/bin/bash

set -eu

# path to list of image paths
DATA_LIST=$1
OUT_DIR=$(dirname $1)
echo $OUT_DIR

# TRAIN=$(seq 19 | awk '{print "chr"$1"_"}')
# TEST=$(printf "20\n21\n22\nX\nY" | awk '{print "chr"$1"_"}')


TRAIN=$(printf "$(seq 4 22)\nX\nY" | awk '{print "chr"$1"_"}')
TEST=$(seq 1 3 | awk '{print "chr"$1"_"}')
echo "$TRAIN"
echo "$TEST"

grep -f <(printf "$TRAIN") $1  > $OUT_DIR/train.txt
grep -f <(printf "$TEST") $1  > $OUT_DIR/val.txt


