#!/bin/bash

MODEL=saved_models/CNN_chrsplit.h5

./run.py predict \
    -h5 \
    -mt CNN \
    -mp $MODEL \
    -i data/giab/giab_list.txt
