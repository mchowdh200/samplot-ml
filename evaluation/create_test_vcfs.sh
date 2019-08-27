#!/bin/bash
set -e

# arg parser -----------------------------------------------------------------
while (( "$#" )); do
    case "$1" in
        -m|--model-path)
            MODEL_PATH=$2
            shift 2;;
        -d|--data-list)
            DATA_LIST=$2
            shift 2;;
        -v|--vcf)
            VCF=$2
            shift 2;;
        -o|--out-dir)
            OUT_DIR=$2
            shift 2;;
        --) # end argument parsing
            shift
            break;;
        -*|--*=) # unsupported flags
            echo "Error: Unsupported flag $1" >&2
            exit 1;;
    esac
done
[[ -z $MODEL_PATH ]] && echo Missing argument --model-path && exit 1
[[ -z $DATA_LIST ]] && echo Missing argument --data-list && exit 1
[[ -z $VCF ]] && echo Missing argument --vcf && exit 1
[[ -z $OUT_DIR ]] && echo Missing argument --out-dir && exit 1


# MODEL_PATH=$1
# DATA_LIST=data/giab/giab_list.txt

# OUT_DIR=data/giab/VCF/$(basename $MODEL_PATH .h5)
# VCF=data/giab/VCF/HG002-smoove.genotyped.vcf.gz
if [ -d $OUT_DIR ]; then
    rm -r $OUT_DIR
fi
mkdir $OUT_DIR
mkdir $OUT_DIR/VCF
mkdir $OUT_DIR/predictions

BED=$OUT_DIR/predictions/$(basename $MODEL_PATH .h5).bed

python3 run.py predict -mp $MODEL_PATH -h5 -i $DATA_LIST \
    | python3 pred2bed.py > $OUT_DIR/predictions/$(basename $MODEL_PATH .h5).bed

python3 annotate.py $VCF $BED | bgzip -c > $OUT_DIR/VCF/ml.vcf.gz
tabix $OUT_DIR/VCF/ml.vcf.gz


# bash filter_vcf.sh $OUT_DIR
# source deactivate




