#!/bin/bash
set -e

# arg parser -----------------------------------------------------------------
while (( "$#" )); do
    case "$1" in
        -h|--help)
            HELP=true
            shift 1;;
        -m|--model-path)
            MODEL_PATH=$2
            shift 2;;
        -d|--data-list)
            DATA_LIST=$2
            shift 2;;
        -a|--augmentation)
            AUGMENTATION=true
            shift 1;;
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

HELP_MESSAGE="
-h/--help:              Display help message
-m/--model-path [PATH]: Path of tf.keras model (.h5)
-d/--data-list [LIST]:  List of file paths of input images
-a/--augmentation:      Whether or not to perform test time augmentation of images
-v/--vcf [VCF]:         VCF to annotate with model predictions
-o/--out-dir [OUT]:     Output directory of annotate VCF
"
[[ ! -z $HELP ]] && echo "$HELP_MESSAGE" && exit
[[ -z $MODEL_PATH ]] && printf "\nMissing argument -m/--model-path\n" && echo "$HELP_MESSAGE" && exit 1
[[ -z $DATA_LIST ]] && printf "\nMissing argument -d/--data-list\n" && echo "$HELP_MESSAGE" && exit 1
[[ -z $VCF ]] && echo "\nMissing argument -v/--vcf\n" && echo "$HELP_MESSAGE" && exit 1
[[ -z $OUT_DIR ]] && printf "\nMissing argument -o/--out-dir\n" && echo "$HELP_MESSAGE" && exit 1


[[ ! -z $AUGMENTATION ]] && AUG_FLAG="-a"
BED=$OUT_DIR/$(basename $MODEL_PATH .h5).bed
echo "storing prediction in $BED"

if [[ ! -f $BED ]]; then
    if [ -d $OUT_DIR ]; then
        rm -r $OUT_DIR
    fi

    mkdir $OUT_DIR
        python3 ../run.py predict -mp $MODEL_PATH -h5 -i $DATA_LIST $AUG_FLAG \
            | python3 pred2bed.py $BED
fi

python3 annotate.py $VCF $BED | bgzip -c > $OUT_DIR/ml.vcf.gz
tabix $OUT_DIR/ml.vcf.gz
