#!/bin/bash
set -e

# arg parser -----------------------------------------------------------------
while (( "$#" )); do
    case "$1" in
        -c|--chrom)
            CHROM=$2
            shift 2;;
        -s|--start)
            START=$2
            shift 2;;
        -e|--end)
            END=$2
            shift 2;;
        -n|--sample)
            SAMPLE=$2
            shift 2;;
        -g|--genotype)
            GENOTYPE=$2
            shift 2;;
        -m|--min-mqual)
            MIN_MQ=$2
            shift 2;;
        -f|--fasta)
            FASTA=$2
            shift 2;;
        -b|--bam-file)
            BAM_FILE=$2
            shift 2;;
        -d|--bam-dir)
            BAM_DIR=$2
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
[[ -z $CHROM ]] && echo Missing argument --chrom && exit 1
[[ -z $START ]] && echo Missing argument --start && exit 1
[[ -z $END ]] && echo Missing argument --end && exit 1
[[ -z $SAMPLE ]] && echo Missing argument --sample && exit 1
[[ -z $GENOTYPE ]] && echo Missing argument --genotype && exit 1
# [[ -z $FASTA ]] && echo Missing argument --fasta && exit 1
# [[ -z $BAM_LIST ]] && echo Missing argument --bam-list && exit 1
# [[ -z $BAM_DIR ]] && echo Missing argument --bam-dir && exit 1
[[ -z $OUT_DIR ]] && echo Missing argument --out-dir && exit 1
[[ -z $MIN_MQ ]] && MIN_MQ=0

# generate image -------------------------------------------------------------
START_DIR=$PWD
cd $BAM_DIR


if [[ -z $BAM_FILE ]]; then # no list provided; use BAM_DIR contents
    # BAMS=$BAM_DIR/$(ls BAM_DIR)
    BAMS=$(find $BAM_DIR -name '*.cram' -or -name '*.bam')
else
    BAMS=$BAM_FILE
fi
# output file
OUT=$OUT_DIR/${CHROM}_${START}_${END}_${SAMPLE}_${GENOTYPE}.png
echo $OUT


if [ ! -f $OUT ]; then
    SVLEN=$(bc <<< "scale=0; ($END-$START)/1")
    # WINDOW=$(bc <<< "scale=0; (1.5*$SVLEN)/1")

    if [[ ! -z $FASTA ]]; then
        FASTA_FLAG="-r $FASTA" # if we didn't provide fasta then the flag var will be unset
    fi
    
    if [[ $SVLEN -gt 5000 ]]; then
        # Too large to plot. Just plot flanking regions and 500 bases around breakpoints
        samplot.py \
            -c $CHROM -s $START -e $END \
            --min_mqual $MIN_MQ -t DEL \
            -b $BAMS -o $OUT $FASTA_FLAG --zoom 1000 --no_haplotypes
    else
        samplot.py \
            -c $CHROM -s $START -e $END \
            --min_mqual $MIN_MQ -t DEL \
            -b $BAMS -o $OUT $FASTA_FLAG --no_haplotypes
            # -b $BAMS -o $OUT $FASTA_FLAG -w $WINDOW --no_haplotypes
    fi
fi
cd $START_DIR
