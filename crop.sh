#!/bin/bash
#SBATCH -p short
#SBATCH --job-name=cropping
#SBATCH --mail-type=ALL
#SBATCH --mail-user=much8161@colorado.edu
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --mem=32gb
#SBATCH --time=24:00:00               # Time limit hrs:min:sec
#SBATCH --output=/Users/much8161/Repositories/samplot-ml/out/crop.out
#SBATCH --error=/Users/much8161/Repositories/samplot-ml/out/crop.err


# script to take samplot image and crop out the
# axes, title, etc.
function crop
{
    # relative path of input image
    p=imgs/$1

    b=`basename $p .png`
    o=crop/${b}.png

    convert \
        -crop 2090x575+175+200 \
        -fill white \
        -draw "rectangle 0,30 500,50" \
        -draw "rectangle 600,0 700,50" \
        $p $o
    echo $o

}

export -f crop

SCRIPT_DIR=$PWD
# DATA_DIR=/scratch/Shares/layer/projects/samplot/ml/data/giab
DATA_DIR=/scratch/Shares/layer/projects/samplot/ml/data/1kg/high_cov

cd $DATA_DIR

ls $DATA_DIR/imgs | gargs -p 64 'crop {}'

cd $SCRIPT_DIR
