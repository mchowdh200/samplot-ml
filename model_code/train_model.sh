#!/bin/bash
#SBATCH -p tesla-k40
#SBATCH --job-name=samplot-ml
#SBATCH --mail-type=ALL
#SBATCH --mail-user=much8161@colorado.edu
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --mem=128gb
#SBATCH --gres=gpu:1
#SBATCH --time=24:00:00
#SBATCH --output=/Users/much8161/Repositories/samplot-ml/model_code/training-log.out
#SBATCH --error=/Users/much8161/Repositories/samplot-ml/model_code/training-log.err

# python3 scratch.py
./run.py train \
    --verbose 1 \
    --processes 16 \
    --batch-size 32 \
    --epochs 50 \
    --model-type CNN \
    --num-classes 3 \
    --data-dir /scratch/Shares/layer/projects/samplot/ml/data/1kg/high_cov/ \
    --learning-rate 0.05 \
    --momentum 0.9 \
    --label-smoothing 0.05 \
    --save-to ../saved_models/CNN.$(date +%m.%d.%H).h5
