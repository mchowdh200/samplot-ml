#!/bin/bash
#SBATCH -p short
#SBATCH --job-name=samplot-spawn
#SBATCH --mail-type=ALL               # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=much8161@colorado.edu
#SBATCH --nodes=1                    # Only use a single node
#SBATCH --ntasks=64                    # Run on a single CPU
#SBATCH --mem=128gb                   # Memory limit
#SBATCH --time=24:00:00               # Time limit hrs:min:sec
#SBATCH --output=/Users/much8161/Repositories/samplot-ml/out/spawn.out
#SBATCH --error=/Users/much8161/Repositories/samplot-ml/out/spawn.err

# LOG_DIR=/scratch/Shares/layer/projects/samplot/ml/data/1kg/high_cov/log/
# LOG_DIR=/dev/null
# cat data/del.sample.bed | gargs 'rj -l '$LOG_DIR' -n {0}_{1}_{2}_{3}_{4} -c "bash gen_img.sh {0} {1} {2} {3} {4}"'
cat data/BED/del.sample.bed | gargs -p 64 "bash gen_img.sh {0} {1} {2} {3} {4}"

