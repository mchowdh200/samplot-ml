#!/bin/bash
#SBATCH -p long
#SBATCH --job-name=generate_training_set
#SBATCH --mail-type=ALL
#SBATCH --mail-user=much8161@colorado.edu
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --mem=128gb
#SBATCH --time=72:00:00
#SBATCH --output=/Users/much8161/Repositories/samplot-ml/data_processing/gen_training.out
#SBATCH --error=/Users/much8161/Repositories/samplot-ml/data_processing/gen_training.err

training_regions=/scratch/Shares/layer/projects/samplot/ml/data/1kg/high_cov/BED/training_regions.bed
fasta=/scratch/Shares/layer/ref/GRCh38_full_analysis_set_plus_decoy_hla/GRCh38_full_analysis_set_plus_decoy_hla.fa
cram_list=~/Repositories/samplot-ml/data_listings/1kg_ftp_cram_list.txt
cram_dir=/scratch/Shares/layer/projects/samplot/ml/data/1kg/high_cov/alignments/
out_dir=/scratch/Shares/layer/projects/samplot/ml/data/1kg/high_cov/imgs/

cat $training_regions | gargs -p 32 \
    "bash gen_img.sh \\
        --chrom chr{0} --start {1} --end {2} --sample {3} --genotype {4} \\
        --min-mqual 5 \\
        --fasta $fasta \\
        --bam-list $cram_list \\
        --bam-dir $cram_dir \\
        --out-dir $out_dir"
