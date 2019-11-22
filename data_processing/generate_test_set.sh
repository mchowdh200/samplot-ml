#!/bin/bash
#SBATCH -p short
#SBATCH --job-name=generate_test
#SBATCH --mail-type=ALL
#SBATCH --mail-user=much8161@colorado.edu
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --mem=128gb
#SBATCH --time=24:00:00
#SBATCH --output=/Users/much8161/Repositories/samplot-ml/data_processing/logs/generate_test.out
#SBATCH --error=/Users/much8161/Repositories/samplot-ml/data_processing/logs/generate_test.err

data_dir=/scratch/Shares/layer/projects/samplot/ml/data/
samples=(HG002 HG00514 HG00733 NA19240)
ref=/scratch/Shares/layer/ref/GRCh38_full_analysis_set_plus_decoy_hla/GRCh38_full_analysis_set_plus_decoy_hla.fa

# smoove data
for sample in ${samples[@]}; do
    [[ $sample != "HG002" ]] && fasta_flag="--fasta $ref"
    bcftools query -f '%CHROM\t%POS\t%INFO/END\n' \
        $data_dir/$sample/VCF/smoove/$sample-smoove.genotyped.del.vcf.gz \
        | gargs -p 64 \
            "bash gen_img.sh \\
                --chrom {0} --start {1} --end {2} \\
                --sample $sample --genotype DEL \\
                --min-mqual 5 \\
                $fasta_flag \\
                --bam-dir $data_dir/$sample/BAM \\
                --out-dir $data_dir/$sample/imgs
                "
    bash crop.sh -p 64 -d $data_dir/$sample
done

# manta data
for sample in ${samples[@]}; do
    [[ $sample != "HG002" ]] && fasta_flag="--fasta $ref"
    bcftools query -f '%CHROM\t%POS\t%INFO/END\n' \
        $data_dir/$sample/VCF/manta/results/variants/diploidSV.del.vcf.gz \
        | gargs -p 64 \
            "bash gen_img.sh \\
                --chrom {0} --start {1} --end {2} \\
                --sample $sample --genotype DEL \\
                --min-mqual 5 \\
                $fasta_flag \\
                --bam-dir $data_dir/$sample/BAM \\
                --out-dir $data_dir/$sample/manta/imgs
                "
    bash crop.sh -p 64 -d $data_dir/$sample/manta
done

