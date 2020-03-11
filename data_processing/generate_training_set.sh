set -u

# s3_destination="s3://layerlab/samplot-ml/1kg/duphold_callset"
# s3_destination="s3://layerlab/samplot-ml/1kg/training_sets/1kg_exclude_sensitive"
n_tasks=$1
s3_destination=$2

data_dir=/mnt/local/data

mkdir $data_dir
mkdir $data_dir/bed
mkdir $data_dir/ref
mkdir $data_dir/cram
mkdir $data_dir/imgs
mkdir $data_dir/crop

# get training regions
aws s3 cp $s3_destination/ref.bed $data_dir/bed
aws s3 cp $s3_destination/het.bed $data_dir/bed
aws s3 cp $s3_destination/alt.bed $data_dir/bed
cat $data_dir/bed/ref.bed $data_dir/bed/het.bed $data_dir/bed/alt.bed > $data_dir/bed/training_regions.bed
training_regions=$data_dir/bed/training_regions.bed

# get reference genome
aws s3 cp s3://1000genomes/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa $data_dir/ref
aws s3 cp s3://1000genomes/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai $data_dir/ref
fasta=$data_dir/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa

# get CRAM indices
cram_list=/home/ubuntu/samplot-ml/data_listings/1kg_high_cov_crams_s3.txt
cat $cram_list | gargs -p $n_tasks "aws s3 cp {}.crai $data_dir/cram"

# generate the images
cat $training_regions | gargs -p $n_tasks \
    "bash gen_img.sh \\
        --chrom {0} --start {1} --end {2} --sample {3} --genotype {6} \\
        --min-mqual 10 \\
        --fasta $fasta \\
        --bam-list $cram_list \\
        --bam-dir $data_dir/cram \\
        --out-dir $data_dir/imgs"

# crop the images
ls $data_dir/imgs | bash crop.sh -p $n_tasks -d $data_dir/imgs -o $data_dir/crop

# upload results to s3
aws s3 cp --recursive $data_dir/imgs $s3_destination/imgs
aws s3 cp --recursive $data_dir/crop $s3_destination/crop


