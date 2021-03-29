## General TODOs
# * provide example scripts on how to programmaticaly create sample: file pairs
#   in the config yaml.
# * or make a rule that looks at the RG tags to get sample names from list of bams

import os
import functools
from glob import glob
from config_utils import Conf


## Setup
################################################################################
configfile: 'conf/samplot-ml-predict.yaml'
conf = Conf(config)

## Rules
################################################################################
rule All:
    input:
        expand(f'{conf.outdir}/{{sample}}-predictions.bed', sample=conf.samples)


gargs = f'{conf.outdir}/bin/gargs'
rule InstallGargs:
    """
    Install system appropriate binary of gargs for use in image generation rules.
    """
    output:
        gargs
    shell:
        f'bash scripts/install_gargs.sh {gargs}'


rule GetReference:
    output:
        fasta = conf.fasta.output,
        fai = conf.fai.output
    run:
        shell(conf.fasta.get_cmd())
        shell(conf.fai.get_cmd())


rule GetBaseVCF:
    output:
        conf.vcf.output
    run:
        shell(conf.vcf.get_cmd())


rule GetDelRegions:
    """
    Get a sample's del regions from the vcf in bed format
    """
    input:
        conf.vcf.output
    output:
        f'{conf.outdir}/bed/{{sample}}-del-regions.bed'
    conda:
        'envs/samplot.yaml'
    shell:
        f"""
        [[ ! -d {conf.outdir}/bed ]] && mkdir {conf.outdir}/bed
        bash scripts/get_del_regions.sh {{input}} {{wildcards.sample}} > {{output}}
        """


### TODO seems a bit messy
def get_images(rule, wildcards):
    """
    Return list of output images from the GenerateImages/CropImages checkpoints
    """
    # gets output of the checkpoint (a directory) and re-evals workflow DAG 
    if rule == "GenerateImages":
        image_dir = checkpoints.GenerateImages.get(sample=wildcards.sample).output[0]
    elif rule == "CropImages":
        image_dir = checkpoints.CropImages.get(sample=wildcards.sample).output[0]
    else:
        raise ValueError(f'Unknown argument for rule: {rule}.'
                             'Must be "GenerateImages" or "CropImages"')
    return glob(f'{image_dir}/*.png')


checkpoint GenerateImages:
    """
    Images from del regions for a given sample.
    """
    threads: workflow.cores
    input:
        # bam/ file: from config. could be a url
        # therefore will not be tracked by snakemake
        gargs_bin = gargs,
        fasta = conf.fasta.output,
        fai = conf.fai.output,
        regions = rules.GetDelRegions.output
    output:
        directory(f'{conf.outdir}/img/{{sample}}')
    params:
        bam = lambda wildcards: conf.alignments[wildcards.sample]
    conda:
        'envs/samplot.yaml'
    shell:
        # TODO put the gen_img.sh script into a function in images_from_regions.sh
        f"""
        [[ ! -d {conf.outdir}/img ]] && mkdir {conf.outdir}/img
        bash scripts/images_from_regions.sh \\
            --gargs-bin {{input.gargs_bin}} \\
            --fasta {{input.fasta}} \\
            --regions {{input.regions}} \\
            --bam {{params.bam}} \\
            --outdir {conf.outdir}/img/{{wildcards.sample}} \\
            --delimiter {conf.delimiter} \\
            --processes {{threads}}
        """


checkpoint CropImages:
    """
    Crop axes and text from images to prepare for samplot-ml input
    """
    threads: workflow.cores
    input:
        gargs_bin = gargs,
        imgs = functools.partial(get_images, 'GenerateImages')
    output:
        # [f'{conf.outdir}/crop/{{sample}}/{os.path.basename(image)}'
        #  for image in get_images()]
        directory(f'{conf.outdir}/crop/{{sample}}')
    conda:
        'envs/samplot.yaml'
    shell:
        f"""
        [[ ! -d {conf.outdir}/crop ]] && mkdir {conf.outdir}/crop
        bash scripts/crop.sh -i {conf.outdir}/img/{{wildcards.sample}} \\
                             -o {{output}} \\
                             -p {{threads}} \\
                             -g {{input.gargs_bin}}
        """

rule CreateImageList:
    """
    Samplot-ml needs list of input images. This rule takes the list
    of a sample's cropped images and puts them in a text file.
    """
    input:
        functools.partial(get_images, 'CropImages')
    output:
        temp(f'{conf.outdir}/{{sample}}-cropped-imgs.txt')
    run:
        with open(output[0], 'w') as out:
            for image_file in input:
                out.write(f'{image_file}\n')


rule PredictImages:
    """
    TODO fixed model path for now
    TODO figure out how to do multithreading cpu with tensorflow (might be automatic?)
    TODO figure out memory usage of cpu based model and how it scales with batch size
    TODO put batch size into config?
    TODO setup config to handle gpu/cpu in the conda and resources directives
    Feed images into samplot-ml to get a bed file of predictions.
    Prediction format (tab separated):
        - chrm start end p_ref p_het p_alt
    """
    threads: workflow.cores
    input:
        f'{conf.outdir}/{{sample}}-cropped-imgs.txt'
    output:
        f'{conf.outdir}/{{sample}}-predictions.bed'
    conda:
        'envs/tensorflow.yaml'
    shell:
        """
        python scripts/predict.py \\
            --image-list {input} \\
            --delimiter {conf.delimiter} \\
            --processes {threads} \\
            --batch-size {threads} \\
            --model-path saved_models/samplot-ml.h5 \\
        > {output}
        """


## TODO
################################################################################
# use annotate.py to output the ml vcf
# change annotate.py to not filter out 0/0 genotypes
# try adding format fields for old gt and p_ref, p_het, p_alt
rule AnnotateVCF:




