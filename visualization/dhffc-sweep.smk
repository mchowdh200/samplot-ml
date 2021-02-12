import numpy as np

configfile: "dhffc-sweep-config.yaml"
input_vcf = config['input_vcf']
truth_set = config['truth_set']
outdir = config['outdir']
logdir = config['logdir']

samples = ["HG002", "HG00514", "HG00733", "NA19240"]
sample = lambda w: w.samples # resolves sample name wildcard
dhffc_range = np.around(np.linspace(0, 1.0, 101), 2)

rule all:
    input:
        # config["outdir"]+"/dhffc-sweep.png"
        expand(
            config["outdir"]+"/{sample}/truvari-{dhffc}.txt",
            sample=samples,
            dhffc=dhffc_range)

rule filter_dhffc:
    input:
        config['outdir'][sample]
    output:
        vcf = temp(config["outdir"]+"/{sample}/filtered-lt-{dhffc}.vcf.gz"),
        index = temp(config["outdir"]+"/{sample}/filtered-lt-{dhffc}.vcf.gz.tbi")
    shell:
        """bcftools view -i 'DHFFC <= {wildcards.dhffc}' {input} |
           bgzip -c > {output.vcf}
           tabix {output.vcf}"""

rule evaluate:
    input:
        filtered = config["outdir"]+"/{sample}/filtered-lt-{dhffc}.vcf.gz",
        filtered_index = config["outdir"]+"/{sample}/filtered-lt-{dhffc}.vcf.gz.tbi",
        truth_set = config["truth_set"]["{sample}"]
    output:
        txt=config["outdir"]+"/{sample}/truvari-{dhffc}.txt",
        dir=temp(directory(config["outdir"]+"/{sample}/truvari-{dhffc}"))
    log:
        config["logdir"]+"/{sample}/evaluate-{dhffc}.log"
    shell:
        """bash ../evaluation/truvari.sh -b {input.truth_set} \
                           -c {input.filtered} \
                           -o {output.dir} 2>&1 |
           tee {log} > {output.txt}"""
        
        
        
