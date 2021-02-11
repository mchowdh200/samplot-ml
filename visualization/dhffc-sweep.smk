
config = "dhffc-sweep-config.yaml"
input_vcf = config['input_vcf']
truth_set = config['truth_set']
outdir = config['outdir']
logdir = config['logdir']
dhffc_range = range(0, 1.1, 0.01)

rule all:
    input:
        # "/dhffc-sweep.png"
        

rule filter_dhffc:
    input:
        config["input_vcf"]
    output:
        vcf = temp("{config[outdir]}/vcf/filtered-lt-{dhffc}.vcf.gz"),
        index = temp("{config[outdir]}/vcf/filtered-lt-{dhffc}.vcf.gz.tbi")
    shell:
        """bcftools view -i 'DHFFC < {wildcards.dhffc}' |
           bgzip -c > {output.vcf}
           tabix {output.vcf}"""

rule evaluate:
    input:
        filtered = "{config[outdir]}/vcf/filtered-lt-{dhffc}.vcf.gz",
        truth_set = config["truth_set"]
    output:
        text="{config[outdir]}/vcf/truvari-{dhffc}.txt",
        dir=temp(directory("{config[outdir]}/vcf/truvari-{dhffc}"))
    shell:
        """bash truvari.sh -b {input.truth_set} \
                           -c {input.filtered} \
                           -o {output.dir} &> {output.txt}"""
        
        
        
