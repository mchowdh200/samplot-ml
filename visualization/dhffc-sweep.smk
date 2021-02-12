import numpy as np

configfile: "dhffc-sweep-config.yaml"
input_vcf = config['input_vcf']
truth_set = config['truth_set']
outdir = config['outdir']
logdir = config['logdir']

samples = ["HG002", "HG00514", "HG00733", "NA19240"]
dhffc_range = np.around(np.linspace(0, 1.0, 101), 2)

rule all:
    input:
        # config["outdir"]+"/dhffc-sweep.png"
        expand(config["outdir"]+"/{sample}/stats-{dhffc}.txt",
               dhffc=dhffc_range, sample=samples)

rule filter_dhffc:
    """
    Filter baseline genotyped SV vcf with DHFFC metric
    """
    input:
        lambda w: config["input_vcf"][w.sample] # resolves sample name wildcard
    output:
        vcf = temp(config["outdir"]+"/{sample}/filtered-lt-{dhffc}.vcf.gz"),
        index = temp(config["outdir"]+"/{sample}/filtered-lt-{dhffc}.vcf.gz.tbi")
    shell:
        """bcftools view -i 'DHFFC <= {wildcards.dhffc}' {input} |
           bgzip -c > {output.vcf}
           tabix {output.vcf}"""

rule evaluate:
    """
    Evaluate each DHFFC filtered vcf with truvari and output to a report file.
    """
    input:
        filtered = config["outdir"]+"/{sample}/filtered-lt-{dhffc}.vcf.gz",
        filtered_index = config["outdir"]+"/{sample}/filtered-lt-{dhffc}.vcf.gz.tbi",
        truth_set = lambda w: config["truth_set"][w.sample]
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

rule get_sample_stats:
    """
    For each sample, DHFFC point, get the TP, FP, and FN stats
    """
    input:
        config["outdir"]+"/{sample}/truvari-{dhffc}.txt"
    output:
        temp(config["outdir"]+"/{sample}/stats-{dhffc}.txt")
    run:
        with open(input[0], 'r') as truvari_report, open(output[0], 'w') as stats:
            for line in truvari_report:
                s = line.strip()
                if '"TP-call":' in s: TP = int(s.split()[1].rstrip(','))
                elif '"FP":' in s: FP = int(s.split()[1].rstrip(','))
                elif '"FN":' in s: FN = int(s.split()[1].rstrip(','))
            stats.write('\t'.join(wildcards.dhffc, TP, FP, FN))
            stats.write('\n')

rule combine_sample_stats:
    """
    Combine the sample, DHFFC stats into a single text file.
    Sorted by DHFFC
    """
    input:
        expand(config["outdir"]+"/{sample}/stats-{dhffc}.txt",
               dhffc=dhffc_range, sample=lambda w: [w.sample])
    output:
        config["outdir"]+"/{sample}/stats-combined.txt"
    shell:
        """
        cat <(printf "dhffc\tTP\tFP\tFN\n") {input} | sort > {output}
        """
        
