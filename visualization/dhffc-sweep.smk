import os
from collections import defaultdict
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

configfile: "dhffc-sweep-config.yaml"
input_vcf = config['input_vcf']
truth_set = config['truth_set']
outdir = config['outdir']
logdir = config['logdir']

samples = ["HG002", "HG00514", "HG00733", "NA19240"]
dhffc_range = np.around(np.linspace(0, 1.0, 101), 2)

rule all:
    input:
        config["outdir"]+"/dhffc-sweep.png"
        # expand(config["outdir"]+"/{sample}/stats-{dhffc}.txt",
        #        dhffc=dhffc_range, sample=samples)
        # expand(config["outdir"]+"/{sample}-stats.txt",
        #        sample=samples)

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

rule eval_baseline:
    """
    get performance of non-filtered vcf for each sample.
    """
    input:
        truth_set = lambda w: config["truth_set"][w.sample]
        baseline_vcf = lambda w: config["input_vcf"][w.sample]
    output:
        txt=config["outdir"]+"/{sample}/baseline.txt",
        dir=temp(directory(config["outdir"]+"/{sample}/baseline-truvari"))
    log:
        config["logdir"]+"/{sample}/eval_baseline.log"
    shell:
        """bash ../evaluation/truvari.sh -b {input.truth_set} \
                           -c {input.baseline_vcf} \
                           -o {output.dir} 2>&1 |
           tee {log} > {output.txt}"""

rule get_baseline_stats:
    """
    get TP, FP, FN info from baseline report
    """
    input:
        config["outdir"]+"/{sample}/baseline.txt"
    output:
        temp(config["outdir"]+"/{sample}-baseline-stats.txt")
    run:
        # TODO put in function
        with open(input[0], 'r') as truvari_report, open(output[0], 'w') as stats:
            for line in truvari_report:
                s = line.strip()
                if '"TP-call":' in s: TP = int(s.split()[1].rstrip(','))
                elif '"FP":' in s: FP = int(s.split()[1].rstrip(','))
                elif '"FN":' in s: FN = int(s.split()[1].rstrip(','))
            stats.write('\t'.join(map(str, [TP, FP, FN])))
            stats.write('\n')
        

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
            stats.write('\t'.join(map(str, [wildcards.dhffc, TP, FP, FN])))
            stats.write('\n')

rule combine_sample_stats:
    """
    Combine the sample, DHFFC stats into a single text file.
    Sorted by DHFFC
    """
    input:
        expand(config["outdir"]+"/{sample}/stats-{dhffc}.txt",
               dhffc=dhffc_range, sample="{sample}")
    output:
        config["outdir"]+"/{sample}-stats.txt"
    shell:
        """
        cat <(printf "dhffc\tTP\tFP\tFN\n") {input} > {output}
        """

rule plot_sample_stats:
    """
    Plot of FP vs TP for each sample on same plot.
    """
    input:
        stats = expand(config["outdir"]+"/{sample}-stats.txt",
                       sample=samples)
        baseline = expand(config["outdir"]+"/{sample}-baseline-stats.txt",
                          sample=samples)
    output:
        roc = config["outdir"]+"/dhffc-sweep.png"
        diff = config["outdir"]+"/dhffc-sweep-diff.png"
        baseline_diff = config["outdir"]+"/dhffc-sweep-baseline-diff.png"
    run:
        # get baseline stats
        baseline_stats = defaultdict(dict)
        for file in input.baseline:
            sample = os.path.basename(file).split('-')[0]
            with open(file) as f:
                (baseline_stats[sample]['TP'],
                 baseline_stats[sample]['FP'],
                 baseline_stats[sample]['FN']) = \
                     map(int, f.readline().rsrtip().split())
                
        # construct plots
        dhffc_sweep = plt.figure()
        dhffc_diff = plt.figure()
        baseline_diff = plt.figure()
        for file in input.stats:
            sample = os.path.basename(file).split('-')[0]
            df = pd.read_csv(file, sep='\t') \
                   .sort_values(by="dhffc")

            # dhffc sweep plot
            df['TPR'] = df.apply(lambda x: x.TP/(x.TP + x.FN), axis=1)
            df['FPR'] = df.apply(lambda x: x.FP/(x.FP + x.FN), axis=1)
            split_point = df.loc[df["dhffc"] == 0.7]

            plt.plot(df.FPR, df.TPR, label=sample, figure=dhffc_sweep)
            plt.plot(split_point.FPR, split_point.TPR,
                     marker='o', color='k', figure=dhffc_sweep)
            plt.xlabel('False Positive Rate', figure=dhffc_sweep)
            plt.ylabel('True Positives Rate', figure=dhffc_sweep)
            plt.title('Duphold DHFFC sweep: Receiver operating characteristic',
                      figure=dhffc_sweep)
            plt.legend(figure=dhffc_sweep)
            plt.savefig(output.roc, figure=dhffc_sweep)

            # dhffc sweep diff
            plt.plot(df.FP.diff(), df.TP.diff(), label=sample, figure=dhffc_diff)

            # baseline diff


