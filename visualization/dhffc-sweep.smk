import os
from collections import Counter, defaultdict
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

callers = ["smoove", "manta"]
samples = ["HG002", "HG00514", "HG00733", "NA19240"]
dhffc_range = np.around(np.linspace(0, 1.0, 1001), 3)

rule all:
    input:
        config["outdir"]+"/dhffc-sweep.png",
        config["outdir"]+"/dhffc-sweep-baseline-diff.png",
        config["outdir"]+"/dhffc-sweep-baseline-diff-zoomed.png"

rule filter_dhffc:
    """
    Filter baseline genotyped SV vcf with DHFFC metric
    """
    input:
        lambda w: config["input_vcf"][w.caller][w.sample]
    output:
        vcf = temp(config["outdir"]+"/{caller}/{sample}/filtered-lt-{dhffc}.vcf.gz"),
        index = temp(config["outdir"]+"/{caller}/{sample}/filtered-lt-{dhffc}.vcf.gz.tbi")
    shell:
        """bcftools view -i 'DHFFC <= {wildcards.dhffc}' {input} |
           bgzip -c > {output.vcf}
           tabix {output.vcf}"""

rule eval_baseline:
    """
    get performance of non-filtered vcf for each sample.
    """
    input:
        truth_set = lambda w: config["truth_set"][w.sample],
        baseline_vcf = lambda w: config["input_vcf"][w.caller][w.sample]
    output:
        txt=config["outdir"]+"/{caller}/{sample}/baseline.txt",
        dir=temp(directory(config["outdir"]+"/{caller}/{sample}/baseline-truvari"))
    log:
        config["logdir"]+"/{caller}/{sample}/eval_baseline.log"
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
        config["outdir"]+"/{caller}/{sample}/baseline.txt"
    output:
        temp(config["outdir"]+"/{caller}/{sample}-{caller}-baseline-stats.txt")
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
        filtered = config["outdir"]+"/{caller}/{sample}/filtered-lt-{dhffc}.vcf.gz",
        filtered_index = config["outdir"]+"/{caller}/{sample}/filtered-lt-{dhffc}.vcf.gz.tbi",
        truth_set = lambda w: config["truth_set"][w.sample]
    output:
        txt=config["outdir"]+"/{caller}/{sample}/truvari-{dhffc}.txt",
        dir=temp(directory(config["outdir"]+"/{caller}/{sample}/truvari-{dhffc}"))
    log:
        config["logdir"]+"/{caller}/{sample}/evaluate-{dhffc}.log"
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
        config["outdir"]+"/{caller}/{sample}/truvari-{dhffc}.txt"
    output:
        temp(config["outdir"]+"/{caller}/{sample}/stats-{dhffc}.txt")
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
        expand(config["outdir"]+"/{caller}/{sample}/stats-{dhffc}.txt",
               dhffc=dhffc_range, sample="{sample}", caller="{caller}")
    output:
        config["outdir"]+"/{caller}/{sample}-{caller}-stats.txt"
    shell:
        """
        cat <(printf "dhffc\tTP\tFP\tFN\n") {input} > {output}
        """

rule plot_sample_stats:
    """
    Plot of FP vs TP for each sample on same plot.
    """
    input:
        stats = expand(config["outdir"]+"/{caller}/{sample}-{caller}-stats.txt",
                       caller=callers, sample=samples),
        baseline = expand(config["outdir"]+"/{caller}/{sample}-{caller}-baseline-stats.txt",
                          caller=callers, sample=samples)
    output:
        roc = config["outdir"]+"/dhffc-sweep.png",
        baseline_diff = config["outdir"]+"/dhffc-sweep-baseline-diff.png",
        baseline_diff_zoomed = config["outdir"]+"/dhffc-sweep-baseline-diff-zoomed.png"
    run:
        # get baseline stats
        baseline_stats = defaultdict(lambda: defaultdict(dict)) # this is gnarly but I love it.
        for file in input.baseline:
            sample = os.path.basename(file).split('-')[0]
            caller = os.path.basename(file).split('-')[1]
            with open(file) as f:
                (baseline_stats[sample][caller]['TP'],
                 baseline_stats[sample][caller]['FP'],
                 baseline_stats[sample][caller]['FN']) = \
                     map(int, f.readline().rstrip().split())

        # sum up baseline stats
        baseline_summed = Counter()
        for sample in samples:
            for caller in callers:
                for stat in ["TP", "FP", "FN"]:
                    baseline_summed[stat] += baseline_stats[sample][caller][stat]
                
        # get dhffc stats
        dhffc_stats = defaultdict(dict) # sample, caller -> dataframe
        for file in input.stats:
            sample = os.path.basename(file).split('-')[0]
            caller = os.path.basename(file).split('-')[1]

            dhffc_stats[sample][caller] = pd.read_csv(file, sep='\t') \
                                            .sort_values(by="dhffc")

        # sum up dhffc stats
        dhffc_summed = {
            stat: np.zeros_like(dhffc_stats[samples[0]][callers[0]][stat].values)
            for stat in ["TP", "FP", "FN"]}
        for sample in samples:
            for caller in callers:
                for stat in ["TP", "FP", "FN"]:
                    dhffc_summed[stat] += dhffc_stats[sample][caller][stat]
            
        # dhffc sweep roc
        for caller in callers:
            plt.gca().set_prop_cycle(None) # reset order of colors
            linestyle = "-" if caller == "smoove" else "--"
            for sample in samples:
                dhffc_stats[sample][caller]["FPR"] = \
                    dhffc_stats[sample][caller].apply(
                        lambda x: x.FP/(x.FP + x.FN), axis=1)

                dhffc_stats[sample][caller]["TPR"] = \
                    dhffc_stats[sample][caller].apply(
                        lambda x: x.TP/(x.TP + x.FN), axis=1)

                split_point = dhffc_stats[sample][caller].loc[
                    dhffc_stats[sample][caller]["dhffc"] == 0.7]

                plt.plot(dhffc_stats[sample][caller]["FPR"],
                         dhffc_stats[sample][caller]["TPR"],
                         linestyle=linestyle,
                         label=sample+" "+caller)
                plt.plot(split_point["FPR"], split_point["TPR"],
                         marker='o', color='k')

        
        plt.gca().spines['right'].set_visible(False)
        plt.gca().spines['top'].set_visible(False)
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positives Rate')
        plt.title('Duphold DHFFC sweep: Receiver operating characteristic')
        plt.legend()
        plt.savefig(output.roc)

        plt.gca().clear()

        # % change plot
        pct_diff_TP = np.abs(
            baseline_summed["TP"] - dhffc_summed["TP"])/baseline_summed["TP"] * 100
        pct_diff_FP = np.abs(
            baseline_summed["FP"] - dhffc_summed["FP"])/baseline_summed["FP"] * 100

        # to find the points 1.1% to 2.4% in TP subtract the TP vector by the target
        # and find the min afterwards
        a = np.argmin(np.abs(pct_diff_TP - 1.1))
        b = np.argmin(np.abs(pct_diff_TP - 2.4))

        plt.plot(pct_diff_FP, pct_diff_TP)
        plt.gca().invert_yaxis()
        plt.xlabel('% Reduction in False Positives')
        plt.ylabel('% Reduction in True Positives')
        plt.savefig(output.baseline_diff)
        plt.gca().clear()
        
        # zoomed % change plot
        plt.gca().invert_yaxis()
        plt.plot(pct_diff_FP, pct_diff_TP)
        plt.plot(pct_diff_FP[a], pct_diff_TP[a], marker='o', color='k')
        plt.plot(pct_diff_FP[b], pct_diff_TP[b], marker='o', color='k')
        plt.annotate(f"({pct_diff_FP[a]:.2f}, {pct_diff_TP[a]:.2f})",
                     (pct_diff_FP[a], pct_diff_TP[a]))
        plt.annotate(f"({pct_diff_FP[b]:.2f}, {pct_diff_TP[b]:.2f})",
                     (pct_diff_FP[b], pct_diff_TP[b]))
        plt.xlim([20, 40])
        plt.ylim([5, 0])
        plt.xlabel('% Reduction in False Positives')
        plt.ylabel('% Reduction in True Positives')

        plt.savefig(output.baseline_diff_zoomed)

                


