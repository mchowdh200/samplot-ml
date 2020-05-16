#!/usr/bin/env python
import sys
import numpy as np
import matplotlib
import pylab
import random
from optparse import OptionParser
import matplotlib.gridspec as gridspec

from matplotlib import rcParams
rcParams['font.family'] = 'Arial'

delim = '\t'
parser = OptionParser()

parser.add_option(
    "-o",
    "--output_file",
    dest="output_file",
    help="Output file")

parser.add_option(
    "-i",
    "--input_file",
    dest="input_file",
    help="Input file")

parser.add_option(
    "--fig_x",
    dest="fig_x",
    type="int",
    default=5,
    help="Figure width")

parser.add_option(
    "--fig_y",
    dest="fig_y",
    type="int",
    default=5,
    help="Figure height")

(options, args) = parser.parse_args()
if not options.output_file:
    parser.error('Output file not given')

if not options.input_file:
    parser.error('Input file not given')

E = {}
header = None
for l in open(options.input_file):
    A = l.rstrip().split(',')
    if header is None:
        header = A
        continue

    sample = A[0]

    if sample not in E:
        E[sample] = {}

    caller = A[1]

    if caller not in  E[sample]:
        E[sample][caller] = {}

    filter_name = A[2]
    
    E[sample][caller][filter_name] = {
        'TP':int(A[3]),
        'FP':int(A[4]),
        'FN':int(A[5]),
        'Precision':float(A[6]),
        'Recall':float(A[7]),
        'F1':float(A[8]) }

fig = matplotlib.pyplot.figure(
    figsize=(options.fig_x,options.fig_y),
    dpi=300)

samples = ['HG002', 'HG00514', 'HG00733', 'NA19240']

outer_grid = gridspec.GridSpec(1, len(samples), wspace=0.2, hspace=0.0)


legend = False
sample_i = 0
for sample in samples:
    print(sample,sample_i)
    inner_grid = gridspec.GridSpecFromSubplotSpec(
        1, 2,
        subplot_spec=outer_grid[sample_i],
        wspace=0.05,
        hspace=0.0)

    ax_i = 0
    axs = []
    ys = []
    title_set = False
    for caller in ['smoove','manta']:
        ax = fig.add_subplot(inner_grid[ax_i])

        if not title_set:
            ax.set_title(sample,loc='left')
            title_set = True

        axs.append(ax)

        tps = []
        fps = []

        for filter_name in ['', 'DHFFC', 'CNN']:
            tps.append(E[sample][caller][filter_name]['TP'])
            fps.append(E[sample][caller][filter_name]['FP'])
        
        ys = tps + fps

        ind = np.arange(len(tps))
        width = 0.35
        rects1 = ax.bar(ind - width/2, tps, width, label='TP')
        rects2 = ax.bar(ind + width/2, fps, width, label='FP')

        ax.set_xlabel(caller)
        ax.set_xticks([0,1,2])
        ax.set_xticklabels(['No filter', 'DHFFC', 'CNN'], rotation=90)

        ax_i+=1

        if not legend:
            ax.legend(
                ncol=2,
                frameon=False,
                loc=1,
                bbox_to_anchor=(1.0, 1.2))
            legend = True

    for ax in axs:
        ax.set_ylim((0,max(ys)))
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
    for ax in axs[1:]:
        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.get_yaxis().set_visible(False)


    sample_i += 1

matplotlib.pyplot.savefig(options.output_file,bbox_inches='tight')

