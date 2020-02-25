import sys
import numpy as np
import pandas as pd

"""
merge the duphold/svtyper annotated 1kg call sets.
* we want the original 1kg genotypes, but annotated with svtyper
  allele balance statistic.
"""

duphold_1kg_bed = sys.argv[1]
svtyper_1kg_bed = sys.argv[2]
out_file = sys.argv[3]

duphold_1kg_bed = pd.read_csv(
    duphold_1kg_bed, sep='\t',
    names=['chrom', 'start', 'end', 
           'sample', 'dhffc', 'gt_orig'],
    dtype=dict(chrom=str, start=int, end=int,
               sample=str, dhffc=str, gt_orig=str))

svtyper_1kg_bed = pd.read_csv(
    svtyper_1kg_bed, sep='\t',
    names=['chrom', 'start', 'end', 
           'sample', 'AB', 'dhffc', 'gt_new'],
    dtype=dict(chrom=str, start=int, end=int,
               sample=str, AB=str, dhffc=str, gt_new=str))

merged = pd.merge(duphold_1kg_bed,
                  svtyper_1kg_bed,
                  how='inner')

merged = merged[['chrom', 'start', 'end', 'sample', 'AB', 'dhffc', 'gt_orig']]
merged.to_csv(out_file, sep='\t', header=False, index=False)
