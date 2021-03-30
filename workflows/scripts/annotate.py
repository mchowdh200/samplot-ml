import sys
import numpy as np
import pysam

bed = sys.argv[1]
sample = sys.argv[2]

predictions = {}
genotypes = {0: (0, 0), 1: (0, 1), 2: (1, 1)}

# TODO create a new format field for OLD_GT
# TODO create a new INFO field for the prediction scores


# go through the predictions bed file
# and get the predictions keyed by region
for l in open(bed, 'r'):
    A = l.rstrip().split()
    region = '\t'.join(A[:3])
    predictions[region] = [float(x) for x in A[3:]] # prediction score

# Pipe the vcf from stdin.
# iterate over each record and replace the genotype
# with the predicted genotypes
with pysam.VariantFile('-', 'rb') as VCF:
    print(str(VCF.header).strip())
    for variant in VCF:
        key = '\t'.join([str(x) for x in [variant.contig,
                                          variant.pos,
                                          variant.stop]]) 
        if key in predictions:
            variant.samples[sample].allele_indices = genotypes[
                np.argmax(predictions[key])]
            if np.argmax(predictions[key]) > 0:
                print(str(variant).rstrip())
