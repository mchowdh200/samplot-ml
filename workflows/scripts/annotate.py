import sys
import numpy as np
import pysam

bed = sys.argv[1]
sample = sys.argv[2]

predictions = {}
genotypes = {0: (0, 0), 1: (0, 1), 2: (1, 1)}


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
    # get the old genotype and prediction scores into the vcf FORMAT fields
    VCF.header.add_meta('FORMAT', items=[('ID', 'OLDGT'),
                                         ('Number', '1'),
                                         ('Type', 'String'),
                                         ('Description', 'Genotype before samplot-ml')])
    VCF.header.add_meta('FORMAT', items=[('ID', 'PREF'),
                                         ('Number', '1'),
                                         ('Type', 'Float'),
                                         ('Description', 'Samplot-ml P(0/0) prediction score')])
    VCF.header.add_meta('FORMAT', items=[('ID', 'PHET'),
                                         ('Number', '1'),
                                         ('Type', 'Float'),
                                         ('Description', 'Samplot-ml P(0/1) prediction score')])
    VCF.header.add_meta('FORMAT', items=[('ID', 'PALT'),
                                         ('Number', '1'),
                                         ('Type', 'Float'),
                                         ('Description', 'Samplot-ml P(1/1) prediction score')])
    print(str(VCF.header).strip())
    for variant in VCF:
        region = '\t'.join([str(x) for x in [variant.contig,
                                          variant.pos,
                                          variant.stop]]) 
        oldgt = variant.samples[sample].allele_indices
        variant.samples[sample]["OLDGT"] = f"{oldgt[0]}/{oldgt[1]}"
        if region in predictions:
            pref, phet, palt = predictions[region]
            variant.samples[sample]["PREF"] = pref
            variant.samples[sample]["PHET"] = phet
            variant.samples[sample]["PALT"] = palt

            # set new GT
            variant.samples[sample].allele_indices = genotypes[
                np.argmax(predictions[region])]
        else:
            variant.samples[sample]["PREF"] = "."
            variant.samples[sample]["PHET"] = "."
            variant.samples[sample]["PALT"] = "."
        print(str(variant).rstrip())
