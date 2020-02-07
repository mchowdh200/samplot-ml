"""
Read the input vcf and use the model predictions to annotate the vcf with the
prediction probabilities as well as change the called genotypes.
"""
import sys
import numpy as np
import pysam
# import cyvcf2

vcf_file=sys.argv[1]
bed_file=sys.argv[2]

predictions = {}
genotypes = {0: (0, 0), 1: (0, 1), 2: (1, 1)}

# go through the predictions bed file
# and get the probability distribution
for l in open(bed_file, 'r'):
    A = l.rstrip().split()
    key = '\t'.join(A[:3]) # region
    predictions[key] = [float(x) for x in A[3:]] # prediction score

# iterate over each record and replace the genotype
# with the predicted genotypes
with pysam.VariantFile(vcf_file, 'rb') as VCF:
    print(str(VCF.header).strip())
    for variant in VCF:
        key = '\t'.join([str(x) for x in [variant.contig,
                                          variant.pos,
                                          variant.stop]]) 
        if key in predictions:
            variant.samples[0].allele_indices = genotypes[
                np.argmax(predictions[key])]
            if np.argmax(predictions[key]) > 0:
                print(str(variant).rstrip())



# vcf = cyvcf2.VCF(vcf_file, gts012=True)
# vcf.add_info_to_header({'ID':'predictions_ref','Number':1,'Type':'Float','Description':'P(ref)'})
# vcf.add_info_to_header({'ID':'predictions_het','Number':1,'Type':'Float','Description':'P(het)'})
# vcf.add_info_to_header({'ID':'predictions_alt','Number':1,'Type':'Float','Description':'P(alt)'})
# vcf.add_info_to_header({'ID':'predictions_argmax','Number':1,'Type':'Integer','Description':'argmax(P(gt))'})
# print(vcf.raw_header.strip())

#for variant in vcf: # or VCF('some.bcf')
#    key = '\t'.join([str(x) for x in [variant.CHROM,
#                                      variant.start + 1, # .start is 0 based bed is 1 based
#                                      variant.INFO.get('END')]])
#    if key in predictions:
#        # variant.INFO['predictions_ref'] = predictions[key][0]
#        # variant.INFO['predictions_het'] = predictions[key][1]
#        # variant.INFO['predictions_alt'] = predictions[key][2]
#        # variant.INFO['predictions_argmax'] = max(x, key=lambda i: x[i])
#        # print(str(variant).rstrip())
