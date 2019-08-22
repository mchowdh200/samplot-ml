
"""
Read the input vcf and use the model predictions to filter vcf by some threshold
"""
import sys
import pysam

vcf_file=sys.argv[1]
bed_file=sys.argv[2]
threshold=sys.argv[3]

predictions = {}
genotypes = {0: (0, 0), 1: (0, 1), 2: (1, 1)}

# go through the predictions bed file
# and get the probability distribution
for l in open(bed_file, 'r'):
    A = l.rstrip().split()
    key = '\t'.join(A[:3])
    predictions[key] = [float(x) for x in A[3:]]

# iterate over each record and replace the genotype
# with the predicted genotypes
with pysam.VariantFile(vcf_file, 'rb') as VCF:
    print(str(VCF.header).strip())
    for variant in VCF:
        key = '\t'.join([str(x) for x in [variant.contig,
                                          variant.pos,
                                          variant.stop]]) 
        if key in predictions:
            if predictions[key][0] < threshold:
                print(str(variant).rstrip())
