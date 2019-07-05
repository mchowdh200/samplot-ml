import sys
import cyvcf2

vcf_file=sys.argv[1]
bed_file=sys.argv[2]

ml = {}

# go through the predictions bed file
# and get the probability distribution
for l in open(bed_file, 'r'):
    A = l.rstrip().split()
    key = '\t'.join(A[:3])
    ml[key] = [float(x) for x in A[3:]]


vcf=cyvcf2.VCF(vcf_file)
vcf.add_info_to_header({'ID':'ml_ref','Number':1,'Type':'Float','Description':'P(ref)'})
vcf.add_info_to_header({'ID':'ml_het','Number':1,'Type':'Float','Description':'P(het)'})
vcf.add_info_to_header({'ID':'ml_alt','Number':1,'Type':'Float','Description':'P(alt)'})
vcf.add_info_to_header({'ID':'ml_argmax','Number':1,'Type':'Integer','Description':'argmax(P(gt))'})


print(vcf.raw_header.strip())

x = [0, 1, 2]

for variant in vcf: # or VCF('some.bcf')
    key = '\t'.join([str(x) for x in [variant.CHROM,
                                      variant.start + 1,
                                      variant.INFO.get('END')]])

    if key in ml:
        variant.INFO['ml_ref'] = ml[key][0]
        variant.INFO['ml_het'] = ml[key][1]
        variant.INFO['ml_alt'] = ml[key][2]

        variant.INFO['ml_argmax'] = max(x, key=lambda i: x[i])

        print(str(variant).rstrip())
