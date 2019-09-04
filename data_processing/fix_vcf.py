# import sys
import pysam

# with pysam.VariantFile(sys.argv[1], 'rb') as VCF:
with pysam.VariantFile("-", 'r') as VCF:
    print(str(VCF.header).strip())
    for variant in VCF:
        variant.stop = variant.pos - variant.info['SVLEN']
        print(str(variant).rstrip())
