import argparse
import matplotlib.pyplot as plt
# from matplotlib_venn import venn3
import pysam

parser = argparse.ArgumentParser()
parser.add_argument('-a', dest='vcf_a', type=str, required=True,
                    help='VCF a')
# parser.add_argument('--a-name', dest='a_name', type=str, required=True,
#                     help='Set name of VCF a')

parser.add_argument('-b', dest='vcf_b', type=str, required=True)
# parser.add_argument('--b-name', dest='b_name', type=str, required=True,
#                     help='Set name of VCF b')

# parser.add_argument('-c', dest='vcf_c', type=str, required=True)
# parser.add_argument('--c-name', dest='c_name', type=str, required=True,
#                     help='Set name of VCF c')

args = parser.parse_args()

vcf_files = [args.vcf_a, args.vcf_b]
vcf_interval_sets = [set(), set()]

for filename, vcf_set in zip(vcf_files, vcf_interval_sets):
    with pysam.VariantFile(filename, 'r') as vcf:
        for variant in vcf:
            vcf_set.add(' '.join((str(variant.contig), 
                                  str(variant.pos), 
                                  str(variant.stop))))

print(len(set.intersection(*vcf_interval_sets)))

# venn3(vcf_interval_sets, set_labels=(args.a_name, 
#                                      args.b_name, 
#                                      args.c_name))
plt.show()
