import sys

with open(sys.argv[1], 'w') as bed_file:
    for pred in sys.stdin:
        chr, start, end = pred.split('_')[:3]
        P = pred.split()[1:]
        print('\t'.join([chr, start, end] + P))
        bed_file.write('\t'.join([chr, start, end] + P) + '\n')
