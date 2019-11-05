import sys

for pred in sys.stdin:
    chr, start, end = pred.split('_')[:3]

    P = pred.split()[1:]
    print(chr, start, end, *P, sep='\t')
