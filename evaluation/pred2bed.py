import sys

for pred in sys.stdin:
    chr, start, end = pred.split('_')[:3]

    P = pred.split()[1:]
    print(chr, start, end, P[0], P[1], P[2], sep='\t')
