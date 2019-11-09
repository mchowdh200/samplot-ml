import sys

gt_map = {
    '0': 'ref',
    '0|0': 'ref',
    '0|1': 'het',
    '1|0': 'het',
    '1|1': 'alt',
    '1': 'alt',
}

for line in sys.stdin:
    # format of a line is: genomic_position,sample\tGT,...
    A = line.rstrip().split(',')
    pos = A[0]
    sample_gt = [x.split() for x in A[3:]]

    for sample, gt in sample_gt:
        if '.' in gt: continue
        print(pos, sample, gt_map[gt], sep='\t')
