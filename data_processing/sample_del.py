
import sys
import numpy as np

for l in sys.stdin:
    A = l.rstrip().split()

    chrm,start,end = A[0:3]

    samples = [x.split(',') for x in A[3:]]

    ref = []
    het = []
    alt = []

    for sample in samples:
        name,gt = sample

        if gt == '0|0':
            ref.append(name)
        elif gt == '1|0' or gt == '0|1':
            het.append(name)
        elif gt == '1|1':
            alt.append(name)

    het_get = min(5, len(het))
    ref_get = min(5, len(ref))
    alt_get = min(5, len(alt))


    het_is = []
    ref_is = []
    alt_is = []

    if het_get > 0:
        het_is = np.random.choice(len(het), het_get)
    if ref_get > 0:
        ref_is = np.random.choice(len(ref), ref_get)
    if alt_get > 0:
        alt_is = np.random.choice(len(alt), alt_get)

    O = [chrm,start,end,0,0]

    for s in ref_is:
        O[3] = ref[s]
        O[4] = 'ref'
        print('\t'.join(O))

    for s in het_is:
        O[3] = het[s]
        O[4] = 'het'
        print('\t'.join(O))

    for s in alt_is:
        O[3] = alt[s]
        O[4] = 'alt'
        print('\t'.join(O))
