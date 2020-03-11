import sys
import math
import random
import collections

"""
sample from the exclude regions bed and produce a bed format output
Note: regardless of the input genotype, we will consider the output
      genotype as 0/0
"""

class Region:
    def __init__(self, chrom, start, end, sample, AB, dhffc, GT):
        gt_map = {
            '0|0': 'ref',
            '0/0': 'ref',
            '0|1': 'het',
            '0/1': 'het',
            '1|0': 'het',
            '1/0': 'het',
            '1|1': 'alt',
            '1/1': 'alt',
        }
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)
        self.sample = sample
        self.AB = float(AB) if AB != '.' else AB
        self.dhffc = float(dhffc) if dhffc != '.' else dhffc
        self.GT = gt_map.get(GT, GT)

    def __str__(self):
        return '\t'.join([
            self.chrom, str(self.start), str(self.end), 
            self.sample, str(self.AB), str(self.dhffc), self.GT
        ])


if __name__ == '__main__':
    # number of entries to sample from each region
    sample_size = int(sys.argv[1])

    # % of regions to sample from (0-1)
    region_percentage = float(sys.argv[2])

    gt_override = None
    if len(sys.argv) == 4:
        gt_override = sys.argv[3]

    # organize entries by region ----------------------------------------------
    region_dict = collections.defaultdict(list)
    for line in sys.stdin:
        # input will be of the form:
        # chrA  startA endA sample AB dhffc GT   chrB startB endB
        # [0    1      2    3      4  5     6]  [7    8      9]
        # -------------- value --------------   ----- key -----
        line = line.rstrip().split()
        region_key = ''.join(line[7:])
        region_dict[region_key].append(
            Region(chrom=line[0],
                   start=line[1],
                   end=line[2],
                   sample=line[3],
                   AB=line[4],
                   dhffc=line[5],
                   GT=gt_override if gt_override else '0/0'))

    # randomly sample certain % of regions from the region_dict
    # note: list(dict) gives us a list of keys, so we are sampling keys
    region_list = random.sample(
        list(region_dict), math.floor(region_percentage * len(region_dict))
    )

    # randomly sample entries from each region --------------------------------
    for region in region_list:
        sampled_regions = random.sample(
            region_dict[region], min(len(region_dict[region]), sample_size))

        for s in sampled_regions:
            if s.dhffc == '.': continue
            print(s)


