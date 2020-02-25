import sys
import random
import collections

# TODO factor this out since I'm using this in sample_exclude_regions.py
class Region:
    def __init__(self, chrom, start, end, sample, AB, dhffc, GT):
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)
        self.sample = sample
        self.AB = float(AB) if AB != '.' else AB
        self.dhffc = float(dhffc) if dhffc != '.' else dhffc
        self.GT = GT

        self.gt_map = {
            '0|0': 'ref',
            '0/0': 'ref',
            '0|1': 'het',
            '0/1': 'het',
            '1|0': 'het',
            '1/0': 'het',
            '1|1': 'alt',
            '1/1': 'alt',
        }
    def __str__(self):
        return '\t'.join([
            self.chrom, str(self.start), str(self.end), 
            self.sample, str(self.AB), str(self.dhffc), self.gt_map[self.GT]
        ])

if __name__ == '__main__':
    sample_size = int(sys.argv[1])

    # organize entries by region ----------------------------------------------
    region_dict = collections.defaultdict(list)
    for line in sys.stdin:
        # input format:
        # chr  start  end  sample  AB  DHFFC  GT
        # 0    1      2    3       4   5      6
        # ----- key -----  ------- value -------
        line = line.rstrip().split()
        region_key = ''.join(line[:3])

        region_dict[region_key].append(
            Region(chrom=line[0],
                   start=line[1],
                   end=line[2],
                   sample=line[3],
                   AB=line[4],
                   dhffc=line[5],
                   GT=line[6]))

    for region in region_dict:
        sampled_regions = random.sample(
            region_dict[region], min(len(region_dict[region]), sample_size))
        for s in sampled_regions:
            if s.dhffc == '.': continue
            print(s)

