import numpy as np
import pandas as pd

df = pd.read_csv('data/BED/del.sample.bed', sep='\t', 
                 names=['chr', 'start', 'end', 'sample', 'genotype']) 

sv_len = df.end - df.start
sv_len_lt1000 = sv_len[sv_len < 1000].values
sv_len_gt1000 = sv_len[sv_len > 1000].values

hist, bin_edges = np.histogram(sv_len_lt1000, bins=10)

print("SV length < 1000")
for i, n in enumerate(hist):
    print(f"{int(bin_edges[i])}-{int(bin_edges[i+1])} : {n}")


hist, bin_edges = np.histogram(sv_len_gt1000, bins=20)

print("SV length > 1000")
for i, n in enumerate(hist):
    print(f"{int(bin_edges[i])}-{int(bin_edges[i+1])} : {n}")
